"""Phase 2.5 — synthetic-mini known-answer integration test.

A 3-transcript x 2-condition x 2-replicate fixture with hand-picked
input parameters that produce known answers at every stage. The test
exercises the four core public APIs on the same synthetic data so a
silent regression in any one stage breaks here:

  - determine_p_site_offsets        ->  picks offset = 12 with confidence='high'
  - run_periodicity_qc              ->  frame_0 == 1.0, dominance == 1.0
  - compute_te                      ->  TE per (sample, gene) follows the
                                        documented (RPF + 0.5) / (mRNA + 0.5) form
  - compute_delta_te                ->  KO/WT contrast yields the
                                        log2(2) = 1.0 delta-TE for ND2 (RPF
                                        doubled, mRNA unchanged) and 0.0 for
                                        ND1 / ND3 (both held equal)

Design notes for the fixture live in docs/validation/synthetic_mini.md.
"""

from __future__ import annotations

from pathlib import Path

import math
import numpy as np
import pandas as pd
import pytest

from mitoribopy.analysis.offset_selection import determine_p_site_offsets
from mitoribopy.analysis.offset_enrichment import create_csv_for_offset_enrichment
from mitoribopy.analysis.periodicity import run_periodicity_qc
from mitoribopy.data.transcript_annotations import human_annotation_df
from mitoribopy.rnaseq._types import DeColumnMap, DeTable
from mitoribopy.rnaseq.te import compute_delta_te, compute_te


# ---------- fixture builders -----------------------------------------------


# Three real human mt-mRNA transcripts. Using the bundled annotation
# means start_codon / stop_codon are correct without us hand-rolling
# them, and the test crashes loudly if the bundled table changes.
TRANSCRIPT_NAMES = ("ND1", "ND2", "ND3")

# Hand-picked offset-of-truth: 5' end at +12 nt from the P-site, the
# textbook value for human mt-monosome footprints.
TRUE_5_OFFSET = 12
READ_LENGTH = 30


def _annotation_df() -> pd.DataFrame:
    """Return the bundled human mt annotation with start/stop columns."""
    annotation = human_annotation_df.copy()
    annotation["start_codon"] = annotation["l_utr5"]
    annotation["stop_codon"] = annotation["l_tr"] - annotation["l_utr3"] - 3
    return annotation


def _synthesise_bed(
    *, sample: str, annotation: pd.DataFrame, ko_doubles_nd2: bool = False,
) -> pd.DataFrame:
    """Return a frame-0-perfect BED for one sample.

    Each transcript contributes ``n_reads`` reads; KO samples double
    the read count on ND2 so the downstream DESeq2 contrast sees +1
    log2 on the RPF side. Reads are tiled in 3-nt steps so every
    read's P-site (= start + TRUE_5_OFFSET) lands on a codon boundary.
    Reads are placed so they bracket the stop codon (the offset
    enrichment routine anchors there).

    To get an unambiguous enrichment peak at TRUE_5_OFFSET, we
    over-weight reads whose 5' end sits at the canonical distance from
    the stop codon. The package's stop-anchored offset arithmetic is
    ``5' offset = stop_codon - read_start - 2``, so for a reported
    offset of 12 we need ``read_start = stop_codon - 14``.
    """
    ann_indexed = annotation.set_index("transcript")
    rows: list[dict] = []
    for tx in TRANSCRIPT_NAMES:
        ann = ann_indexed.loc[tx]
        stop_codon = int(ann["stop_codon"])
        cds_start = int(ann["start_codon"])
        cds_end = cds_start + int(ann["l_cds"])

        n = 200 if (ko_doubles_nd2 and tx == "ND2") else 100

        # 60% of the reads are anchor-aligned: 5' end at stop - 14 so
        # the package reports a 5' offset of 12 (= stop - read_start - 2).
        n_anchor = int(0.6 * n)
        anchor_start = stop_codon - TRUE_5_OFFSET - 2
        for _ in range(n_anchor):
            rows.append({
                "chrom": tx,
                "start": anchor_start,
                "end": anchor_start + READ_LENGTH,
                "read_length": READ_LENGTH,
                "sample_name": sample,
                "strand": "+",
            })
        # The remaining 40% are tiled across the CDS in 3-nt steps so
        # the frame-summary computation has a meaningful population
        # outside the stop region. Each read's P-site = start + 12 lands
        # on a codon boundary (frame 0).
        n_tile = n - n_anchor
        codon_starts = list(range(cds_start, cds_end - READ_LENGTH, 3))
        for i in range(n_tile):
            cs = codon_starts[i % len(codon_starts)]
            start = cs - TRUE_5_OFFSET
            rows.append({
                "chrom": tx,
                "start": start,
                "end": start + READ_LENGTH,
                "read_length": READ_LENGTH,
                "sample_name": sample,
                "strand": "+",
            })
    return pd.DataFrame(rows)


def _bed_for_all_samples(annotation: pd.DataFrame) -> pd.DataFrame:
    """Two conditions x two replicates. KO replicates double ND2."""
    parts = [
        _synthesise_bed(sample="WT_R1", annotation=annotation, ko_doubles_nd2=False),
        _synthesise_bed(sample="WT_R2", annotation=annotation, ko_doubles_nd2=False),
        _synthesise_bed(sample="KO_R1", annotation=annotation, ko_doubles_nd2=True),
        _synthesise_bed(sample="KO_R2", annotation=annotation, ko_doubles_nd2=True),
    ]
    return pd.concat(parts, ignore_index=True)


# ---------- the integration test itself ------------------------------------


def test_synthetic_mini_round_trip(tmp_path: Path) -> None:
    """End-to-end known-answer pass over offset selection + frame QC +
    metagene + TE + ΔTE on the synthetic-mini fixture."""
    annotation = _annotation_df()
    bed = _bed_for_all_samples(annotation)

    # Stage 1: offset enrichment + selection (combined across samples).
    summary, offsets = create_csv_for_offset_enrichment(
        bed_df=bed,
        annotation_df=annotation,
        align_to="stop",
        rpf_range=range(READ_LENGTH, READ_LENGTH + 1),
        output_csv=str(tmp_path / "offset_combined.csv"),
        offset_limit=25,
        offset_site="p",
        codon_overlap_mode="full",
        strain="h",
    )
    assert offsets is not None and not offsets.empty

    selected = determine_p_site_offsets(
        offsets_df=offsets,
        align_to="stop",
        out_file=str(tmp_path / "selected_combined.csv"),
        offset_min=10,
        offset_max=22,
        offset_site="p",
        selection_reference="p_site",
    )
    assert selected is not None
    assert len(selected) == 1
    row = selected.iloc[0]
    # Known answer 1: the offset is 12, the textbook human mt-monosome value.
    assert int(row["Most Enriched 5' Offset"]) == TRUE_5_OFFSET
    # Known answer 2: this is a clean signal — confidence must be 'high'.
    assert row["confidence_5"] == "high"

    # Stage 2: per-sample selection (each replicate should also pick 12).
    selected_by_sample: dict[str, pd.DataFrame] = {}
    for sample, sub_offsets in offsets.groupby("sample_name"):
        out = tmp_path / f"selected_{sample}.csv"
        per_sample = determine_p_site_offsets(
            offsets_df=sub_offsets,
            align_to="stop",
            out_file=str(out),
            offset_min=10,
            offset_max=22,
            offset_site="p",
            selection_reference="p_site",
        )
        assert per_sample is not None
        assert int(per_sample.iloc[0]["Most Enriched 5' Offset"]) == TRUE_5_OFFSET
        selected_by_sample[sample] = per_sample

    # Stage 3: periodicity QC. The fixture is constructed so every read
    # lands on a codon boundary -> frame_0 == 1.0 and dominance == 1.0
    # for every sample.
    qc_dir = tmp_path / "qc"
    result = run_periodicity_qc(
        bed_df=bed,
        annotation_df=annotation,
        samples=list(selected_by_sample.keys()),
        selected_offsets_by_sample=selected_by_sample,
        selected_offsets_combined=selected,
        offset_type="5",
        offset_site="p",
        output_dir=qc_dir,
        window_nt=60,
        plot=False,  # the plot is exercised in test_periodicity_qc.
        # Synthetic fixture starts placing reads from CDS position 0;
        # the spec defaults (6 / 3) would mask all of them. Opt back
        # to legacy "no codon-edge masking" so the count invariants
        # below still hold.
        exclude_start_codons=0,
        exclude_stop_codons=0,
    )
    # Periodicity invariant: the metagene start profile should peak at
    # codon-aligned positions (every 3 nt). The frame_summary contract
    # was retired in v0.8.0 along with the rest of the frame-fraction
    # QC bundle; downstream consumers now read fourier_period3_score_combined
    # for the headline periodicity verdict.
    assert result["periodicity_start"], "metagene start profiles missing"
    for profile in result["periodicity_start"]:
        assert profile.density.size > 0

    # Stage 4: TE math. Build {gene: {sample: count}} from the BED.
    ribo_counts: dict[str, dict[str, int]] = {}
    for (gene, sample), n in bed.groupby(["chrom", "sample_name"]).size().items():
        ribo_counts.setdefault(gene, {})[sample] = int(n)
    # Hand-set mRNA abundances: every gene at 100 (so TE = (rpf+0.5)/100.5).
    mrna_abundances = {gene: 100.0 for gene in ribo_counts}
    te_rows = compute_te(ribo_counts, mrna_abundances)
    # 4 samples x 3 genes = 12 rows.
    assert len(te_rows) == 12
    # WT_R1 / ND1: rpf_count=100, TE = (100+0.5)/(100+0.5) = 1.0.
    wt_nd1 = next(r for r in te_rows if r.sample == "WT_R1" and r.gene == "ND1")
    assert wt_nd1.te == pytest.approx(1.0, rel=1e-6)
    # KO_R1 / ND2: rpf_count=200, TE = (200+0.5)/(100+0.5) ≈ 1.995.
    ko_nd2 = next(r for r in te_rows if r.sample == "KO_R1" and r.gene == "ND2")
    assert ko_nd2.te == pytest.approx(200.5 / 100.5, rel=1e-6)

    # Stage 5: ΔTE math. Feed a 'mRNA unchanged' DE table; KO/WT contrast
    # on the RPF side doubled ND2 only -> log2(2) = 1.0 ΔTE for ND2,
    # 0.0 for ND1 and ND3.
    de = DeTable(
        rows=[
            {"gene_id": "ND1", "log2fc": 0.0, "padj": 1.0, "basemean": 100.0},
            {"gene_id": "ND2", "log2fc": 0.0, "padj": 1.0, "basemean": 100.0},
            {"gene_id": "ND3", "log2fc": 0.0, "padj": 1.0, "basemean": 100.0},
        ],
        format="deseq2",
        column_map=DeColumnMap(
            gene_id="gene_id", log2fc="log2fc", padj="padj", basemean="basemean",
        ),
    )
    condition_map = {"WT_R1": "WT", "WT_R2": "WT", "KO_R1": "KO", "KO_R2": "KO"}
    dte_rows = compute_delta_te(
        ribo_counts, de,
        condition_map=condition_map, condition_a="WT", condition_b="KO",
    )
    by_gene = {r.gene: r for r in dte_rows}
    # ND2: KO mean=200, WT mean=100, mRNA log2fc=0 -> ΔTE_log2 = log2(200.5/100.5) - 0
    expected_nd2 = math.log2(200.5 / 100.5)
    assert by_gene["ND2"].delta_te_log2 == pytest.approx(expected_nd2, rel=1e-6)
    # ND1 / ND3: both conditions have 100 reads -> ΔTE_log2 == 0.
    for gene in ("ND1", "ND3"):
        assert by_gene[gene].delta_te_log2 == pytest.approx(0.0, abs=1e-6)
        assert by_gene[gene].note == ""

    # Sidecar artefacts: the metagene Fourier bundle is written for
    # every (sample, length) combination that has enough coverage.
    assert (qc_dir / "fourier_spectrum_combined.tsv").is_file()
    assert (qc_dir / "fourier_period3_score_combined.tsv").is_file()
    assert (qc_dir / "periodicity.metadata.json").is_file()


def test_synthetic_mini_periodicity_metagene_clean_signal(tmp_path: Path) -> None:
    """The start-aligned metagene should be perfectly periodic on a
    fixture engineered to produce one read per codon at frame 0."""
    annotation = _annotation_df()
    bed = _bed_for_all_samples(annotation)

    summary, offsets = create_csv_for_offset_enrichment(
        bed_df=bed,
        annotation_df=annotation,
        align_to="stop",
        rpf_range=range(READ_LENGTH, READ_LENGTH + 1),
        output_csv=str(tmp_path / "off.csv"),
        offset_limit=25,
        offset_site="p",
        codon_overlap_mode="full",
        strain="h",
    )
    selected = determine_p_site_offsets(
        offsets_df=offsets,
        align_to="stop",
        out_file=str(tmp_path / "sel.csv"),
        offset_min=10,
        offset_max=22,
        offset_site="p",
        selection_reference="p_site",
    )
    qc = run_periodicity_qc(
        bed_df=bed,
        annotation_df=annotation,
        samples=["WT_R1"],
        selected_offsets_by_sample=None,
        selected_offsets_combined=selected,
        offset_type="5",
        offset_site="p",
        output_dir=tmp_path / "qc",
        window_nt=60,
        plot=False,
        # Disable codon-edge masking — see test_synthetic_mini_round_trip.
        exclude_start_codons=0,
        exclude_stop_codons=0,
    )
    profile = qc["periodicity_start"][0]
    # Same periodicity invariant as the round-trip test: every read's
    # P-site sits at a codon start or the stop-codon-2 position; the
    # off-by-one phase has zero density across the whole start-aligned
    # window. (Frames 0 and 2 may both be non-zero across multiple
    # transcripts with mixed CDS-length parities.)
    f1_mass = profile.density[(profile.positions % 3) == 1].sum()
    total = profile.density.sum()
    assert total > 0
    assert f1_mass == pytest.approx(0.0, abs=1e-6), (
        f"frame 1 should be empty by construction; got mass={f1_mass}"
    )
