"""IGV-compatible BedGraph exports of footprint density.

Reads the per-transcript footprint_density CSVs written by
``run_translation_profile_analysis`` and emits one BedGraph per sample
per requested site. The output directory layout is::

    igv_tracks/
        <sample>/
            <sample>_p_site.bedgraph
            <sample>_a_site.bedgraph

Each file starts with a ``track`` header that gives IGV a default name
and color (P-site forest green, A-site dark orange). Consecutive
positions with equal values are collapsed into a single span, which is
the canonical BedGraph compression and keeps the files small even at
nucleotide resolution.

BedGraph was chosen over BigWig to avoid the pyBigWig binary
dependency. mt-genome scale is small enough that the readability gain
of plain text outweighs the size penalty.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from ..console import log_info, log_warning


_SITE_LABEL = {"p": "P-site", "a": "A-site"}
_SITE_COLOR = {"p": "34,139,34", "a": "255,140,0"}
_SITE_COLUMN = {"p": "P_site", "a": "A_site"}
_SITE_FILENAME = {"p": "p_site", "a": "a_site"}


def _emit_bedgraph_for_sample(
    sample_dir: Path,
    site: str,
    out_path: Path,
) -> int:
    """Write one BedGraph for ``sample`` covering every transcript.

    Returns the number of non-zero spans written. Returns 0 when the
    sample's footprint_density folder is missing or empty so the caller
    can decide whether to log a skip.
    """
    foot_dir = sample_dir / "footprint_density"
    if not foot_dir.is_dir():
        return 0

    sample_name = sample_dir.name
    site_column = _SITE_COLUMN[site]
    spans_written = 0

    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", encoding="utf-8") as out:
        out.write(
            f'track type=bedGraph name="{sample_name}_{_SITE_FILENAME[site]}" '
            f'description="MitoRiboPy footprint density ({_SITE_LABEL[site]})" '
            f'color={_SITE_COLOR[site]} visibility=full\n'
        )
        for csv_path in sorted(foot_dir.glob("*_footprint_density.csv")):
            transcript = csv_path.name[: -len("_footprint_density.csv")]
            df = pd.read_csv(csv_path)
            if site_column not in df.columns or "Position" not in df.columns:
                continue
            # CSV ``Position`` column is 1-based (Position 1 = nt 0
            # in BED-style 0-based coords). Convert to BED 0-based
            # half-open spans and collapse runs of equal value.
            positions = df["Position"].to_numpy()
            values = df[site_column].to_numpy()
            n = len(values)
            i = 0
            while i < n:
                v = values[i]
                if v <= 0:
                    i += 1
                    continue
                j = i + 1
                while j < n and values[j] == v:
                    j += 1
                start_bed = int(positions[i]) - 1
                end_bed = int(positions[j - 1])
                out.write(f"{transcript}\t{start_bed}\t{end_bed}\t{v}\n")
                spans_written += 1
                i = j
    return spans_written


def run_igv_export(
    translation_profile_dir: str | Path,
    output_dir: str | Path,
    sites: list[str] | None = None,
) -> None:
    """Walk ``translation_profile_dir`` and write a BedGraph per sample.

    Parameters
    ----------
    translation_profile_dir
        Flat ``translation_profile/`` root containing one subdir per
        sample, each with a ``footprint_density/`` folder.
    output_dir
        Where to write ``<sample>/<sample>_<site>.bedgraph``.
    sites
        Sites to export. Defaults to ``["p", "a"]``.
    """
    translation_profile_dir = Path(translation_profile_dir)
    output_dir = Path(output_dir)
    if sites is None:
        sites = ["p", "a"]
    sites = [s for s in (str(s).lower() for s in sites) if s in _SITE_COLUMN]
    if not sites:
        log_warning("IGV", "No valid sites requested; skipping IGV export.")
        return
    if not translation_profile_dir.is_dir():
        log_warning(
            "IGV",
            f"Translation-profile dir not found => {translation_profile_dir}; skipping.",
        )
        return

    sample_dirs = sorted(p for p in translation_profile_dir.iterdir() if p.is_dir())
    if not sample_dirs:
        log_warning("IGV", "No sample subdirs under translation_profile/; skipping.")
        return

    for sample_dir in sample_dirs:
        sample_name = sample_dir.name
        sample_out = output_dir / sample_name
        for site in sites:
            out_path = sample_out / f"{sample_name}_{_SITE_FILENAME[site]}.bedgraph"
            spans = _emit_bedgraph_for_sample(sample_dir, site, out_path)
            log_info(
                "IGV",
                f"{sample_name}/{_SITE_FILENAME[site]}: wrote {spans} spans => {out_path}",
            )

    log_info("IGV", f"Wrote BedGraph tracks under {output_dir}.")
