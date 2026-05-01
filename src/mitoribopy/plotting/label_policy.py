"""Shared label-density policy for MitoRiboPy plots.

The policy follows assessment §5: never label every point on a
whole-transcriptome volcano, but always label every point on an
mt-mRNA-only plot. The policy is deterministic, depends only on the
visible point count + a per-point ranking score, and can be inspected
by the figure validator (the chosen policy name is recorded in the
plot's ``.metadata.json`` sidecar).

Buckets (default thresholds):

* ``n_points <= 20``           → label every point (``policy="all"``)
* ``20 < n_points <= 200``     → label the top ``small_top`` (default 15)
                                  by absolute ranking score
                                  (``policy="top_n_small"``)
* ``n_points > 200``           → label the top ``large_top`` (default 10)
                                  by absolute ranking score and write
                                  the full ranked list into
                                  ``label_candidates.tsv`` so a reviewer
                                  can audit the cut
                                  (``policy="top_n_large"``)

The ranking score is supplied by the caller (e.g. ``-log10(padj)`` for
volcanoes, ``residual_abs`` for codon-correlation scatters). When no
score is supplied, the policy falls back to plot-order; that's fine
for small mt-only plots where every point is labelled anyway.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Sequence


__all__ = [
    "LabelDecision",
    "decide_labels",
    "write_label_candidates",
]


@dataclass(frozen=True)
class LabelDecision:
    """The outcome of :func:`decide_labels`.

    * ``policy`` is the string recorded in the plot sidecar.
    * ``indices`` are positions into the original ``points`` list that
      should be labelled. The order matches input order so the caller
      can iterate ``points`` once.
    * ``n_points`` is the total visible point count (== ``len(points)``).
    """

    policy: str
    indices: tuple[int, ...]
    n_points: int

    @property
    def n_labelled(self) -> int:
        return len(self.indices)


def decide_labels(
    points: Sequence,
    *,
    scores: Sequence[float] | None = None,
    annotate: str | bool = "auto",
    small_threshold: int = 20,
    medium_threshold: int = 200,
    small_top: int = 15,
    large_top: int = 10,
) -> LabelDecision:
    """Decide which point indices to label according to the auto policy.

    ``annotate`` accepts:

    * ``"auto"`` (default) — apply the bucket policy above.
    * ``"all"`` or ``True`` — label every point regardless of count.
    * ``"none"`` or ``False`` — label nothing.
    * an integer ``N`` — label the top ``N`` by score (or the first ``N``
      when no score is supplied).

    ``scores`` should be the same length as ``points``; each entry is
    interpreted as a ranking magnitude (the higher the more important).
    Pass already-signed ``-log10(padj)`` for volcanoes, ``abs(residual)``
    for codon-correlation scatters, etc. ``None`` is treated as 0 and
    therefore deprioritised.
    """
    n = len(points)
    if annotate is True or annotate == "all":
        return LabelDecision("all", tuple(range(n)), n)
    if annotate is False or annotate == "none":
        return LabelDecision("none", (), n)
    if isinstance(annotate, int) and not isinstance(annotate, bool):
        keep = max(0, int(annotate))
        return LabelDecision(
            f"top_{keep}",
            _top_indices(scores, keep, n),
            n,
        )

    # Bucketed auto policy.
    if n <= small_threshold:
        return LabelDecision("all", tuple(range(n)), n)
    if n <= medium_threshold:
        return LabelDecision(
            "top_n_small",
            _top_indices(scores, small_top, n),
            n,
        )
    return LabelDecision(
        "top_n_large",
        _top_indices(scores, large_top, n),
        n,
    )


def _top_indices(
    scores: Sequence[float] | None, k: int, n: int
) -> tuple[int, ...]:
    """Return the indices of the top-*k* entries in *scores*.

    Falls back to the first ``k`` indices when ``scores`` is None / too
    short. Ties break by original index (stable). Output is sorted by
    original index so the caller can keep one pass over points.
    """
    if k <= 0 or n == 0:
        return ()
    if scores is None or len(scores) < n:
        return tuple(range(min(k, n)))

    def _key(i: int) -> tuple:
        s = scores[i]
        rank = -float("inf") if s is None else float(s)
        # Negative for descending sort.
        return (-rank, i)

    ordered = sorted(range(n), key=_key)
    chosen = ordered[: min(k, n)]
    return tuple(sorted(chosen))


def write_label_candidates(
    path: Path | str,
    *,
    items: Sequence[tuple[str, float | None]],
) -> Path:
    """Write a TSV of ``label\\tscore`` rows in descending score order.

    Useful for the ``policy="top_n_large"`` case where the plot only
    shows the top ``large_top`` labels but a reviewer wants to audit
    the full ranked list.
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    # None scores sort AFTER any real number (lowest priority); among
    # real numbers, descending magnitude wins.
    ordered = sorted(
        items,
        key=lambda item: (item[1] is None, -float(item[1]) if item[1] is not None else 0.0),
    )
    with path.open("w", encoding="utf-8") as handle:
        handle.write("label\tscore\n")
        for label, score in ordered:
            score_str = "" if score is None else f"{float(score):.6g}"
            handle.write(f"{label}\t{score_str}\n")
    return path
