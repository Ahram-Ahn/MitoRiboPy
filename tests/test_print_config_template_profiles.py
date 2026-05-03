"""Tests for ``mitoribopy all --print-config-template --profile X``.

The v0.7.0 contract:

* ``--profile minimal`` (default) — the curated short template; backwards-
  compatible with pre-v0.7.0 callers that did not pass a profile.
* ``--profile publication`` — minimal + a publication-readiness overlay
  declaring `strict: true`, the periodicity statistical-hardening
  defaults, and a commented `rnaseq_mode: de_table` block.
* ``--profile exhaustive`` — the full annotated example from
  ``examples/templates/pipeline_config.example.yaml`` (every flag with
  its default and a one-line comment).
"""

from __future__ import annotations

import io
import sys
from contextlib import redirect_stdout, redirect_stderr
from pathlib import Path

import pytest

from mitoribopy.cli.all_ import (
    _CONFIG_TEMPLATE,
    _exhaustive_template_path,
    _print_config_template,
)


def _capture(profile: str | None = None) -> tuple[str, str]:
    out = io.StringIO()
    err = io.StringIO()
    with redirect_stdout(out), redirect_stderr(err):
        if profile is None:
            _print_config_template()
        else:
            _print_config_template(profile=profile)
    return out.getvalue(), err.getvalue()


# ---------------------------------------------------------------------------
# Default + minimal
# ---------------------------------------------------------------------------


def test_default_profile_is_minimal_and_backwards_compatible() -> None:
    """Calling _print_config_template() with no arg must keep the
    pre-v0.7.0 behaviour: emit the minimal template verbatim."""
    out_default, err_default = _capture()
    out_minimal, err_minimal = _capture("minimal")
    assert out_default == out_minimal
    assert out_default == _CONFIG_TEMPLATE
    assert err_default == ""
    assert err_minimal == ""


def test_minimal_profile_includes_every_stage_section() -> None:
    out, _ = _capture("minimal")
    assert "align:" in out
    assert "rpf:" in out
    assert "rnaseq:" in out  # commented but present
    # Should NOT carry the publication overlay markers.
    assert "publication-readiness profile" not in out


# ---------------------------------------------------------------------------
# Publication overlay
# ---------------------------------------------------------------------------


def test_publication_profile_includes_minimal_then_overlay() -> None:
    out, err = _capture("publication")
    assert err == ""
    # Minimal template appears first in full.
    assert _CONFIG_TEMPLATE in out
    # Publication overlay markers.
    assert "publication-readiness profile" in out
    assert "strict: true" in out
    # Statistical-hardening knobs surfaced explicitly.
    assert "fourier_bootstrap_n: 200" in out
    assert "fourier_permutations_n: 200" in out
    assert "metagene_normalize: per_gene_unit_mean" in out
    # Recommends the de_table flow.
    assert "rnaseq_mode: de_table" in out


def test_publication_overlay_is_appended_not_prepended() -> None:
    out, _ = _capture("publication")
    minimal_idx = out.index("# ---- align ")
    overlay_idx = out.index("publication-readiness profile")
    assert minimal_idx < overlay_idx


# ---------------------------------------------------------------------------
# Exhaustive template
# ---------------------------------------------------------------------------


def test_exhaustive_template_path_resolves_in_repo_layout() -> None:
    path = _exhaustive_template_path()
    # The repo-layout copy should always exist when running tests from
    # the repo. Pip-installed checkouts can return None — guarded
    # separately in test_exhaustive_profile_falls_back_when_template_missing.
    assert path is not None
    assert path.is_file()
    assert path.name == "pipeline_config.example.yaml"


def test_exhaustive_profile_emits_every_documented_section() -> None:
    out, err = _capture("exhaustive")
    assert err == ""
    # The exhaustive template carries the long header comment about
    # listing every flag.
    assert "exhaustive pipeline_config example" in out
    # All three stage sections are present.
    assert "align:" in out
    assert "rpf:" in out
    assert "rnaseq" in out
    # And several keys that are NOT in the minimal template (so we
    # know we got the long version).
    assert "max_parallel_samples" in out or "max_threads_per_sample" in out


def test_unknown_profile_falls_back_to_minimal_with_warning() -> None:
    out, err = _capture("not_a_real_profile")
    assert out == _CONFIG_TEMPLATE
    assert "unknown --profile" in err
    assert "not_a_real_profile" in err
