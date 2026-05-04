# Release checklist

Mandatory gate for every PyPI cut. Skipping any step is a release-day
incident; if a step cannot be done now, postpone the release rather
than ship.

The fastest way to use this is to copy the boxes into a release issue
and tick them as you go.

---

## 0. Pick the version

* Patch (`0.6.x`) for bug fixes / docs / CI hygiene.
* Minor (`0.x.0`) for additive features (new flags, new subcommands).
* Major (`x.0.0`) for breaking changes (removed flags, manifest
  schema with renamed fields).

> v0.6.2 (2026-05-01) is a **patch** release: fourth-edit cleanup +
> reproducibility scaffolding, no biological logic change.

## 1. Version-source agreement

Every place the version is recorded MUST read the same string.

```bash
grep -n 'version = ' pyproject.toml
grep -n '__version__' src/mitoribopy/__init__.py
grep -n 'Version:'    src/MitoRiboPy.egg-info/PKG-INFO   # rebuilt in step 5
grep -n 'Manuscript-target version' README.md
grep -n '^## \['                       CHANGELOG.md | head -3
```

* [ ] `pyproject.toml` `[project] version`
* [ ] `src/mitoribopy/__init__.py` `__version__` fallback
* [ ] `README.md` "Manuscript-target version" + `pip install` lower bound
* [ ] `CHANGELOG.md` has a dated entry for the new version
* [ ] `docs/release-notes/v<X.Y.Z>.md` exists
* [ ] CLI quick-start examples mention the right version

```bash
python -m pip install -e .
mitoribopy --version          # MUST print the new version
```

## 2. Template + config validation

Every shipped YAML template parses with `validate-config --strict`.

```bash
for f in examples/templates/*.yaml; do
  echo "VALIDATING $f"
  mitoribopy validate-config "$f" --strict || exit 1
done
```

* [ ] every `examples/templates/*.yaml` passes `validate-config --strict`
* [ ] no public template uses `strain: h`, `strain: y`, `merge_density`,
  `fastq_dir`, or pre-canonical kit names as the recommended form
* [ ] `dedup_strategy:` examples use `auto` or `umi_coordinate`, not
  `umi-tools` (the legacy alias still parses; we just don't lead with it)
* [ ] `allow_pseudo_replicates_for_demo_not_publication:` uses the
  long canonical key — the short form must round-trip via
  `migrate-config` to the long form

## 3. Tests + linting

```bash
python -m pytest -q
python -m pytest -q -m requires_tools   # only when external tools are on PATH
```

* [ ] full unit-test suite green on the local machine
* [ ] CI green on Python 3.10 / 3.11 / 3.12
* [ ] no new `FutureWarning`/`DeprecationWarning` from numpy / pandas /
  matplotlib in the captured stderr

## 4. Toy data smoke test

The toy bundle under `examples/toy_data/` (when present) MUST run end
to end on a laptop in under ~10 minutes.

```bash
mitoribopy all \
  --config examples/toy_data/pipeline_config.yaml \
  --output /tmp/mitoribopy_toy_release_test \
  --strict
```

* [ ] `run_manifest.json` written, `schema_version` ≥ 1.3.0
* [ ] `canonical_config.yaml` written
* [ ] `resource_plan.json` written
* [ ] `warnings.tsv` written (header-only is OK)
* [ ] `outputs_index.tsv` written
* [ ] `SUMMARY.md` written
* [ ] `align/umi_qc.tsv` written
* [ ] `rpf/rpf_counts.tsv` matches `examples/toy_data/expected/expected_rpf_counts.tsv` (when published)

## 5. Build artifacts

```bash
rm -rf dist/ build/ src/MitoRiboPy.egg-info
python -m pip install --upgrade build twine
python -m build
twine check dist/*
```

* [ ] `dist/MitoRiboPy-<X.Y.Z>-py3-none-any.whl` builds
* [ ] `dist/MitoRiboPy-<X.Y.Z>.tar.gz` builds
* [ ] `twine check dist/*` passes

## 6. Clean-environment install (TestPyPI)

```bash
conda create -n mitoribopy_release_test python=3.12 -y
conda activate mitoribopy_release_test
python -m pip install dist/MitoRiboPy-<X.Y.Z>-py3-none-any.whl
mitoribopy --version
mitoribopy --help
mitoribopy rpf --help
mitoribopy all --help
mitoribopy all --print-config-template > /tmp/template.yaml
mitoribopy validate-config /tmp/template.yaml --strict --no-path-checks
```

* [ ] wheel installs cleanly into a fresh conda env
* [ ] `mitoribopy --version` matches the release tag
* [ ] every subcommand's `--help` renders without traceback
* [ ] `--print-config-template > template.yaml; validate-config --strict
  --no-path-checks template.yaml` passes

```bash
twine upload --repository testpypi dist/*
python -m pip install --index-url https://test.pypi.org/simple/ \
                     --extra-index-url https://pypi.org/simple/ \
                     "mitoribopy==<X.Y.Z>"
```

* [ ] uploaded to TestPyPI
* [ ] installs from TestPyPI in a clean env

## 7. Git tag + GitHub release

```bash
git tag -a v<X.Y.Z> -m "MitoRiboPy v<X.Y.Z>"
git push origin v<X.Y.Z>
```

* [ ] release branch / `main` is at the exact commit you tested
* [ ] tag pushed
* [ ] GitHub release created and links to the relevant
  `docs/release-notes/v<X.Y.Z>.md`

## 8. PyPI upload

```bash
twine upload dist/*
```

* [ ] uploaded to PyPI
* [ ] PyPI project page lists the new version as latest
* [ ] PyPI long description matches the README "Manuscript compatibility"
  table for the new version
* [ ] PyPI classifier reflects `Development Status :: 4 - Beta` (or
  later) — do NOT regress to `Alpha`

## 9. Documentation link check

`README.md` and the docs under `docs/` link out heavily — the
release-day check is that none of those links rot before the version
is the latest on PyPI.

```bash
# Internal links: every relative .md / .yaml / .tsv / .png / .svg
# / .csv link target must exist somewhere in the tree.
python - <<'PY'
import re, sys
from pathlib import Path
ROOT = Path.cwd()
LINK = re.compile(
    r"\[[^\]]+\]\("
    r"([^)#]+\.(?:md|svg|png|yaml|yml|tsv|csv|py|json))"
    r"\)"
)
broken: list[str] = []
for md in ROOT.rglob("*.md"):
    if any(p in md.parts for p in ("dist", "build", ".git", ".pytest_cache")):
        continue
    text = md.read_text(encoding="utf-8")
    for target in LINK.findall(text):
        if target.startswith(("http://", "https://", "mailto:")):
            continue
        candidate = (md.parent / target).resolve()
        if not candidate.exists():
            broken.append(f"{md.relative_to(ROOT)} -> {target}")
sys.exit(
    "Broken internal links:\n  " + "\n  ".join(broken)
    if broken else 0
)
PY
```

* [ ] internal Markdown links resolve on a clean checkout (every
  relative `*.md`, image, and example file referenced from any doc
  exists at the path the link claims)
* [ ] the README's CHANGELOG / docs / examples links land on the
  same commit as the release tag (no `main`-pinned URLs that will
  drift after the next merge)
* [ ] `docs/release-notes/v<X.Y.Z>.md` links to the merged
  `CHANGELOG.md` entry, not a draft

## 10. Reproducibility / archival

* [ ] Zenodo archive minted for the tag (or explicitly deferred with
  a comment in the release notes — do not ship a manuscript pinning
  a version that is not archived)
* [ ] DOI / citation block in `README.md` updated, or the reason for
  not updating is in the release notes
* [ ] `CITATION.cff` `version` and `date-released` fields match the
  release tag

## 11. Post-release sanity

```bash
python -m pip install --upgrade 'mitoribopy>=<X.Y.Z>'
mitoribopy --version
```

* [ ] `pip install --upgrade 'mitoribopy>=<X.Y.Z>'` from PyPI works
  on a clean env
* [ ] one quick toy run succeeds against the freshly-installed copy
* [ ] CHANGELOG entry's "Unreleased" section moved to the new version

## 12. Manuscript version pin

The PyPI release exists so a paper's Methods can pin a single,
installable, archived version. Before the manuscript ships:

* [ ] the Methods section names the exact PyPI version
  (`mitoribopy==<X.Y.Z>`), not `>=` and not the GitHub `main` branch
* [ ] every command in the Methods reproduces in a clean conda env
  with that exact pin (the [docs/rnaseq_te.md](../rnaseq_te.md)
  publication-boundary gates are honoured: `--strict`,
  `rnaseq_mode: de_table` for any TE / ΔTE claim)
* [ ] the Methods cites the Zenodo DOI minted in step 10
* [ ] supplementary materials include the exact YAML configs used
  (these validate cleanly under
  `mitoribopy validate-config --strict`)
* [ ] supplementary materials include `run_manifest.json`,
  `canonical_config.yaml`, and `outputs_index.tsv` from the
  publication run (the auditable provenance triplet)

---

**Do not submit a manuscript that pins this version until every box
above is ticked.** The PyPI / GitHub / manuscript triplet must agree.
