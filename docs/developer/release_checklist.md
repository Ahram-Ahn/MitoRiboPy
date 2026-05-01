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

## 9. Reproducibility / archival

* [ ] Zenodo archive minted for the tag (or explicitly deferred with
  a comment in the release notes — do not ship a manuscript pinning
  a version that is not archived)
* [ ] DOI / citation block in `README.md` updated, or the reason for
  not updating is in the release notes

## 10. Post-release sanity

```bash
python -m pip install --upgrade 'mitoribopy>=<X.Y.Z>'
mitoribopy --version
```

* [ ] `pip install --upgrade 'mitoribopy>=<X.Y.Z>'` from PyPI works
  on a clean env
* [ ] one quick toy run succeeds against the freshly-installed copy
* [ ] CHANGELOG entry's "Unreleased" section moved to the new version

---

**Do not submit a manuscript that pins this version until every box
above is ticked.** The PyPI / GitHub / manuscript triplet must agree.
