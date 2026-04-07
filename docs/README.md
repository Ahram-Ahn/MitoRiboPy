# Documentation

This directory collects package documentation that is useful after installation,
for development, and for release tracking.

## Start here

- [../README.md](../README.md): package overview, installation, CLI usage, and input/output formats
- [release-notes/README.md](release-notes/README.md): version-by-version release notes
- [PACKAGE_REFACTOR_SCHEME.md](PACKAGE_REFACTOR_SCHEME.md): historical refactor design notes with current-status annotations

## Architecture

- [diagrams/mitoribopy_package_flow.mmd](diagrams/mitoribopy_package_flow.mmd): authoritative Mermaid source for the current package flow
- [diagrams/render_mitoribopy_diagram.sh](diagrams/render_mitoribopy_diagram.sh): helper script to render the Mermaid source to SVG

## Historical reports

The files under [reports](reports) are archived migration, validation, and audit
notes from the packaging effort. They are kept for scientific traceability and
project history. Some of those reports intentionally use historical module names
such as `inframe_analysis`, `igv`, or `VARNA` because they describe older
development stages.
