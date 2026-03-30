from mitoribopy.cli import _normalize_args


def test_normalize_args_strips_run_and_separator() -> None:
    assert _normalize_args(["run", "--", "-s", "h"]) == ["-s", "h"]


def test_normalize_args_leaves_standard_argv() -> None:
    assert _normalize_args(["-s", "y", "--align", "start"]) == ["-s", "y", "--align", "start"]
