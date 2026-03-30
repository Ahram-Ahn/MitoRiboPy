from __future__ import annotations

import json

import pytest

from mitoribopy.config.runtime import load_user_config, resolve_rpf_range


def test_resolve_rpf_range_from_cli_values() -> None:
    assert resolve_rpf_range("h", [29, 34]) == [29, 30, 31, 32, 33, 34]


def test_resolve_rpf_range_defaults_by_strain() -> None:
    assert resolve_rpf_range("y", None) == [37, 38, 39, 40, 41]
    assert resolve_rpf_range("h", None) == [28, 29, 30, 31, 32, 33, 34]


def test_resolve_rpf_range_rejects_invalid_interval() -> None:
    with pytest.raises(ValueError):
        resolve_rpf_range("h", [34, 29])


def test_load_user_config_keeps_known_keys_and_ignores_unknown(tmp_path, capsys) -> None:
    config_path = tmp_path / "config.json"
    config_path.write_text(
        json.dumps({"strain": "h", "align": "stop", "unknown_key": 123}),
        encoding="utf-8",
    )

    cfg = load_user_config(str(config_path))

    assert cfg == {"strain": "h", "align": "stop"}
    captured = capsys.readouterr()
    assert "Ignoring unknown keys: unknown_key" in captured.out
