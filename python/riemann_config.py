"""Minimal parser for the euler1d ``.cfg`` format (``key = value``, ``#`` comments)."""
from __future__ import annotations

from pathlib import Path

_INT_KEYS = {"n_cells", "max_steps", "wavenumber"}


def load_config(path: str | Path) -> dict[str, float | int | str]:
    """Parse a case file. Numeric values become float (or int for count-like keys);
    other values stay strings. ``gamma`` defaults to 1.4 when absent."""
    cfg: dict[str, float | int | str] = {}
    for lineno, raw in enumerate(Path(path).read_text().splitlines(), start=1):
        line = raw.split("#", 1)[0].strip()
        if not line:
            continue
        if "=" not in line:
            raise ValueError(f"{path}:{lineno}: expected 'key = value'")
        key, value = (part.strip() for part in line.split("=", 1))
        if key in _INT_KEYS:
            cfg[key] = int(value)
        else:
            try:
                cfg[key] = float(value)
            except ValueError:
                cfg[key] = value
    cfg.setdefault("gamma", 1.4)
    return cfg
