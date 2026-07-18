#!/usr/bin/env python3
"""Compatibility launcher for the packaged ``group-ops`` command.

New code should import :mod:`group_theory_operations` or use ``group-ops``.
"""

from __future__ import annotations

from pathlib import Path
import sys


SOURCE = Path(__file__).resolve().parent / "src"
if SOURCE.is_dir() and str(SOURCE) not in sys.path:
    sys.path.insert(0, str(SOURCE))

from group_theory_operations import *  # noqa: F401,F403,E402
from group_theory_operations.cli import main  # noqa: E402


if __name__ == "__main__":
    raise SystemExit(main())
