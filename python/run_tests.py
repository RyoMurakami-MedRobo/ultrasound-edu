#!/usr/bin/env python3
"""Run the Python-port test suite without needing pytest installed.

Prefers pytest when available; otherwise falls back to the module runners.
"""

import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent


def main() -> int:
    try:
        import pytest  # noqa: F401
        return subprocess.call([sys.executable, "-m", "pytest", "-q", str(HERE / "tests")])
    except ImportError:
        rc = 0
        for name in ("test_contract.py", "test_parity.py", "test_trace.py"):
            print(f"=== {name} ===")
            rc |= subprocess.call([sys.executable, str(HERE / "tests" / name)])
        return rc


if __name__ == "__main__":
    raise SystemExit(main())
