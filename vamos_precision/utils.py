from __future__ import annotations
import logging
import subprocess
from pathlib import Path

def setup_logging(verbosity: int = 1):
    level = logging.WARNING if verbosity <= 0 else logging.INFO if verbosity == 1 else logging.DEBUG
    logging.basicConfig(level=level, format="%(levelname)s:%(name)s:%(message)s")

def ensure_dir(path: str | Path) -> Path:
    p = Path(path)
    p.mkdir(parents=True, exist_ok=True)
    return p

def run_r_script(script_path: str | Path, args: list[str] | None = None) -> int:
    """Run an R script using Rscript."""
    cmd = ["Rscript", str(script_path)] + (args or [])
    return subprocess.call(cmd)
