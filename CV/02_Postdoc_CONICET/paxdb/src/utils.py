# paxdb/src/utils.py

from __future__ import annotations

import logging
from pathlib import Path
from typing import Optional


def setup_logging(log_dir: Path, verbose: bool = True) -> None:
    """
    Configure a basic logging setup (both console and file).

    Parameters
    ----------
    log_dir : Path
        Directory where the log file will be written.
    verbose : bool
        If True, log INFO to console; otherwise only WARNING+.
    """
    log_dir.mkdir(parents=True, exist_ok=True)
    log_file = log_dir / "paxdb_analysis.log"

    # Root logger
    logger = logging.getLogger()
    logger.setLevel(logging.DEBUG)

    # Console handler
    ch = logging.StreamHandler()
    ch.setLevel(logging.INFO if verbose else logging.WARNING)
    ch_formatter = logging.Formatter("[%(levelname)s] %(message)s")
    ch.setFormatter(ch_formatter)

    # File handler
    fh = logging.FileHandler(log_file, mode="w")
    fh.setLevel(logging.DEBUG)
    fh_formatter = logging.Formatter(
        "%(asctime)s [%(levelname)s] %(name)s - %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )
    fh.setFormatter(fh_formatter)

    logger.handlers = []
    logger.addHandler(ch)
    logger.addHandler(fh)


def resolve_path(p):
    from pathlib import Path
    p = Path(p)
    return p.expanduser().resolve()
