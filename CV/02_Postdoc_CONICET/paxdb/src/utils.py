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

    # Clear existing handlers (useful if re-running in the same process)
    if logger.handlers:
        for h in list(logger.handlers):
            logger.removeHandler(h)

    # File handler (always DEBUG)
    fh = logging.FileHandler(str(log_file), mode="w", encoding="utf-8")
    fh.setLevel(logging.DEBUG)
    fh_formatter = logging.Formatter(
        fmt="%(asctime)s [%(levelname)s] %(name)s - %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )
    fh.setFormatter(fh_formatter)
    logger.addHandler(fh)

    # Console handler
    ch = logging.StreamHandler()
    ch.setLevel(logging.INFO if verbose else logging.WARNING)
    ch_formatter = logging.Formatter("%(levelname)s - %(message)s")
    ch.setFormatter(ch_formatter)
    logger.addHandler(ch)

    logging.getLogger(__name__).info("Logging initialized. Log file: %s", log_file)


def resolve_path(base: Path, relative: str) -> Path:
    """
    Resolve a path relative to some base directory.

    This is helpful when the species metadata gives paths relative
    to its own location.

    Parameters
    ----------
    base : Path
        Base directory (e.g. the directory containing species.tsv).
    relative : str
        Relative path string from that base.

    Returns
    -------
    Path
        Absolute path.
    """
    p = Path(relative)
    if not p.is_absolute():
        p = base / p
    return p.resolve()
