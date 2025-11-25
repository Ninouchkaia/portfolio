# paxdb/src/fasta_parser.py

from __future__ import annotations

from pathlib import Path
from typing import Dict, Iterator, Tuple


def parse_fasta(fasta_path: Path) -> Iterator[Tuple[str, str]]:
    """
    Simple FASTA parser that yields (sequence_id, sequence) pairs.

    This avoids external dependencies (no Biopython needed) and is
    sufficient for proteome-scale amino acid counting.

    The sequence identifier is taken as the first token after '>'.

    Parameters
    ----------
    fasta_path : Path
        Path to the FASTA file.

    Yields
    ------
    (seq_id, sequence) : (str, str)
        Sequence identifier and amino acid sequence (no whitespace).
    """
    seq_id = None
    seq_chunks = []

    with fasta_path.open("r") as f:
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                # Emit previous sequence if any
                if seq_id is not None:
                    yield seq_id, "".join(seq_chunks)
                # New sequence
                header = line[1:].strip()
                seq_id = header.split()[0]
                seq_chunks = []
            else:
                seq_chunks.append(line.strip())

    # Emit last record
    if seq_id is not None:
        yield seq_id, "".join(seq_chunks)


def load_fasta_as_dict(fasta_path: Path) -> Dict[str, str]:
    """
    Convenience wrapper to load a FASTA into a dict: seq_id -> sequence.
    """
    return dict(parse_fasta(fasta_path))
