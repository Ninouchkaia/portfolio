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
        Sequence identifier and sequence string (no whitespace).
    """
    with fasta_path.open("r", encoding="utf-8") as handle:
        seq_id = None
        chunks = []

        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                # Finish previous sequence
                if seq_id is not None:
                    yield seq_id, "".join(chunks)
                # Start new sequence
                header = line[1:]
                seq_id = header.split()[0]
                chunks = []
            else:
                chunks.append(line)

        # Last sequence
        if seq_id is not None:
            yield seq_id, "".join(chunks)


def load_fasta_as_dict(fasta_path: Path) -> Dict[str, str]:
    """
    Load a FASTA file into a simple dict {seq_id: sequence}.

    Parameters
    ----------
    fasta_path : Path

    Returns
    -------
    dict
        Mapping from sequence ID to protein sequence.
    """
    return {seq_id: seq for seq_id, seq in parse_fasta(fasta_path)}
