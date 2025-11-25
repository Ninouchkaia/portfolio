# paxdb/src/protein.py

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict
from collections import Counter


@dataclass
class Protein:
    """
    Representation of a protein with sequence and abundance.

    This is a minimal, explicit version of what the original
    newdef_protein* scripts were doing.

    Attributes
    ----------
    seq_id : str
        Sequence identifier (e.g. UniProt id).
    sequence : str
        Amino acid sequence (one-letter codes).
    abundance : float
        Relative abundance (e.g. PaxDB score). Defaults to 1.0 for
        unweighted analyses.
    """

    seq_id: str
    sequence: str
    abundance: float = 1.0

    def aa_counts(self) -> Dict[str, int]:
        """
        Return a dictionary mapping amino acid -> count in this sequence.
        """
        return dict(Counter(self.sequence))
