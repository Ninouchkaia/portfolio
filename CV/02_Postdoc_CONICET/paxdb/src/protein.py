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
        unweighted counts.
    """

    seq_id: str
    sequence: str
    abundance: float = 1.0

    def length(self) -> int:
        """Return the sequence length."""
        return len(self.sequence)

    def amino_acid_counts(self) -> Dict[str, int]:
        """
        Return a raw count of each amino acid in this protein.

        Ambiguous residues (e.g. B, Z, X) are counted separately and can
        optionally be discarded in downstream analyses.
        """
        return dict(Counter(self.sequence))

    def amino_acid_frequencies(self) -> Dict[str, float]:
        """
        Return fractional composition (per protein).

        Frequencies sum to 1.0 over all characters present in the sequence.
        """
        counts = self.amino_acid_counts()
        total = float(sum(counts.values())) or 1.0
        return {aa: c / total for aa, c in counts.items()}

    def weighted_amino_acid_counts(self) -> Dict[str, float]:
        """
        Return abundance-weighted amino acid counts.

        Each raw count is multiplied by the protein abundance.
        This is the key link between PaxDB abundances and proteome-level
        amino acid usage.
        """
        counts = self.amino_acid_counts()
        return {aa: c * self.abundance for aa, c in counts.items()}
