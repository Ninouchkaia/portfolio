# paxdb/src/relationships.py

from __future__ import annotations

import logging
from dataclasses import dataclass
from typing import Dict, Iterable, List, Tuple

import numpy as np
import pandas as pd
from scipy.stats import pearsonr

from .aa_metrics import (
    STANDARD_AA,
    aa_vector_from_costs,
    aa_vector_from_freqs,
)

logger = logging.getLogger(__name__)


@dataclass
class AARelationshipResult:
    """
    Container for the relationship between amino-acid usage
    and amino-acid cost for a given species.
    """

    species_id: str
    r: float
    pvalue: float
    used_aas: List[str]

    def as_dict(self) -> Dict[str, object]:
        return {
            "species_id": self.species_id,
            "pearson_r": self.r,
            "pvalue": self.pvalue,
            "n_amino_acids": len(self.used_aas),
            "amino_acids": "".join(self.used_aas),
        }


def correlate_usage_with_cost(
    species_id: str,
    freqs: Dict[str, float],
    costs: Dict[str, float],
    amino_acids: Iterable[str] = STANDARD_AA,
    drop_nan_costs: bool = True,
) -> AARelationshipResult:
    """
    Compute Pearson correlation between amino-acid usage and cost.

    This corresponds conceptually to parts of the MBE paper where
    proteome-level frequencies are compared to energetic costs.

    Parameters
    ----------
    species_id : str
        Species identifier for reporting.
    freqs : dict
        AA -> relative frequency.
    costs : dict
        AA -> energetic cost (e.g. ATP units).
    amino_acids : iterable of str
        Amino acids to consider (default: 20 standard).
    drop_nan_costs : bool
        If True, ignore amino acids that have missing cost values.

    Returns
    -------
    AARelationshipResult
        Pearson r, p-value, and which AAs were used.
    """
    aa_list = list(amino_acids)

    x = aa_vector_from_freqs(freqs, aa_list)
    y = aa_vector_from_costs(costs, aa_list)

    if drop_nan_costs:
        mask = ~np.isnan(y)
        x = x[mask]
        y = y[mask]
        used_aas = [aa for aa, keep in zip(aa_list, mask) if keep]
    else:
        used_aas = aa_list

    if x.size < 3:
        raise ValueError("Not enough amino acids with defined cost for correlation.")

    r, p = pearsonr(x, y)
    logger.info(
        "Species %s: Pearson r(usage, cost) = %.3f (p=%.2e) over %d amino acids",
        species_id,
        r,
        p,
        len(used_aas),
    )
    return AARelationshipResult(species_id, float(r), float(p), used_aas)


def results_to_dataframe(results: Iterable[AARelationshipResult]) -> pd.DataFrame:
    """
    Convert a list of AARelationshipResult to a tidy pandas DataFrame.
    """
    rows = [res.as_dict() for res in results]
    return pd.DataFrame(rows)
