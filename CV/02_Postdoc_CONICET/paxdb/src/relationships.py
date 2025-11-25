# paxdb/src/relationships.py

from __future__ import annotations

import logging
from dataclasses import dataclass
from typing import Dict, Iterable, List, Tuple

import numpy as np
import pandas as pd
from scipy.stats import pearsonr

from paxdb.src.aa_metrics import (
    STANDARD_AA,
    aa_vector_from_costs,
    aa_vector_from_freqs,
)

logger = logging.getLogger(__name__)


@dataclass
class AARelationshipResult:
    """
    Store correlation between amino-acid usage and biosynthetic cost
    for a given species.
    """

    species_id: str
    pearson_r: float
    p_value: float
    used_aas: List[str]

    def as_dict(self) -> Dict[str, object]:
        return {
            "species_id": self.species_id,
            "pearson_r": self.pearson_r,
            "p_value": self.p_value,
            "n_amino_acids": len(self.used_aas),
        }


def correlate_usage_with_cost(
    species_id: str,
    aa_usage: pd.DataFrame,
    aa_costs: Dict[str, float],
) -> AARelationshipResult:
    """
    Compute Pearson correlation between AA usage and AA cost for a species.

    Parameters
    ----------
    species_id : str
        Identifier of the species (e.g. 4932-S.cerevisiae).
    aa_usage : DataFrame
        Table with index=AA, columns including 'freq' or 'weighted_freq'.
    aa_costs : dict
        Mapping AA -> cost.

    Returns
    -------
    AARelationshipResult
    """
    # Align usage and cost vectors over the STANDARD_AA list
    cost_vec = aa_vector_from_costs(aa_costs, STANDARD_AA)
    usage_vec = aa_vector_from_freqs(aa_usage, column="freq", amino_acids=STANDARD_AA)

    # Filter out NaNs
    mask = ~np.isnan(cost_vec) & ~np.isnan(usage_vec)
    used_costs = cost_vec[mask]
    used_usage = usage_vec[mask]
    used_aas = [aa for aa, m in zip(STANDARD_AA, mask) if m]

    if len(used_aas) < 3:
        logger.warning(
            "Species %s: only %d amino acids usable for correlation; returning NaNs.",
            species_id,
            len(used_aas),
        )
        return AARelationshipResult(species_id, float("nan"), float("nan"), used_aas)

    r, p = pearsonr(used_usage, used_costs)

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
