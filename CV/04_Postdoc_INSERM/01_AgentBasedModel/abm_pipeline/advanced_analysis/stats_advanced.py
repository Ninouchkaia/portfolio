# abm_pipeline/advanced_analysis/stats_advanced.py

from __future__ import annotations
import pandas as pd
from pathlib import Path
from scipy.stats import mannwhitneyu

from abm_pipeline.advanced_analysis.utils import (
    load_pareto_dataframe,
    ensure_dir,
)
from abm_pipeline.config import ADVANCED_RESULTS_DIR
from abm_pipeline.parameter_exploration.utils import logger

def run_parameter_stats_tests(
    pareto_file_1: str,
    pareto_file_2: str,
    save_dir: str | None = None,
    exclude_cols=None,
):
    """
    Test statistique entre deux jeux de paramètres optimisés
    (ex: classe 1 vs classe 2)
    Equivalent de stat_test_paretos.py.

    Sauvegarde résultats dans TSV.
    """
    if save_dir is None:
        save_path = ensure_dir(ADVANCED_RESULTS_DIR / "stats")
    else:
        save_path = ensure_dir(save_dir)

    df1 = load_pareto_dataframe(pareto_file_1)
    df2 = load_pareto_dataframe(pareto_file_2)

    if exclude_cols is None:
        exclude_cols = ["delta_fitness_via", "delta_fitness_conc"]

    params = [c for c in df1.columns if c not in exclude_cols]

    results = []
    for param in params:
        stat, pval = mannwhitneyu(df1[param], df2[param], alternative="two-sided")
        results.append([param, stat, pval])

    out_df = pd.DataFrame(results, columns=["parameter", "MannWhitneyU", "pvalue"])
    out_df.to_csv(save_path / "stats_comparison.tsv", sep="\t", index=False)

    logger.info(f"Statistical comparison written: {save_path / 'stats_comparison.tsv'}")
