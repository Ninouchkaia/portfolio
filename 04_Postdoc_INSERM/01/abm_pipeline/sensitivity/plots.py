# abm_pipeline/sensitivity/plots.py

from __future__ import annotations
from pathlib import Path
from typing import List, Dict

import matplotlib.pyplot as plt
import pandas as pd

from abm_pipeline.sensitivity.utils import load_sensitivity_csv
from abm_pipeline.parameter_exploration.utils import logger


def plot_sensitivity_for_param(
    exp_name: str,
    csv_dir: str,
    save_dir: str,
) -> None:
    """
    Produit une figure (viabilité + concentration) pour UN paramètre perturbé.

    exp_name est du type 'perturb-gui-apo-mov',
    CSV attendu : {csv_dir}/ABM_2D_sensitivity_{exp_name}.csv
    """
    csv_path = Path(csv_dir) / f"ABM_2D_sensitivity_{exp_name}.csv"
    logger.info(f"Loading sensitivity CSV: {csv_path}")

    df_viability, df_remaining = load_sensitivity_csv(str(csv_path))

    # moyennes sur les runs
    via_mean = df_viability.mean(axis=1)
    conc_mean = df_remaining.mean(axis=1)

    # figure
    fig, axes = plt.subplots(2, 1, figsize=(6, 7), sharex=True)

    axes[0].plot(via_mean.index, via_mean, label="Simulation", linewidth=2)
    axes[0].set_ylabel("Viability (%)")
    axes[0].set_title(f"Perturbation: {exp_name}")

    axes[1].plot(conc_mean.index, conc_mean, label="Simulation", linewidth=2)
    axes[1].set_ylabel("Concentration Ratio (%)")
    axes[1].set_xlabel("Step")

    fig.tight_layout()

    save_dir_path = Path(save_dir)
    save_dir_path.mkdir(parents=True, exist_ok=True)
    out_path = save_dir_path / f"sensitivity_{exp_name}.png"
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    logger.info(f"Saved {out_path}")


def plot_sensitivity_all_params(
    exp_list: List[str],
    csv_dir: str,
    save_dir: str,
) -> None:
    """Produit une figure par paramètre dans exp_list."""
    for exp_name in exp_list:
        plot_sensitivity_for_param(exp_name, csv_dir, save_dir)
