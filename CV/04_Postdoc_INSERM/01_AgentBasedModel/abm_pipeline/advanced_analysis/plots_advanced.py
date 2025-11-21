# abm_pipeline/advanced_analysis/plots_advanced.py

from __future__ import annotations
from pathlib import Path
from typing import List, Optional

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

from abm_pipeline.advanced_analysis.utils import (
    load_pareto_dataframe,
    ensure_dir,
)
from abm_pipeline.config import ADVANCED_RESULTS_DIR
from abm_pipeline.parameter_exploration.utils import logger


def make_violinplot(
    df: pd.DataFrame,
    param: str,
    save_dir: Optional[str] = None,
    title: Optional[str] = None,
):
    """
    Produit le violin plot d’un paramètre.
    """
    if save_dir is None:
        save_path = ensure_dir(ADVANCED_RESULTS_DIR / "violin")
    else:
        save_path = ensure_dir(save_dir)

    fig, ax = plt.subplots(figsize=(6, 5))

    sns.violinplot(
        data=df,
        x=param,
        palette="Set2",  
        ax=ax,
    )

    if title:
        ax.set_title(title)

    out_path = save_path / f"violin_{param}.png"
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    logger.info(f"Saved violin plot: {out_path}")


def make_violinplots_all_parameters(
    pareto_file: str,
    save_dir: Optional[str] = None,
    exclude_cols: Optional[List[str]] = None,
):
    """
    Produit un violin plot pour chaque paramètre du fichier Pareto.
    Equivalent de parameter_violin_plots.py + all_pareto_df_to_violin.py
    """
    df = load_pareto_dataframe(pareto_file)

    if exclude_cols is None:
        # À adapter selon ton fichier Pareto exact
        exclude_cols = [
            "delta_fitness_via",
            "delta_fitness_conc",
        ]

    params = [col for col in df.columns if col not in exclude_cols]

    for param in params:
        make_violinplot(df, param, save_dir=save_dir, title=f"Distribution of {param}")




from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA


def run_pca_analysis(
    pareto_file: str,
    save_dir: Optional[str] = None,
):
    """
    Effectue une PCA sur les paramètres optimisés.
    Equivalent de pca_analysis.py
    (Seulement la partie ‘projection 2D’ + scree plot.)
    """
    df = load_pareto_dataframe(pareto_file)

    if save_dir is None:
        save_path = ensure_dir(ADVANCED_RESULTS_DIR / "pca")
    else:
        save_path = ensure_dir(save_dir)

    # On enlève fitness + éventuelles col. non-param
    df_params = df.drop(columns=["delta_fitness_via", "delta_fitness_conc"], errors="ignore")

    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(df_params)

    pca = PCA(n_components=2)
    comps = pca.fit_transform(X_scaled)

    fig, ax = plt.subplots(figsize=(6, 6))
    ax.scatter(comps[:, 0], comps[:, 1], s=30, alpha=0.7)
    ax.set_xlabel("PC1")
    ax.set_ylabel("PC2")
    ax.set_title("PCA of Optimized Parameters")

    fig.savefig(save_path / "pca_scatter.png", dpi=300, bbox_inches="tight")

    # Scree plot
    fig2, ax2 = plt.subplots(figsize=(6, 4))
    ax2.plot(pca.explained_variance_ratio_, marker="o")
    ax2.set_title("Explained Variance Ratio")
    ax2.set_xlabel("Principal Component")
    ax2.set_ylabel("Variance Ratio")

    fig2.savefig(save_path / "pca_scree.png", dpi=300, bbox_inches="tight")

    logger.info(f"PCA results saved to {save_path}")
