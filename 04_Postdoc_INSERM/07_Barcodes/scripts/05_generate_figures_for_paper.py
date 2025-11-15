#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Generate paper figures for the barcoding analysis:
- Fig.3: drug clustering based on barcode signatures (log2FC)
- Fig.4: drug–drug correlation clustered heatmap
- Fig.5: drug network (correlation > threshold)
"""

import os
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import networkx as nx

# -----------------------------
# Configuration
# -----------------------------

DATA_DIR = os.path.join("data")
OUT_DIR = os.path.join("results", "paper_figures")
os.makedirs(OUT_DIR, exist_ok=True)

LOGFC_FILE = os.path.join(DATA_DIR, "merged_pval_filtered_deseq2_fillna.csv")
CORR_FILE  = os.path.join(DATA_DIR, "correl_matrix1_merged_pval_filtered_deseq2.csv")
COLOR_MAP_FILE = os.path.join(DATA_DIR, "color_mapping.tsv")

CORR_THRESHOLD = 0.8  # for the network edges

# -----------------------------
# Helpers
# -----------------------------

def load_logfc_matrix(path: str) -> pd.DataFrame:
    df = pd.read_csv(path, sep=";", header=0, index_col=0)
    # should already be log2FC from DESeq2
    # enforce float + drop all-NaN cols/rows
    df = df.astype(float)
    df = df.dropna(axis=0, how="all").dropna(axis=1, how="all")
    # assert no duplicate cols
    if df.columns.duplicated().any():
        dup = df.columns[df.columns.duplicated()].tolist()
        raise ValueError(f"Duplicate condition columns in logFC matrix: {dup}")
    return df

def load_corr_matrix(logfc_df: pd.DataFrame, corr_file: str | None = None) -> pd.DataFrame:
    if corr_file and os.path.exists(corr_file):
        corr = pd.read_csv(corr_file, sep=";", header=0, index_col=0)
        return corr.astype(float)
    # else compute from logfc
    corr = logfc_df.corr(method="pearson", min_periods=1)
    return corr

def load_color_mapping(path: str) -> dict:
    """
    color_mapping.tsv:
       <int_color_code>\t<drug_label>
    Example: 1  Gefitinib_006u
    Returned mapping: drug_name (without replicate) -> color float in [0,1]
    """
    cmap = {}
    if not os.path.exists(path):
        return cmap

    with open(path, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) != 2:
                continue
            color_int, drug_label = parts
            color_int = int(color_int)

            # harmonize some special cases as in networks.py
            if drug_label == "Fluor_006u":
                drug_key = "5Fluor_006u"
            elif drug_label == "Azacyt_1,5u":
                drug_key = "Azacyt_1.5u"
            elif drug_label == "Bafilo_1,2n":
                drug_key = "Bafilo_1.2n"
            else:
                drug_key = drug_label

            # scale to [0, 1] for colormap
            cmap[drug_key] = (color_int + 1) / 100.0

    return cmap

def condition_to_drug_label(cond_name: str) -> str:
    """
    Reduce a condition name like 'GefitinibA_006u_exp200921_run1_sample1'
    to a drug label useful for coloring:
      - keep drug + dose, ignore exp/run/sample/replicate
    """
    parts = cond_name.split("_")
    if len(parts) < 2:
        return cond_name

    head = parts[0]
    dose = parts[1]

    # Drop replicate letter for drugs (last char)
    if head.startswith("Ctrl") or "Contro" in head:
        drug = "CtrlMs"
    elif head.startswith("Temps0") or "Temps0" in head:
        drug = "Temps0"
    else:
        # e.g. GefitinibA → Gefitinib
        drug = head[:-1]

    return f"{drug}_{dose}"

# -----------------------------
# Figure 3 – log2FC heatmap (drug signatures)
# -----------------------------

def figure3_drug_signature_heatmap(logfc_df: pd.DataFrame, out_dir: str):
    """
    Clustering of conditions based on barcode log2FC signatures.
    """
    # Optionally, center by barcode (row-wise) to normalize
    # Here: we keep raw log2FC, but clip extreme values for visual clarity.
    clipped = logfc_df.clip(lower=-4, upper=4)

    # Column order = clustered; row order = clustered
    sns.set(style="white")
    g = sns.clustermap(
        clipped,
        method="complete",
        metric="euclidean",
        cmap="RdBu_r",
        center=0.0,
        xticklabels=False,
        yticklabels=False,
        robust=False,
        figsize=(10, 12),
    )
    g.fig.suptitle("Annotated drug clustering based on barcode signatures (log2FC)", y=1.02)

    pdf_path = os.path.join(out_dir, "Fig3_drug_signature_heatmap.pdf")
    png_path = os.path.join(out_dir, "Fig3_drug_signature_heatmap.png")
    g.fig.savefig(pdf_path, bbox_inches="tight")
    g.fig.savefig(png_path, dpi=300, bbox_inches="tight")
    plt.close(g.fig)

# -----------------------------
# Figure 4 – Correlation heatmap
# -----------------------------

def figure4_correlation_heatmap(corr: pd.DataFrame, out_dir: str):
    sns.set(style="white")
    # make sure diagonal is exactly 1
    np.fill_diagonal(corr.values, 1.0)

    g = sns.clustermap(
        corr,
        method="complete",
        metric="euclidean",
        cmap="RdBu_r",
        vmin=-1.0,
        vmax=1.0,
        center=0.0,
        xticklabels=False,
        yticklabels=False,
        figsize=(10, 10),
    )
    g.fig.suptitle("Drug–drug correlation clustered heatmap (Pearson r)", y=1.02)

    pdf_path = os.path.join(out_dir, "Fig4_drug_correlation_heatmap.pdf")
    png_path = os.path.join(out_dir, "Fig4_drug_correlation_heatmap.png")
    g.fig.savefig(pdf_path, bbox_inches="tight")
    g.fig.savefig(png_path, dpi=300, bbox_inches="tight")
    plt.close(g.fig)

# -----------------------------
# Figure 5 – Correlation network
# -----------------------------

def figure5_correlation_network(
    corr: pd.DataFrame,
    color_mapping: dict,
    out_dir: str,
    threshold: float = 0.8,
):
    # Build edge list from upper triangle
    edges = []
    conds = corr.index.tolist()
    for i, s in enumerate(conds):
        for j in range(i + 1, len(conds)):
            t = conds[j]
            r = corr.iloc[i, j]
            if np.isnan(r):
                continue
            if r >= threshold:
                edges.append((s, t, float(r)))

    G = nx.Graph()
    for s, t, r in edges:
        G.add_edge(s, t, weight=r)

    # Node colors based on drug classes
    node_colors = []
    for n in G.nodes():
        drug_label = condition_to_drug_label(n)
        c = color_mapping.get(drug_label, 0.5)  # default mid-grey
        node_colors.append(c)

    # Positions – Fruchterman–Reingold (spring layout)
    pos = nx.spring_layout(G, k=0.15, iterations=200, seed=42)

    plt.figure(figsize=(10, 10))
    nx.draw_networkx_edges(
        G,
        pos,
        edgelist=G.edges(),
        width=0.5,
        alpha=0.6,
        edge_color="black",
    )

    nx.draw_networkx_nodes(
        G,
        pos,
        node_size=300,
        node_color=node_colors,
        cmap="viridis",
        alpha=0.9,
    )

    # Labels: drug only (no exp/run/sample)
    labels = {n: condition_to_drug_label(n).split("_")[0] for n in G.nodes()}
    nx.draw_networkx_labels(G, pos, labels=labels, font_size=6)

    plt.title(f"Drug network (Pearson r ≥ {threshold})")
    plt.axis("off")

    pdf_path = os.path.join(out_dir, "Fig5_drug_network.pdf")
    png_path = os.path.join(out_dir, "Fig5_drug_network.png")
    plt.savefig(pdf_path, bbox_inches="tight")
    plt.savefig(png_path, dpi=300, bbox_inches="tight")
    plt.close()

# -----------------------------
# Main
# -----------------------------

def main():
    print("Loading log2FC matrix…")
    logfc = load_logfc_matrix(LOGFC_FILE)

    print("Loading / computing correlation matrix…")
    corr = load_corr_matrix(logfc, CORR_FILE)

    print("Loading color mapping…")
    cmap = load_color_mapping(COLOR_MAP_FILE)

    print("Generating Fig.3 (heatmap of barcode signatures)…")
    figure3_drug_signature_heatmap(logfc, OUT_DIR)

    print("Generating Fig.4 (correlation clustered heatmap)…")
    figure4_correlation_heatmap(corr, OUT_DIR)

    print("Generating Fig.5 (correlation network)…")
    figure5_correlation_network(corr, cmap, OUT_DIR, threshold=CORR_THRESHOLD)

    print("Done. Figures saved in:", OUT_DIR)

if __name__ == "__main__":
    main()
