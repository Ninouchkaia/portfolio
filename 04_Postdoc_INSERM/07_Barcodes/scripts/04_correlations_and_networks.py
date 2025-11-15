#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import glob
import os

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import networkx as nx


def parse_args():
    p = argparse.ArgumentParser(
        description=(
            "Fusionne les résultats DESeq2 (log2FC) en une matrice, "
            "calcule une matrice de corrélation condition-condition, "
            "et génère une heatmap + un réseau."
        )
    )
    p.add_argument(
        "--deseq2_dir",
        default="results/deseq2",
        help="Dossier contenant les sorties DESeq2 (fichiers .tsv)",
    )
    p.add_argument(
        "--pattern",
        default="*.tsv",
        help="Pattern des fichiers DESeq2 (par défaut *.tsv)",
    )
    p.add_argument(
        "--output_prefix",
        default="results/networks/drug_signatures",
        help="Préfixe des fichiers de sortie",
    )
    p.add_argument(
        "--corr_threshold",
        type=float,
        default=0.8,
        help="Seuil de corrélation pour dessiner une arête dans le réseau",
    )
    return p.parse_args()


def build_logfc_matrix(deseq2_dir, pattern):
    pattern_path = os.path.join(deseq2_dir, pattern)
    files = sorted(glob.glob(pattern_path))
    if not files:
        raise FileNotFoundError(f"Aucun fichier ne matche {pattern_path}")

    dfs = []
    colnames = []

    for f in files:
        print(f"[INFO] Lecture DESeq2 : {f}")
        # nom de colonne = nom du fichier sans extension
        base = os.path.basename(f)
        cond_name = os.path.splitext(base)[0]

        df = pd.read_csv(f, sep="\t", header=0, index_col=0)
        if "log2FoldChange" not in df.columns:
            print(f"[WARN] Pas de colonne 'log2FoldChange' dans {f}, ignoré.")
            continue

        logfc = df["log2FoldChange"]
        dfs.append(logfc)
        colnames.append(cond_name)

    if not dfs:
        raise RuntimeError("Aucun fichier n'a fourni de log2FoldChange exploitable.")

    merged = pd.concat(dfs, axis=1)
    merged.columns = colnames
    merged = merged.astype(float)
    print(f"[INFO] Matrice log2FC : {merged.shape[0]} barcodes x {merged.shape[1]} conditions")
    return merged


def plot_correlation_heatmap(corr, out_png):
    sns.set(style="white")
    plt.figure(figsize=(10, 8))
    g = sns.clustermap(
        corr,
        cmap="RdBu_r",
        vmin=-1,
        vmax=1,
        linewidths=0,
        xticklabels=True,
        yticklabels=True,
    )
    plt.setp(g.ax_heatmap.get_xticklabels(), rotation=90, fontsize=6)
    plt.setp(g.ax_heatmap.get_yticklabels(), fontsize=6)
    g.fig.suptitle("Drug-drug Pearson correlation (log2FC signatures)", y=1.02)
    g.savefig(out_png, dpi=300)
    plt.close()
    print(f"[INFO] Heatmap de corrélation sauvegardée : {out_png}")


def build_and_plot_network(corr, threshold, out_png):
    links = corr.stack().reset_index()
    links.columns = ["source", "target", "correlation"]

    links_filtered = links[
        (links["source"] != links["target"]) & (links["correlation"] >= threshold)
    ]

    print(f"[INFO] # arêtes avec corr >= {threshold} : {links_filtered.shape[0]}")

    G = nx.from_pandas_edgelist(
        links_filtered, "source", "target", edge_attr="correlation"
    )

    plt.figure(figsize=(10, 10))
    pos = nx.spring_layout(G, k=0.2, seed=42)

    # couleur des noeuds en fonction du degré (juste pour avoir un gradient)
    degrees = dict(G.degree())
    deg_vals = np.array(list(degrees.values()), dtype=float)
    if len(deg_vals) > 0:
        deg_norm = (deg_vals - deg_vals.min()) / (deg_vals.ptp() + 1e-9)
    else:
        deg_norm = deg_vals

    node_colors = deg_norm

    nx.draw_networkx_nodes(G, pos, node_size=300, node_color=node_colors, cmap="viridis")
    nx.draw_networkx_edges(G, pos, width=0.5, alpha=0.7)
    nx.draw_networkx_labels(G, pos, font_size=6)

    plt.axis("off")
    plt.title(f"Drug network (corr >= {threshold})", fontsize=12)
    plt.tight_layout()
    plt.savefig(out_png, dpi=300)
    plt.close()
    print(f"[INFO] Réseau sauvegardé : {out_png}")


def main():
    args = parse_args()
    os.makedirs(os.path.dirname(args.output_prefix), exist_ok=True)

    logfc = build_logfc_matrix(args.deseq2_dir, args.pattern)

    out_logfc = f"{args.output_prefix}_logfc_matrix.csv"
    logfc.to_csv(out_logfc, sep=";")
    print(f"[INFO] Matrice log2FC sauvegardée : {out_logfc}")

    corr = logfc.corr(method="pearson", min_periods=1)

    out_corr = f"{args.output_prefix}_correlation_matrix.csv"
    corr.to_csv(out_corr, sep=";")
    print(f"[INFO] Matrice de corrélation sauvegardée : {out_corr}")

    out_heatmap = f"{args.output_prefix}_correlation_heatmap.png"
    plot_correlation_heatmap(corr, out_heatmap)

    out_network = f"{args.output_prefix}_network.png"
    build_and_plot_network(corr, args.corr_threshold, out_network)

    print("[DONE] 04_correlations_and_networks terminé.")


if __name__ == "__main__":
    main()
