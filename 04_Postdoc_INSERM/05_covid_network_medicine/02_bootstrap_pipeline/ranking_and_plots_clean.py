#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Input: multilayer network + mock networks provided by M. De Domenico (CoMuNe Lab)

Étape 6 du pipeline :
    - Lit les fichiers p-values (un par focal_type)
    - Classe les entités selon :
          * p_shapiro
          * p_dagostino
          * zscore (descendant)
    - Génère :
          * un fichier TSV "topN"
          * un fichier TSV complet trié
    - Génère des barplots horizontaux (png)
          * abs values: absolute degrees / observed
          * rel values: zscore ou -log10(pvalue), selon option

Ce module remplace entièrement :
    - random-network-analysis-*-ranking.py
    - random-network-analysis-*-filtering.py
    - pulling_results.py / pulling_results_undirected.py
    - plots.py
"""

import argparse
import csv
import os
import math

import numpy as np
import matplotlib.pyplot as plt


# ---------------------------------------------------------------------
# Verbose print
# ---------------------------------------------------------------------

VERBOSE = True
def vprint(*args, **kwargs):
    if VERBOSE:
        print(*args, **kwargs)



# ---------------------------------------------------------------------
# I/O helpers
# ---------------------------------------------------------------------

def load_pvalues(path):
    """
    Charge le TSV produit par compute_pvalues_clean.py.
    Colonnes :
        node_id, sample_size, min, max, mean, median, sd,
        zscore, observed, isNormalShapiro, p_shapiro,
        isNormalDagostino, p_dagostino
    """
    vprint(f"[load_pvalues] Loading {path}")
    table = {}
    with open(path, "r", encoding="utf-8") as f:
        reader = csv.reader(f, delimiter="\t")
        header = next(reader)
        for row in reader:
            node = row[0]
            sample_size = int(row[1])
            vmin = float(row[2]) if row[2] != "None" else None
            vmax = float(row[3]) if row[3] != "None" else None
            mean = float(row[4]) if row[4] != "None" else None
            med  = float(row[5]) if row[5] != "None" else None
            sd   = float(row[6]) if row[6] != "None" else None
            zscore = float(row[7]) if row[7] != "None" else None
            observed = float(row[8]) if row[8] != "None" else None
            isN_S = row[9]
            p_S = float(row[10]) if row[10] not in ("None", "") else None
            isN_D = row[11]
            p_D = float(row[12]) if row[12] not in ("None", "") else None

            table[node] = {
                "sample_size": sample_size,
                "min": vmin,
                "max": vmax,
                "mean": mean,
                "median": med,
                "sd": sd,
                "observed": observed,
                "zscore": zscore,
                "isNormalShapiro": isN_S,
                "p_shapiro": p_S,
                "isNormalDagostino": isN_D,
                "p_dagostino": p_D
            }
    return table



# ---------------------------------------------------------------------
# Ranking
# ---------------------------------------------------------------------

def rank_table(table, metric="p_shapiro"):
    """
    Classe les entités selon un metric donné (p_shapiro / p_dagostino / zscore).
    """
    vprint(f"[rank_table] Ranking by {metric}")

    def sortkey(node):
        val = table[node][metric]
        if val is None or (isinstance(val, float) and math.isnan(val)):
            return float("inf")   # p-values None = en bas
        return val

    reverse = False  # p-values → tri ascendant
    if metric == "zscore":
        reverse = True   # zscore → tri descendant

    ranked = sorted(table.keys(), key=sortkey, reverse=reverse)
    return ranked



def write_ranked_tsv(table, ranked_nodes, outfile):
    """
    Construit un TSV trié.
    """
    vprint(f"[write_ranked_tsv] Writing {outfile}")

    os.makedirs(os.path.dirname(outfile), exist_ok=True)
    with open(outfile, "w", encoding="utf-8", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow([
            "node_id", "observed", "mean", "sd", "median",
            "zscore", "p_shapiro", "p_dagostino",
            "isNormalShapiro", "isNormalDagostino"
        ])
        for node in ranked_nodes:
            d = table[node]
            writer.writerow([
                node,
                d["observed"],
                d["mean"],
                d["sd"],
                d["median"],
                d["zscore"],
                d["p_shapiro"],
                d["p_dagostino"],
                d["isNormalShapiro"],
                d["isNormalDagostino"]
            ])



def write_topN(ranked_nodes, table, N, outfile):
    """
    Exporte le top N sous forme TSV.
    """
    vprint(f"[write_topN] Writing top{N}: {outfile}")

    os.makedirs(os.path.dirname(outfile), exist_ok=True)
    with open(outfile, "w", encoding="utf-8", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow([
            "rank", "node_id", "observed", "zscore", "p_shapiro", "p_dagostino"
        ])
        for i, node in enumerate(ranked_nodes[:N], start=1):
            d = table[node]
            writer.writerow([
                i,
                node,
                d["observed"],
                d["zscore"],
                d["p_shapiro"],
                d["p_dagostino"]
            ])



# ---------------------------------------------------------------------
# Plots
# ---------------------------------------------------------------------

def plot_topN(
    table,
    ranked_nodes,
    N,
    outfile,
    value_field="zscore",   # "zscore", "observed", "-log10(p_shapiro)"
    title="Top entities"
):
    """
    Génère un barplot horizontal simple.
    Compatible avec :
      - zscore
      - observed
      - -log10(p_shapiro)
    """
    vprint(f"[plot_topN] Plot to {outfile}")

    nodes = ranked_nodes[:N]
    values = []

    if value_field == "-log10(p_shapiro)":
        for n in nodes:
            p = table[n]["p_shapiro"]
            if p is None or p <= 0:
                values.append(0)
            else:
                values.append(-math.log10(p))
    else:
        for n in nodes:
            v = table[n][value_field]
            values.append(v if v is not None else 0)

    plt.figure(figsize=(8, max(4, N * 0.35)))
    y = np.arange(len(nodes))
    plt.barh(y, values, color="grey", edgecolor="black")

    plt.yticks(y, nodes)
    plt.gca().invert_yaxis()
    plt.title(title)
    plt.xlabel(value_field)
    plt.tight_layout()

    os.makedirs(os.path.dirname(outfile), exist_ok=True)
    plt.savefig(outfile, dpi=300)
    plt.close()



# ---------------------------------------------------------------------
# Multi-focal orchestrator
# ---------------------------------------------------------------------

def process_focal_types(
    focal_types,
    pvalue_files,
    outdir,
    topN=20,
    ranking_metric="p_shapiro",
    plot_metric="zscore"
):
    for ftype in focal_types:
        vprint("\n===================================================")
        vprint(f"[process_focal_types] {ftype}")
        vprint("===================================================\n")

        pfile = pvalue_files[ftype]
        table = load_pvalues(pfile)

        ranked_nodes = rank_table(table, metric=ranking_metric)

        # fichiers TSV
        ranked_out = os.path.join(outdir, f"{ftype}_ranking.tsv")
        write_ranked_tsv(table, ranked_nodes, ranked_out)

        top_out = os.path.join(outdir, f"{ftype}_top{topN}.tsv")
        write_topN(ranked_nodes, table, topN, top_out)

        # plot
        plot_out = os.path.join(outdir, f"{ftype}_top{topN}_{plot_metric}.png")
        plot_title = f"Top {topN} {ftype} by {plot_metric}"
        plot_topN(table, ranked_nodes, topN, plot_out, value_field=plot_metric, title=plot_title)



# ---------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description="Ranking + plots from p-values.")

    parser.add_argument("--focal_types",
                        required=True,
                        help="Comma-separated list of focal types (GO,disease,drug,...)")

    parser.add_argument("--pvalue_files",
                        required=True,
                        help="Comma-separated list of p-value TSVs (same order as focal_types)")

    parser.add_argument("--outdir", required=True,
                        help="Directory for ranking TSVs + plots")

    parser.add_argument("--topN", default="20", help="Top N entities")
    parser.add_argument("--ranking_metric", default="p_shapiro",
                        help="Metric: p_shapiro | p_dagostino | zscore")
    parser.add_argument("--plot_metric", default="zscore",
                        help="Metric for plots: zscore | observed | -log10(p_shapiro)")
    parser.add_argument("--verbose", default="True")

    args = parser.parse_args()

    global VERBOSE
    VERBOSE = (args.verbose.lower() == "true")

    focal_types = [x.strip() for x in args.focal_types.split(",")]
    pfiles = [x.strip() for x in args.pvalue_files.split(",")]

    if len(focal_types) != len(pfiles):
        raise ValueError("Mismatch: focal_types and pvalue_files length differ.")

    p_dict = dict(zip(focal_types, pfiles))

    process_focal_types(
        focal_types=focal_types,
        pvalue_files=p_dict,
        outdir=args.outdir,
        topN=int(args.topN),
        ranking_metric=args.ranking_metric,
        plot_metric=args.plot_metric
    )


if __name__ == "__main__":
    main()
