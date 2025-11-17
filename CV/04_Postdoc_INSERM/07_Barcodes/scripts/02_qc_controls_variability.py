#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import os
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


def parse_args():
    p = argparse.ArgumentParser(
        description=(
            "Calcule la variabilité des barcodes dans les contrôles par expérience, "
            "et génère un TSV + un violonplot global."
        )
    )
    p.add_argument(
        "--counts",
        default="data/processed/barcodes_filtered_counts.csv",
        help="Matrice de counts filtrée (sortie de 01_preprocess_barcodes.py)",
    )
    p.add_argument(
        "--output_prefix",
        default="results/qc_controls/controls_variability",
        help="Préfixe des fichiers de sortie",
    )
    p.add_argument(
        "--sep",
        default=";",
        help="Séparateur du fichier de counts (par défaut ;) ",
    )
    return p.parse_args()


def extract_experiment(colname):
    parts = colname.split("_")
    if len(parts) >= 3:
        return parts[2]
    return "unknown"


def compute_controls_stats(df):
    records = []
    for col in df.columns:
        if "Contro" not in col:
            continue
        exp = extract_experiment(col)
        # on groupera par exp + barcode ensuite
    # On re-boucle de manière structurée
    experiments = sorted({extract_experiment(c) for c in df.columns if "Contro" in c})

    for exp in experiments:
        ctrl_cols = [c for c in df.columns
                    if ("Contro" in c) and (extract_experiment(c) == exp)]
        if not ctrl_cols:
            continue
        sub = df[ctrl_cols]

        for barcode, row in sub.iterrows():
            values = row.values.astype(float)
            avg = float(np.mean(values))
            std = float(np.std(values, ddof=1)) if len(values) > 1 else 0.0
            if avg > 0:
                max_min_ratio = float((values.max() - values.min()) / avg)
            else:
                max_min_ratio = 0.0
            records.append(
                {
                    "experiment": exp,
                    "barcode": barcode,
                    "n_controls": len(values),
                    "mean_reads": avg,
                    "std_reads": std,
                    "max_min_over_mean": max_min_ratio,
                }
            )

    df_stats = pd.DataFrame.from_records(records)
    return df_stats


def plot_violin(df_stats, output_png):
    sns.set(style="whitegrid")
    plt.figure(figsize=(8, 4))
    ax = sns.violinplot(x="max_min_over_mean", data=df_stats, bw=0.05)
    ax.set_xlabel("(max - min) / moyenne des reads (contrôles)")
    ax.set_title("Variabilité des barcodes dans les contrôles (tous expériences confondues)")
    plt.tight_layout()
    plt.savefig(output_png, dpi=300)
    plt.close()
    print(f"[INFO] Violonplot sauvegardé : {output_png}")


def main():
    args = parse_args()
    os.makedirs(os.path.dirname(args.output_prefix), exist_ok=True)

    df = pd.read_csv(args.counts, sep=args.sep, header=0, index_col=0)
    print(f"[INFO] Matrice counts : {df.shape[0]} barcodes x {df.shape[1]} échantillons")

    df_stats = compute_controls_stats(df)

    out_tsv = f"{args.output_prefix}_per_barcode.tsv"
    df_stats.to_csv(out_tsv, sep="\t", index=False)
    print(f"[INFO] Stats sauvegardées : {out_tsv}")

    out_png = f"{args.output_prefix}_violin.png"
    plot_violin(df_stats, out_png)

    print("[DONE] 02_qc_controls_variability terminé.")


if __name__ == "__main__":
    main()
