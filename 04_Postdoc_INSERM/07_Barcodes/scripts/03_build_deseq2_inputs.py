#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import os
import pandas as pd


def parse_args():
    p = argparse.ArgumentParser(
        description=(
            "Construit les fichiers d'entrée pour DESeq2 : "
            "matrice de counts filtrée + table de design."
        )
    )
    p.add_argument(
        "--counts",
        default="data/processed/barcodes_filtered_counts.csv",
        help="Matrice de counts filtrée (sortie de 01_preprocess_barcodes.py)",
    )
    p.add_argument(
        "--output_dir",
        default="results/deseq2_inputs",
        help="Dossier de sortie pour counts + design",
    )
    p.add_argument(
        "--sep",
        default=";",
        help="Séparateur du fichier de counts (par défaut ;) ",
    )
    return p.parse_args()


def parse_sample_name(sample):
    """
    Suppose des noms du type:
    DrugX1_006u_exp200921_run1_xxx
    Contro1_000u_exp200921_run1_xxx
    Temps01_000u_exp200921_run1_xxx
    Adapte ce parsing si besoin.
    """
    parts = sample.split("_")
    drug_repl = parts[0] if len(parts) > 0 else ""
    dose = parts[1] if len(parts) > 1 else ""
    exp = parts[2] if len(parts) > 2 else ""
    run = parts[3] if len(parts) > 3 else ""

    replicate = drug_repl[-1] if len(drug_repl) > 0 else ""
    drug = drug_repl[:-1] if len(drug_repl) > 1 else drug_repl

    is_control = "Contro" in drug_repl
    is_timezero = ("Temps0" in sample) or ("Temps" in sample)

    if is_control:
        condition = "control"
    elif is_timezero:
        condition = "time0"
    else:
        condition = f"{drug}_{dose}"

    return {
        "sample": sample,
        "drug_raw": drug_repl,
        "drug": drug,
        "dose": dose,
        "exp": exp,
        "run": run,
        "replicate": replicate,
        "condition": condition,
        "is_control": is_control,
        "is_timezero": is_timezero,
    }


def main():
    args = parse_args()
    os.makedirs(args.output_dir, exist_ok=True)

    df_counts = pd.read_csv(args.counts, sep=args.sep, header=0, index_col=0)
    print(f"[INFO] Matrice counts filtrée : {df_counts.shape}")

    # On réécrit une version TSV pour DESeq2
    counts_path = os.path.join(args.output_dir, "counts_for_deseq2.tsv")
    df_counts.to_csv(counts_path, sep="\t")
    print(f"[INFO] Counts DESeq2 : {counts_path}")

    # Table de design
    meta_records = []
    for sample in df_counts.columns:
        meta_records.append(parse_sample_name(sample))

    df_design = pd.DataFrame.from_records(meta_records)
    design_path = os.path.join(args.output_dir, "design_for_deseq2.tsv")
    df_design.to_csv(design_path, sep="\t", index=False)
    print(f"[INFO] Design DESeq2 : {design_path}")

    print("[DONE] 03_build_deseq2_inputs terminé.")


if __name__ == "__main__":
    main()
