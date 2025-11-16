#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import glob
import os
import pandas as pd


def parse_args():
    p = argparse.ArgumentParser(
        description=(
            "Fusionne les fichiers *run*.csv de comptage barcodes, "
            "filtre les barcodes peu détectés, et produit les matrices filtrées "
            "en counts bruts et normalisés à 1e6 par échantillon."
        )
    )
    p.add_argument(
        "--input_dir",
        default="data/raw",
        help="Dossier contenant les fichiers *run*.csv",
    )
    p.add_argument(
        "--pattern",
        default="*run*.csv",
        help="Pattern des fichiers de comptage (par défaut *run*.csv)",
    )
    p.add_argument(
        "--output_prefix",
        default="data/processed/barcodes",
        help="Préfixe des fichiers de sortie (sans extension)",
    )
    p.add_argument(
        "--sep",
        default=";",
        help="Séparateur des CSV d'entrée (par défaut ;) ",
    )
    p.add_argument(
        "--min_reads",
        type=float,
        default=1.0,
        help="Seuil minimal de reads pour considérer qu'un barcode est détecté dans un échantillon",
    )
    p.add_argument(
        "--min_controls",
        type=int,
        default=5,
        help="Nombre minimal d'échantillons contrôle où le barcode doit être détecté",
    )
    p.add_argument(
        "--min_timezeros",
        type=int,
        default=5,
        help="Nombre minimal d'échantillons Temps0 où le barcode doit être détecté",
    )
    return p.parse_args()


def load_and_merge_runs(input_dir, pattern, sep):
    pattern_path = os.path.join(input_dir, pattern)
    files = sorted(glob.glob(pattern_path))
    if not files:
        raise FileNotFoundError(f"Aucun fichier ne matche {pattern_path}")

    dfs = []
    for f in files:
        print(f"[INFO] Lecture de {f}")
        df = pd.read_csv(f, sep=sep, header=0, index_col=0)
        dfs.append(df)

    # jointure sur les barcodes (index)
    merged = pd.concat(dfs, axis=1)
    merged = merged.fillna(0)
    # s'assurer que c'est bien numérique
    merged = merged.apply(pd.to_numeric, errors="coerce").fillna(0)
    print(f"[INFO] Matrice fusionnée : {merged.shape[0]} barcodes x {merged.shape[1]} échantillons")
    return merged


def filter_barcodes(df_raw, min_reads, min_controls, min_timezeros):
    cols_control = [c for c in df_raw.columns if "Contro" in c]
    cols_time0 = [c for c in df_raw.columns if "Temps" in c]

    if not cols_control:
        print("[WARN] Aucun échantillon contrôle trouvé (cols contenant 'Contro').")
    if not cols_time0:
        print("[WARN] Aucun échantillon Temps0 trouvé (cols contenant 'Temps').")

    print(f"[INFO] # contrôles : {len(cols_control)}, # temps0 : {len(cols_time0)}")

    detected_controls = (df_raw[cols_control] >= min_reads).sum(axis=1) if cols_control else 0
    detected_time0 = (df_raw[cols_time0] >= min_reads).sum(axis=1) if cols_time0 else 0

    mask = (detected_controls >= min_controls) & (detected_time0 >= min_timezeros)
    df_filt = df_raw.loc[mask].copy()

    print(f"[INFO] Après filtrage : {df_filt.shape[0]} barcodes conservés.")
    return df_filt


def normalize_cpm(df_counts):
    # normalisation par colonne (échantillon)
    lib_sizes = df_counts.sum(axis=0)
    cpm = df_counts.div(lib_sizes, axis=1) * 1_000_000
    return cpm


def main():
    args = parse_args()

    os.makedirs(os.path.dirname(args.output_prefix), exist_ok=True)

    df_raw = load_and_merge_runs(args.input_dir, args.pattern, args.sep)
    out_raw = f"{args.output_prefix}_combined_raw_counts.csv"
    df_raw.to_csv(out_raw, sep=";")
    print(f"[INFO] Sauvegardé : {out_raw}")

    df_filt = filter_barcodes(
        df_raw,
        min_reads=args.min_reads,
        min_controls=args.min_controls,
        min_timezeros=args.min_timezeros,
    )

    out_filt = f"{args.output_prefix}_filtered_counts.csv"
    df_filt.to_csv(out_filt, sep=";")
    print(f"[INFO] Sauvegardé : {out_filt}")

    df_cpm = normalize_cpm(df_filt)
    out_cpm = f"{args.output_prefix}_filtered_cpm.csv"
    df_cpm.to_csv(out_cpm, sep=";")
    print(f"[INFO] Sauvegardé : {out_cpm}")

    print("[DONE] 01_preprocess_barcodes terminé.")


if __name__ == "__main__":
    main()

