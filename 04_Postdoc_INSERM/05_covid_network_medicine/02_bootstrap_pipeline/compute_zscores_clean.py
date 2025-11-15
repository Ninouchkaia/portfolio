#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Input: multilayer network + mock networks provided by M. De Domenico (CoMuNe Lab)

Étape 4 du pipeline :
    - Lit les distributions mock (µ, σ, …) générées à l’Étape 3
    - Lit les valeurs observées (de l’Étape 2 : *withZeros.tsv*)
    - Pour chaque entité (GO, disease, drug, symptom, protein):
           Z = (observé - µ) / σ
    - Supporte n types d’entités dans un même run (multi-focal)
    - Exporte un TSV propre contenant :
           node_id, observed, mean, sd, min, max, median, zscore
    - Laisse Étape 5 traiter la conversion en p-values

Le module remplace :
    - toutes les variantes de fichiers *_NormTest_alpha_0.05.tsv
    - les anciens scripts de calcul de Z-scores dans /bootstrap/
"""

import argparse
import csv
import os
import numpy as np

# ---------------------------------------------------------------------
# Verbose print
# ---------------------------------------------------------------------

VERBOSE = True
def vprint(*args, **kwargs):
    if VERBOSE:
        print(*args, **kwargs)


# ---------------------------------------------------------------------
# Reading helpers : distributions & observed
# ---------------------------------------------------------------------

def load_mock_distribution_tsv(path):
    """
    Lecture du TSV créé par compute_mock_distributions_clean.py

    Colonnes attendues :
        node_id, n_mocks, mean, sd, min, max, median, values

    Retourne :
        stats[node_id] = {
            'mean': ...,
            'sd': ...,
            'min': ...,
            'max': ...,
            'median': ...,
            'n': ...,
            'values': [...]
        }
    """
    stats = {}
    vprint(f"[load_mock_distribution] Loading {path}")
    with open(path, "r", encoding="utf-8") as f:
        reader = csv.reader(f, delimiter="\t")
        header = next(reader)

        for row in reader:
            node = row[0]
            n = int(row[1])
            mean = float(row[2])
            sd = float(row[3])
            vmin = float(row[4])
            vmax = float(row[5])
            med = float(row[6])
            vals = [float(x) for x in row[7].split(",")]

            stats[node] = {
                "n": n,
                "mean": mean,
                "sd": sd,
                "min": vmin,
                "max": vmax,
                "median": med,
                "values": vals
            }
    return stats


def load_observed_withzeros(path):
    """
    Lecture du observed_withZeros TSV créé à l’Étape 2.
    Format attendu :
        source_node, degree, relative_degree

    On lit seulement "degree" comme valeur observée brute.
    """
    obs = {}
    vprint(f"[load_observed] Loading {path}")
    with open(path, "r", encoding="utf-8") as f:
        reader = csv.reader(f, delimiter="\t")
        header = next(reader)
        for row in reader:
            node = row[0]
            deg = float(row[1])
            obs[node] = deg
    return obs


# ---------------------------------------------------------------------
# Z-score computation
# ---------------------------------------------------------------------

def compute_zscores_for_focal_type(
    observed_dict,
    dist_stats_dict
):
    """
    observed_dict :  {node: observed_degree}
    dist_stats_dict : {node: {'mean':..., 'sd':..., ...}}

    Retourne :
        result[node] = {
            'observed': ...,
            'mean': ...,
            'sd': ...,
            'min': ...,
            'max': ...,
            'median': ...,
            'zscore': ...
        }
    """

    result = {}

    for node, obs_val in observed_dict.items():
        if node not in dist_stats_dict:
            # le noeud n'existe pas dans les mocks → sd=0 → pas de Z
            result[node] = {
                "observed": obs_val,
                "mean": None,
                "sd": None,
                "min": None,
                "max": None,
                "median": None,
                "zscore": None
            }
            continue

        stats = dist_stats_dict[node]
        mu = stats["mean"]
        sd = stats["sd"]

        if sd == 0 or sd is None:
            z = None
        else:
            z = (obs_val - mu) / sd

        result[node] = {
            "observed": obs_val,
            "mean": mu,
            "sd": sd,
            "min": stats["min"],
            "max": stats["max"],
            "median": stats["median"],
            "zscore": z
        }

    return result


# ---------------------------------------------------------------------
# Export
# ---------------------------------------------------------------------

def export_zscore_table(zdict, outfile):
    """
    Exporte un TSV propre :
        node_id, observed, mean, sd, min, max, median, zscore
    """
    rows = []
    for node, d in zdict.items():
        rows.append([
            node,
            d["observed"],
            d["mean"],
            d["sd"],
            d["min"],
            d["max"],
            d["median"],
            d["zscore"]
        ])

    # tri optionnel par z-score décroissant (zscore None en bas)
    def sort_key(row):
        z = row[7]
        if z is None:
            return -999999
        return z

    rows.sort(key=sort_key, reverse=True)

    os.makedirs(os.path.dirname(outfile), exist_ok=True)
    with open(outfile, "w", encoding="utf-8", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(["node_id", "observed", "mean", "sd", "min", "max", "median", "zscore"])
        writer.writerows(rows)

    vprint(f"[export_zscores] Wrote {len(rows)} rows to {outfile}")


# ---------------------------------------------------------------------
# Multiple focal types in one run
# ---------------------------------------------------------------------

def process_multiple_focal_types(focal_types, observed_files, distribution_files, outdir):
    """
    focal_types : ["GO", "disease", "drug", ...]
    observed_files :  {focal_type: path_to_observed_withZeros}
    distribution_files : {focal_type: path_to_mock_distribution}

    outdir : dossier où écrire les Z-scores.

    Pour chaque type, on calcule et exporte :
        outdir/{focal_type}_zscores.tsv
    """

    for ftype in focal_types:
        vprint("\n===================================================")
        vprint(f"[process] Focal type: {ftype}")
        vprint("===================================================\n")

        obs_file = observed_files.get(ftype)
        dist_file = distribution_files.get(ftype)

        if obs_file is None or dist_file is None:
            vprint(f"[WARNING] Missing file for {ftype}, skipping.")
            continue

        obs = load_observed_withzeros(obs_file)
        dstat = load_mock_distribution_tsv(dist_file)

        zdict = compute_zscores_for_focal_type(obs, dstat)

        outfile = os.path.join(outdir, f"{ftype}_zscores.tsv")
        export_zscore_table(zdict, outfile)


# ---------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description="Compute Z-scores from observed + mock distributions.")

    parser.add_argument("--focal_types",
                        required=True,
                        help="Comma-separated list of focal types (e.g. GO,disease,drug)")

    parser.add_argument("--observed_files",
                        required=True,
                        help="Comma-separated list of observed_withZeros files, "
                             "same order as focal_types")

    parser.add_argument("--distribution_files",
                        required=True,
                        help="Comma-separated list of mock-distribution TSVs, "
                             "same order as focal_types")

    parser.add_argument("--outdir", required=True,
                        help="Output folder for Z-score tables")

    parser.add_argument("--verbose", default="True",
                        help="True/False")

    args = parser.parse_args()

    global VERBOSE
    VERBOSE = (args.verbose.lower() == "true")

    focal_types = [x.strip() for x in args.focal_types.split(",")]
    obs_paths = [x.strip() for x in args.observed_files.split(",")]
    dist_paths = [x.strip() for x in args.distribution_files.split(",")]

    if len(focal_types) != len(obs_paths) or len(focal_types) != len(dist_paths):
        raise ValueError("Mismatch in lengths: focal_types, observed_files, distribution_files")

    observed_dict = dict(zip(focal_types, obs_paths))
    dist_dict = dict(zip(focal_types, dist_paths))

    process_multiple_focal_types(
        focal_types=focal_types,
        observed_files=observed_dict,
        distribution_files=dist_dict,
        outdir=args.outdir
    )


if __name__ == "__main__":
    main()
