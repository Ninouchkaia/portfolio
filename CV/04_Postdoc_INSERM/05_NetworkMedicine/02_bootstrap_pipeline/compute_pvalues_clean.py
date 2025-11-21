#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Input: multilayer network + mock networks provided by M. De Domenico (CoMuNe Lab)

Étape 5 du pipeline :
    - Lit les distributions mock (µ, σ, valeurs) de l’Étape 3
    - Lit les Z-scores de l’Étape 4
    - Pour chaque entité (GO, disease, drug, symptom, protein) :
        * teste la normalité (Shapiro, D’Agostino)
        * si normal -> p = 1 - erf(|Z| / sqrt(2))
        * sinon -> p = (sd^2) / (n * |x_covid - µ|)  (borne de Chebyshev utilisée dans les scripts originaux)
    - Gère les cas particuliers :
        * sd = 0  -> z = 0, p = 0.5 si normal, 1.0 sinon
        * p > 1   -> p tronqué à 1
        * peu de points -> p = NaN, isNormal = 'NaN'

Le module remplace :
    - test_normal_distributions.py
    - transform_Zscore_in_pvalue.py
    - convert_pval_manlio.py
    - convert_pval_manlio_adjust.py
    (logique intégrée ici de façon propre et modulaire).

Usage (CLI example):
    python compute_pvalues_clean.py \
      --focal_types GO,drug \
      --zscores_files results/GO_zscores.tsv,results/drug_zscores.tsv \
      --distribution_files results/covid_GO_degree_to_proteins_distribution_in_mock_networks.tsv,results/covid_drug_degree_to_proteins_distribution_in_mock_networks.tsv \
      --outdir results/pvalues/
"""

import argparse
import csv
import os
import math

import numpy as np
from scipy.stats import shapiro, normaltest
from math import erf

from compute_zscores_clean import load_mock_distribution_tsv  # même dossier


# ---------------------------------------------------------------------
# Verbose print
# ---------------------------------------------------------------------

VERBOSE = True

def vprint(*args, **kwargs):
    if VERBOSE:
        print(*args, **kwargs)


# ---------------------------------------------------------------------
# Load Z-score TSV
# ---------------------------------------------------------------------

def load_zscore_table(path):
    """
    Lecture d'un TSV de Z-scores produit par compute_zscores_clean.py

    Colonnes attendues :
        node_id, observed, mean, sd, min, max, median, zscore

    Retourne :
        ztab[node_id] = {
            'observed': ...,
            'mean': ...,
            'sd': ...,
            'min': ...,
            'max': ...,
            'median': ...,
            'zscore': ...
        }
    """
    ztab = {}
    vprint(f"[load_zscore_table] Loading {path}")
    with open(path, "r", encoding="utf-8") as f:
        reader = csv.reader(f, delimiter="\t")
        header = next(reader)
        for row in reader:
            node = row[0]
            obs = float(row[1])
            mean = float(row[2]) if row[2] not in ("", "None") else None
            sd   = float(row[3]) if row[3] not in ("", "None") else None
            vmin = float(row[4]) if row[4] not in ("", "None") else None
            vmax = float(row[5]) if row[5] not in ("", "None") else None
            med  = float(row[6]) if row[6] not in ("", "None") else None
            z    = row[7]
            z    = float(z) if z not in ("", "None") else None

            ztab[node] = {
                "observed": obs,
                "mean": mean,
                "sd": sd,
                "min": vmin,
                "max": vmax,
                "median": med,
                "zscore": z
            }
    return ztab


# ---------------------------------------------------------------------
# p-value computation (Shapiro + D’Agostino)
# ---------------------------------------------------------------------

def _pvalue_from_z_chebyshev(sd, n, observed, mean):
    """
    Variante Chebyshev, cohérente avec les anciens scripts :
        p = (sd^2) / (n * |observed - mean|)
    avec traitement spécial deg_covid == mean -> p = 0.5
    (mais ce cas est géré à l'extérieur).
    """
    if n <= 0:
        return float("nan")
    if observed == mean:
        return 0.5
    return (sd * sd) / (n * abs(observed - mean))


def _pvalues_for_node(distribution, z, observed, mean, sd,
                      alpha=0.05):
    """
    distribution : liste de valeurs (mock)
    z           : Z-score (float ou None)
    observed    : valeur dans le réseau réel
    mean, sd    : µ et σ de la distribution mock

    Retourne :
        {
          'isNormalShapiro': <str>,
          'p_shapiro': float,
          'isNormalDagostino': <str>,
          'p_dagostino': float
        }
    suivant la logique originale :
        - test Shapiro si n > 3
        - test D’Agostino si n > 8
        - normal -> erf
        - non-normal -> Chebyshev
    """
    n = len(distribution)

    # Valeurs par défaut
    isN_S = "NaN"
    p_S = float("nan")
    isN_D = "NaN"
    p_D = float("nan")

    # Shapiro
    if n > 3:
        stat, p = shapiro(distribution)
        if p < alpha:
            isNormalShapiro = False
        else:
            isNormalShapiro = True
        isN_S = str(isNormalShapiro)

        if isNormalShapiro:
            # cas normal
            if z is None or math.isnan(z):
                p_S = 0.5
            else:
                p_S = 1.0 - erf(abs(z) / math.sqrt(2.0))
        else:
            # Chebyshev
            if (observed == mean):
                p_S = 0.5
            else:
                p_S = _pvalue_from_z_chebyshev(sd, n, observed, mean)
    else:
        p_S = float("nan")
        isN_S = "NaN"

    # D’Agostino
    if n > 8:
        stat, p = normaltest(distribution)
        if p < alpha:
            isNormalDagostino = False
        else:
            isNormalDagostino = True
        isN_D = str(isNormalDagostino)

        if isNormalDagostino:
            if z is None or math.isnan(z):
                p_D = 0.5
            else:
                p_D = 1.0 - erf(abs(z) / math.sqrt(2.0))
        else:
            if (observed == mean):
                p_D = 0.5
            else:
                p_D = _pvalue_from_z_chebyshev(sd, n, observed, mean)
    else:
        p_D = float("nan")
        isN_D = "NaN"

    return {
        "isNormalShapiro": isN_S,
        "p_shapiro": p_S,
        "isNormalDagostino": isN_D,
        "p_dagostino": p_D
    }


def _adjust_special_cases(sd, z, isN_S, p_S, isN_D, p_D):
    """
    Implémente la logique de convert_pval_manlio_adjust.py :
    - si sd == 0 → z = 0, p_shap/p_dago = 0.5 (si normal) ou 1 (si non-normal)
    - p > 1 → tronqué à 1
    """
    # sd == 0
    if sd == 0.0:
        z = 0.0
        # Shapiro
        if isN_S == "True":
            p_S = 0.5
        elif isN_S == "False":
            p_S = 1.0
        # D’Agostino
        if isN_D == "True":
            p_D = 0.5
        elif isN_D == "False":
            p_D = 1.0

    # p > 1 → tronquer
    if p_S is not None and not math.isnan(p_S) and p_S > 1.0:
        p_S = 1.0
    if p_D is not None and not math.isnan(p_D) and p_D > 1.0:
        p_D = 1.0

    return z, p_S, p_D


# ---------------------------------------------------------------------
# Process focal type
# ---------------------------------------------------------------------

def compute_pvalues_for_focal_type(zscore_file, distribution_file, outfile, alpha=0.05):
    """
    Un focal_type (GO, disease, drug, symptom, etc.) :

    - lit distributions mock (mean, sd, val list,...)
    - lit observed + zscores
    - calcule p_shapiro & p_dagostino
    - applique l'ajustement sd==0 & p>1
    - écrit outfile TSV
    """
    ztab = load_zscore_table(zscore_file)
    dist_stats = load_mock_distribution_tsv(distribution_file)

    rows = []

    for node, zd in ztab.items():
        obs = zd["observed"]
        mean = zd["mean"]
        sd = zd["sd"]
        z = zd["zscore"]

        if node not in dist_stats or mean is None or sd is None:
            # manque d’info mock → tout NaN
            sample_size = dist_stats.get(node, {}).get("n", 0)
            vmin = zd["min"]
            vmax = zd["max"]
            med = zd["median"]
            isN_S = "NaN"
            p_S = float("nan")
            isN_D = "NaN"
            p_D = float("nan")
        else:
            dstat = dist_stats[node]
            sample_size = dstat["n"]
            vmin = dstat["min"]
            vmax = dstat["max"]
            med = dstat["median"]
            distribution = dstat["values"]

            res = _pvalues_for_node(
                distribution=distribution,
                z=z,
                observed=obs,
                mean=mean,
                sd=sd,
                alpha=alpha
            )

            isN_S = res["isNormalShapiro"]
            p_S = res["p_shapiro"]
            isN_D = res["isNormalDagostino"]
            p_D = res["p_dagostino"]

            z, p_S, p_D = _adjust_special_cases(sd, z, isN_S, p_S, isN_D, p_D)

        # ligne finale
        rows.append([
            node,
            sample_size,
            vmin,
            vmax,
            mean,
            med,
            sd,
            z,
            obs,
            isN_S,
            p_S,
            isN_D,
            p_D
        ])

    # tri par p_shapiro croissant (les plus significatifs en haut)
    def sort_key(row):
        pS = row[10]
        if pS is None or (isinstance(pS, float) and math.isnan(pS)):
            return 1e9
        return pS

    rows.sort(key=sort_key)

    os.makedirs(os.path.dirname(outfile), exist_ok=True)
    with open(outfile, "w", encoding="utf-8", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow([
            "node_id",
            "sample_size",
            "min",
            "max",
            "mean",
            "median",
            "sd",
            "zscore",
            "observed",
            "isNormalShapiro",
            "p_shapiro",
            "isNormalDagostino",
            "p_dagostino"
        ])
        writer.writerows(rows)

    vprint(f"[compute_pvalues_for_focal_type] Wrote {len(rows)} rows to {outfile}")


# ---------------------------------------------------------------------
# Multi-focal CLI
# ---------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description="Compute p-values from Z-scores + mock distributions (Shapiro + D’Agostino, erf/Chebyshev).")

    parser.add_argument("--focal_types",
                        required=True,
                        help="Comma-separated list of focal types (e.g. GO,disease,drug)")

    parser.add_argument("--zscores_files",
                        required=True,
                        help="Comma-separated list of Z-score TSV files, "
                             "same order as focal_types")

    parser.add_argument("--distribution_files",
                        required=True,
                        help="Comma-separated list of mock-distribution TSVs, "
                             "same order as focal_types")

    parser.add_argument("--outdir", required=True,
                        help="Output folder for p-value tables")

    parser.add_argument("--alpha", default="0.05",
                        help="Alpha for normality tests (default 0.05)")

    parser.add_argument("--verbose", default="True",
                        help="True/False")

    args = parser.parse_args()

    global VERBOSE
    VERBOSE = (args.verbose.lower() == "true")
    alpha = float(args.alpha)

    focal_types = [x.strip() for x in args.focal_types.split(",")]
    z_paths = [x.strip() for x in args.zscores_files.split(",")]
    dist_paths = [x.strip() for x in args.distribution_files.split(",")]

    if len(focal_types) != len(z_paths) or len(focal_types) != len(dist_paths):
        raise ValueError("Mismatch in lengths: focal_types, zscores_files, distribution_files")

    z_dict = dict(zip(focal_types, z_paths))
    dist_dict = dict(zip(focal_types, dist_paths))

    for ftype in focal_types:
        vprint("\n===================================================")
        vprint(f"[main] Focal type: {ftype}")
        vprint("===================================================\n")

        zfile = z_dict[ftype]
        dfile = dist_dict[ftype]

        outfile = os.path.join(args.outdir, f"{ftype}_pvalues.tsv")
        compute_pvalues_for_focal_type(zfile, dfile, outfile, alpha=alpha)


if __name__ == "__main__":
    main()
