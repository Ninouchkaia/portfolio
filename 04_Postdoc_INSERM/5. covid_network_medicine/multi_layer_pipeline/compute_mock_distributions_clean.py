#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Input: multilayer network + mock networks provided by M. De Domenico (CoMuNe Lab)

Étape 3 du pipeline :
    - Parcourt l'ensemble des réseaux mock (nodes_XXXX.csv / edges_XXXX.csv)
    - Calcule, pour chaque nœud d'un type donné (par ex. GO, disease, drug, symptom),
      son degré vers un autre type de nœud (souvent protein), dans chaque mock.
    - Construit une distribution de degrés par nœud (liste de valeurs sur tous les mocks).
    - Calcule des statistiques descriptives (n, mean, sd, min, max, median).
    - Exporte un TSV de distributions pour usage ultérieur (Z-scores, tests de normalité, p-values).

Ce module remplace les anciens scripts de type :
    - get_structural_degrees_distributions_across_mock_networks.py
    - covid_*_degree_to_proteins_distribution_in_mock_networks.tsv
    - toutes les variantes par type d'entité (GO, disease, drug, symptom, etc.)

Usage (CLI example):
    python compute_mock_distributions_clean.py \
      --mock_basepath Mock_networks_21Apr/ \
      --node_template_file COVID19_GDDS_nodes.csv \
      --focal_type GO \
      --neighbor_type protein \
      --directed True \
      --outfile results/covid_GO_degree_to_proteins_distribution_in_mock_networks.tsv
"""

import argparse
import csv
import os
from collections import defaultdict

import numpy as np
import networkx as nx

from io_networks_clean import load_mock_network_ids   # même dossier


# -------------------------------------------------------------------------
# Verbose print
# -------------------------------------------------------------------------

VERBOSE = True

def vprint(*args, **kwargs):
    if VERBOSE:
        print(*args, **kwargs)


# -------------------------------------------------------------------------
# I/O helpers
# -------------------------------------------------------------------------

def load_node_types(node_file):
    """
    Lit un fichier de noeuds (format CovMulNet19) et renvoie :
        node_types : dict {node_id: type}
    On suppose une structure du style :
        node_id, type, description, ...
    """
    node_types = {}
    with open(node_file, "r", encoding="utf-8") as f:
        reader = csv.reader(f)
        header = next(reader)
        for row in reader:
            node_id = row[0]
            node_type = row[1].strip()
            node_types[node_id] = node_type
    vprint(f"[load_node_types] Loaded {len(node_types)} nodes from {node_file}")
    return node_types


def load_edges_csv(edge_file, weighted=False):
    """
    Lit un fichier d'arêtes CSV.
    Format attendu minimal : source, target [, weight, ...]
    Si weighted=True, on suppose que la 3e colonne est un poids numérique.
    """
    edges = []
    with open(edge_file, "r", encoding="utf-8") as f:
        reader = csv.reader(f)
        header = next(reader)
        for row in reader:
            src = row[0]
            tgt = row[1]
            if weighted and len(row) > 2:
                try:
                    w = float(row[2])
                except ValueError:
                    w = 1.0
                edges.append((src, tgt, w))
            else:
                edges.append((src, tgt))
    vprint(f"[load_edges_csv] Loaded {len(edges)} edges from {edge_file}")
    return edges


# -------------------------------------------------------------------------
# Core: degree / strength per mock
# -------------------------------------------------------------------------

def compute_metric_for_mock(
    node_types,
    edges,
    focal_type="GO",
    neighbor_type="protein",
    directed=True,
    weighted=False
):
    """
    Calcule pour un mock donné :
      - pour chaque nœud de type focal_type (ex: GO),
      - son degré (ou strength) vers des voisins de type neighbor_type (ex: protein).

    Retourne : dict {focal_node: valeur}

    Remarque :
    - weighted=True => utilise l'attribut 'weight' des arêtes (somme des poids)
    - weighted=False => simple degree (compte du nombre de voisins)
    """
    # Sélection des sets de noeuds
    focal_nodes = {n for n, t in node_types.items() if t.lower() == focal_type.lower()}
    neighbor_nodes = {n for n, t in node_types.items() if t.lower() == neighbor_type.lower()}

    if directed:
        G = nx.DiGraph()
    else:
        G = nx.Graph()

    # Ajout des noeuds & arêtes
    G.add_nodes_from(node_types.keys())
    if weighted:
        # edges: (src, tgt, weight)
        for src, tgt, w in edges:
            G.add_edge(src, tgt, weight=w)
    else:
        # edges: (src, tgt)
        G.add_edges_from(edges)

    metric = {n: 0.0 for n in focal_nodes}

    for node in focal_nodes:
        if directed:
            neigh_iter = G.successors(node)
        else:
            neigh_iter = G.neighbors(node)

        if weighted:
            total_weight = 0.0
            for neigh in neigh_iter:
                if neigh in neighbor_nodes:
                    data = G.get_edge_data(node, neigh)
                    if data is None:
                        continue
                    # si plusieurs arêtes, on les somme
                    if isinstance(data, dict) and "weight" in data:
                        total_weight += data["weight"]
                    elif isinstance(data, dict):
                        # multigraph-like, ou dict de dict
                        # on somme tous les 'weight' rencontrés
                        for sub in data.values():
                            if isinstance(sub, dict) and "weight" in sub:
                                total_weight += sub["weight"]
            metric[node] = total_weight
        else:
            count = 0
            for neigh in neigh_iter:
                if neigh in neighbor_nodes:
                    count += 1
            metric[node] = float(count)

    return metric


# -------------------------------------------------------------------------
# Distribution aggregator
# -------------------------------------------------------------------------

def aggregate_distributions_over_mocks(
    mock_basepath,
    node_template_file,
    focal_type="GO",
    neighbor_type="protein",
    directed=True,
    weighted=False,
    exclude_ids=None,
    verbose=True
):
    """
    Parcourt l'ensemble des mock networks dans mock_basepath
    (avec fichiers nodes_XXXX.csv / edges_XXXX.csv),
    calcule le degré/strength pour chaque noeud focal_type,
    et construit une distribution sur tous les mocks.

    Retourne :
        distributions : dict {node_id: [valeurs sur chaque mock]}
    """
    global VERBOSE
    VERBOSE = verbose

    # Types de nœuds (on les prend depuis un template, supposé commun à tous les mocks)
    node_types = load_node_types(node_template_file)

    # Liste des mocks
    mock_ids, nodes_dict, edges_dict = load_mock_network_ids(
        mock_basepath,
        exclude_ids=exclude_ids,
        verbose=verbose
    )

    distributions = defaultdict(list)

    vprint(f"[aggregate_distributions] Processing {len(mock_ids)} mock networks...")

    for mock_id in mock_ids:
        node_file = nodes_dict.get(mock_id, None)
        edge_file = edges_dict.get(mock_id, None)
        if node_file is None or edge_file is None:
            vprint(f"[aggregate_distributions][WARNING] Missing files for mock {mock_id}, skipping.")
            continue

        vprint(f"[aggregate_distributions] Mock {mock_id}...")

        edges = load_edges_csv(edge_file, weighted=weighted)

        metric = compute_metric_for_mock(
            node_types=node_types,
            edges=edges,
            focal_type=focal_type,
            neighbor_type=neighbor_type,
            directed=directed,
            weighted=weighted
        )

        for n, val in metric.items():
            distributions[n].append(val)

    vprint(f"[aggregate_distributions] Done. {len(distributions)} focal nodes with distributions.")
    return distributions


# -------------------------------------------------------------------------
# Export
# -------------------------------------------------------------------------

def export_distributions(distributions, outfile):
    """
    Exporte les distributions mock par nœud dans un TSV.

    Format de sortie :
        node_id, n_mocks, mean, sd, min, max, median, values

    Où "values" est une liste des valeurs séparées par des virgules.
    """
    rows = []
    for node, vals in distributions.items():
        arr = np.array(vals, dtype=float)
        n = len(arr)
        if n == 0:
            continue
        mean = float(np.mean(arr))
        # ddof=1 pour sd de type "sample", si n>1
        sd = float(np.std(arr, ddof=1)) if n > 1 else 0.0
        vmin = float(np.min(arr))
        vmax = float(np.max(arr))
        med = float(np.median(arr))
        values_str = ",".join(str(v) for v in vals)
        rows.append((node, n, mean, sd, vmin, vmax, med, values_str))

    # tri optionnel par mean décroissant (utile en debug)
    rows.sort(key=lambda r: r[2], reverse=True)

    os.makedirs(os.path.dirname(outfile), exist_ok=True)
    with open(outfile, "w", encoding="utf-8", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(["node_id", "n_mocks", "mean", "sd", "min", "max", "median", "values"])
        writer.writerows(rows)

    vprint(f"[export_distributions] Wrote {len(rows)} rows to {outfile}")


# -------------------------------------------------------------------------
# CLI
# -------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Compute mock degree distributions for focal nodes."
    )

    parser.add_argument("--mock_basepath", required=True,
                        help="Path to folder containing nodes_XXXX.csv and edges_XXXX.csv")

    parser.add_argument("--node_template_file", required=True,
                        help="Node file (CSV) with node types, used as template for mocks")

    parser.add_argument("--focal_type", required=True,
                        help="Node type whose distribution we compute (e.g., GO, disease, drug, symptom)")

    parser.add_argument("--neighbor_type", default="protein",
                        help="Neighbor node type used to define degree (default: protein)")

    parser.add_argument("--directed", default="True",
                        help="True/False for directed graphs")

    parser.add_argument("--weighted", default="False",
                        help="True/False: use edge weights if available (3rd column)")

    parser.add_argument("--exclude_ids", default="",
                        help="Comma-separated list of mock IDs to exclude")

    parser.add_argument("--outfile", required=True,
                        help="Output TSV file for distributions")

    parser.add_argument("--verbose", default="True",
                        help="True/False for verbose output")

    args = parser.parse_args()

    global VERBOSE
    VERBOSE = (args.verbose.lower() == "true")
    directed = (args.directed.lower() == "true")
    weighted = (args.weighted.lower() == "true")

    if args.exclude_ids.strip():
        exclude_ids = [int(x) for x in args.exclude_ids.split(",") if x.strip().isdigit()]
    else:
        exclude_ids = None

    dists = aggregate_distributions_over_mocks(
        mock_basepath=args.mock_basepath,
        node_template_file=args.node_template_file,
        focal_type=args.focal_type,
        neighbor_type=args.neighbor_type,
        directed=directed,
        weighted=weighted,
        exclude_ids=exclude_ids,
        verbose=VERBOSE
    )

    export_distributions(dists, args.outfile)


if __name__ == "__main__":
    main()
