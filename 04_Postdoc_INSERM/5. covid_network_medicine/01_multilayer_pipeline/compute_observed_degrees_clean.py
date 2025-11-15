#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Input: multilayer network + mock networks provided by M. De Domenico (CoMuNe Lab)

This module computes the OBSERVED structural metrics of CovMulNet19:
- degree / strength for each node of a given source type to a given target type
- supports directed or undirected mode
- exports:
    * observed degrees
    * observed relative degrees
    * observed withZeros tables
    * rankings

This clean version replaces the following historical scripts:
- protein-protein-degrees-rankings.py
- protein-drug-degrees-(directed|undirected)-rankings.py
- protein-disease-degrees-directed-rankings-*.py
- protein-symptom-degrees-directed-rankings-*.py
- protein-GO-degrees-directed-rankings-*.py
- protein-target-degrees-directed/undirected corrected generic scripts
- all “withZeros.tsv” generation blocks
- all degree ranking scripts

Usage (CLI example):
    python compute_observed_degrees_clean.py \
        --nodes COVID19_GDDS_nodes.csv \
        --edges COVID19_GDDS_edges.csv \
        --source_type protein \
        --target_type disease \
        --directed True \
        --output outdir/
"""

import argparse
import csv
import os
import networkx as nx


# -------------------------------------------------------------------------
# Verbose print
# -------------------------------------------------------------------------

VERBOSE = True

def vprint(*a, **k):
    if VERBOSE:
        print(*a, **k)


# -------------------------------------------------------------------------
# I/O Helpers
# -------------------------------------------------------------------------

def load_nodes(node_file):
    """
    Load node list from CovMulNet19 node file.
    Expected format: node_id, node_type, description, ...
    """
    nodes = {}
    with open(node_file, "r", encoding="utf-8") as f:
        reader = csv.reader(f)
        header = next(reader)
        for row in reader:
            node_id = row[0]
            node_type = row[1].strip()
            nodes[node_id] = node_type
    vprint(f"[load_nodes] Loaded {len(nodes)} nodes.")
    return nodes


def load_edges(edge_file):
    """
    Load edges from CovMulNet19 edge file.
    Expected format: source, target, ...
    """
    edges = []
    with open(edge_file, "r", encoding="utf-8") as f:
        reader = csv.reader(f)
        header = next(reader)
        for row in reader:
            src = row[0]
            tgt = row[1]
            edges.append((src, tgt))
    vprint(f"[load_edges] Loaded {len(edges)} edges.")
    return edges


# -------------------------------------------------------------------------
# Core Computation
# -------------------------------------------------------------------------

def compute_observed_degrees(
    nodes,
    edges,
    source_type="protein",
    target_type="disease",
    directed=True
):
    """
    Build graph, filter by source-target types, return degree mapping.

    Returns:
        observed: {source_node: count of edges to target nodes}
        all_source_nodes: list of all nodes of source type
        all_target_nodes: list of nodes of target type
    """
    # Select nodes of interest
    source_nodes = {n for n, t in nodes.items() if t.lower() == source_type.lower()}
    target_nodes = {n for n, t in nodes.items() if t.lower() == target_type.lower()}

    vprint(f"[compute_observed_degrees] {len(source_nodes)} source nodes ({source_type})")
    vprint(f"[compute_observed_degrees] {len(target_nodes)} target nodes ({target_type})")

    # Build graph
    if directed:
        G = nx.DiGraph()
    else:
        G = nx.Graph()

    G.add_nodes_from(nodes.keys())
    G.add_edges_from(edges)

    observed = {}

    # Count connections source -> target
    for s in source_nodes:
        count = 0
        if directed:
            neighbors = G.successors(s)
        else:
            neighbors = G.neighbors(s)

        for neigh in neighbors:
            if neigh in target_nodes:
                count += 1

        observed[s] = count

    return observed, sorted(source_nodes), sorted(target_nodes)


# -------------------------------------------------------------------------
# Helpers: relative degrees, zeros, ranking
# -------------------------------------------------------------------------

def compute_relative_degrees(observed_dict):
    total = sum(observed_dict.values())
    if total == 0:
        return {k: 0.0 for k in observed_dict}
    return {k: v / total for k, v in observed_dict.items()}


def write_table(path, header, rows):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", encoding="utf-8", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(header)
        writer.writerows(rows)
    vprint(f"[write] Wrote: {path}")


# -------------------------------------------------------------------------
# Export Full Results (observed, relative, withZeros, ranking)
# -------------------------------------------------------------------------

def export_observed_results(
    observed,
    source_nodes,
    source_type,
    target_type,
    outdir,
    directed
):
    mode = "directed" if directed else "undirected"

    # 1. observed degrees (only nodes w/ deg > 0)
    rows_obs = [(n, observed[n]) for n in source_nodes if observed[n] > 0]
    write_table(
        os.path.join(outdir, f"{source_type}_to_{target_type}_{mode}_observed.tsv"),
        ["source_node", "degree"],
        rows_obs
    )

    # 2. observed relative degrees
    rel = compute_relative_degrees(observed)
    rows_rel = [(n, rel[n]) for n in source_nodes if observed[n] > 0]
    write_table(
        os.path.join(outdir, f"{source_type}_to_{target_type}_{mode}_relative.tsv"),
        ["source_node", "relative_degree"],
        rows_rel
    )

    # 3. withZeros (all source nodes)
    rows_zero = [(n, observed[n], rel[n]) for n in source_nodes]
    write_table(
        os.path.join(outdir, f"{source_type}_to_{target_type}_{mode}_withZeros.tsv"),
        ["source_node", "degree", "relative_degree"],
        rows_zero
    )

    # 4. ranking (sorted descending)
    ranked = sorted(rows_obs, key=lambda x: x[1], reverse=True)
    write_table(
        os.path.join(outdir, f"{source_type}_to_{target_type}_{mode}_ranking.tsv"),
        ["source_node", "degree"],
        ranked
    )


# -------------------------------------------------------------------------
# Main (CLI)
# -------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser()

    parser.add_argument("--nodes", required=True, help="COVID19 nodes CSV")
    parser.add_argument("--edges", required=True, help="COVID19 edges CSV")

    parser.add_argument("--source_type", required=True,
                        help="Type of source nodes (e.g., protein)")
    parser.add_argument("--target_type", required=True,
                        help="Type of target nodes (e.g., disease)")

    parser.add_argument("--directed", default="True",
                        help="True/False")

    parser.add_argument("--output", required=True,
                        help="Output directory")

    parser.add_argument("--verbose", default="True",
                        help="True/False")

    args = parser.parse_args()

    global VERBOSE
    VERBOSE = (args.verbose.lower() == "true")
    directed = (args.directed.lower() == "true")

    nodes = load_nodes(args.nodes)
    edges = load_edges(args.edges)

    observed, source_nodes, target_nodes = compute_observed_degrees(
        nodes,
        edges,
        source_type=args.source_type,
        target_type=args.target_type,
        directed=directed
    )

    export_observed_results(
        observed,
        source_nodes,
        args.source_type,
        args.target_type,
        args.output,
        directed
    )


if __name__ == "__main__":
    main()
