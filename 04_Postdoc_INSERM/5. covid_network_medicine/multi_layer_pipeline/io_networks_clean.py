#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Input: multilayer network + mock networks provided by M. De Domenico (CoMuNe Lab)

This module handles:
- discovery of mock network files (nodes_*.csv, edges_*.csv)
- consistency between node files and edge files
- filtering out corrupted or unwanted networks
- returning a clean, sorted list of available mock network IDs

Used in:
- compute_mock_distributions.py
- downstream bootstrap and z-score computations
"""

import os
import re

VERBOSE = True

def vprint(*args, **kwargs):
    if VERBOSE:
        print(*args, **kwargs)


# ----------------------------------------------------------------------
# Utilities
# ----------------------------------------------------------------------

def _extract_id(filename, prefix="nodes_"):
    """
    Extract numeric ID from filenames like: nodes_00342.csv
    """
    num = filename.replace(".csv", "").replace(prefix, "")
    if num.isdigit():
        return int(num)
    return None


def _detect_files(basepath):
    """
    Returns two dicts:
    - nodes: {id: filepath}
    - edges: {id: filepath}
    """
    nodes = {}
    edges = {}

    for fname in os.listdir(basepath):
        if fname.startswith("nodes_") and fname.endswith(".csv"):
            node_id = _extract_id(fname, prefix="nodes_")
            if node_id is not None:
                nodes[node_id] = os.path.join(basepath, fname)

        elif fname.startswith("edges_") and fname.endswith(".csv"):
            edge_id = _extract_id(fname, prefix="edges_")
            if edge_id is not None:
                edges[edge_id] = os.path.join(basepath, fname)

    return nodes, edges


def _detect_inconsistencies(nodes_dict, edges_dict):
    """
    Detect IDs present in one set and missing in the other.
    """
    node_ids = set(nodes_dict.keys())
    edge_ids = set(edges_dict.keys())
    diff_nodes = node_ids - edge_ids
    diff_edges = edge_ids - node_ids
    return diff_nodes, diff_edges


# ----------------------------------------------------------------------
# Main function
# ----------------------------------------------------------------------

def load_mock_network_ids(
    basepath,
    exclude_ids=None,
    verbose=True,
):
    """
    Detects all usable mock network IDs in a folder.

    Parameters
    ----------
    basepath : str
        Path containing nodes_*.csv and edges_*.csv.
    exclude_ids : list[int]
        List of network IDs to exclude (corrupted or invalid).
    verbose : bool

    Returns
    -------
    mock_ids : list[int]
        Sorted list of valid mock network IDs.
    nodes_dict : dict
        Mapping id -> node file path
    edges_dict : dict
        Mapping id -> edge file path
    """

    global VERBOSE
    VERBOSE = verbose

    vprint(f"[io_networks] Scanning folder: {basepath}")

    nodes_dict, edges_dict = _detect_files(basepath)
    vprint(f"[io_networks] Found {len(nodes_dict)} node files.")
    vprint(f"[io_networks] Found {len(edges_dict)} edge files.")

    # Check consistency
    diff_nodes, diff_edges = _detect_inconsistencies(nodes_dict, edges_dict)
    if diff_nodes:
        vprint(f"[io_networks][WARNING] IDs missing edges: {sorted(diff_nodes)}")
    if diff_edges:
        vprint(f"[io_networks][WARNING] IDs missing nodes: {sorted(diff_edges)}")

    # Intersection = usable
    usable_ids = set(nodes_dict.keys()).intersection(edges_dict.keys())
    vprint(f"[io_networks] IDs with both nodes & edges: {len(usable_ids)}")

    # Exclusions
    if exclude_ids:
        vprint(f"[io_networks] Excluding {len(exclude_ids)} IDs explicitly.")
        usable_ids = usable_ids - set(exclude_ids)

    mock_ids = sorted(list(usable_ids))
    vprint(f"[io_networks] Final usable mock networks: {len(mock_ids)}")

    return mock_ids, nodes_dict, edges_dict


# ----------------------------------------------------------------------
# Example usage (not executed when imported)
# ----------------------------------------------------------------------

if __name__ == "__main__":
    # Example provided path
    base = "A:/Downloads/Projects/workFromHome/Projects/Covid/random_networks/Mock_networks_21Apr"

    # Example known-bad networks (from your scripts)
    exclude = [494, 525, 532, 569, 592, 629, 658, 662, 674, 676, 
               701, 705, 712, 736, 757, 788, 796, 825, 827, 838, 
               873, 876, 908, 910, 916, 917, 918, 938, 958, 965, 
               1088, 1125, 1128, 1145, 1146, 1165, 1171, 1183, 1185, 
               1192, 1195, 1211, 1244, 1254, 1257, 1300, 1353, 1354, 
               1364, 1372, 1394, 1395, 1401, 1410, 1427, 1454, 1477, 
               1497, 1517, 1532, 1549, 1551, 1575, 1578, 1611, 1624, 
               1629, 1637, 1650, 1655, 1697, 1700, 1713, 1736, 1745, 
               1762, 1790, 1800, 1815, 1824, 1827, 1866, 1872, 1886, 
               1887, 2053, 2072, 2106, 2115, 2128, 2160, 2217, 2343, 
               2383, 2388, 2413, 2458, 2465]

    ids, nd, ed = load_mock_network_ids(base, exclude_ids=exclude, verbose=True)
    print(ids[:10])  # first 10 usable ids
