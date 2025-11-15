
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import csv
import re
from collections import Counter
import networkx as nx
import pandas as pd

VERBOSE = False
def vprint(msg):
    if VERBOSE:
        print(msg)

def load_nodes(path):
    with open(path, 'r') as f:
        reader = csv.reader(f)
        data = [n for n in reader][1:]
    names = [n[0] for n in data]
    type_dict = {n[0]: n[1] for n in data}
    descr_dict = {n[0]: n[2] for n in data}
    vprint(f"Loaded {len(names)} nodes")
    return names, type_dict, descr_dict

def load_edges(path):
    with open(path, 'r') as f:
        reader = csv.reader(f)
        edges = [tuple(e) for e in reader][1:]
    vprint(f"Loaded {len(edges)} edges")
    return edges

def build_graph(nodes, edges, type_dict, descr_dict):
    G = nx.Graph()
    G.add_nodes_from(nodes)
    G.add_edges_from(edges)
    nx.set_node_attributes(G, type_dict, 'type')
    nx.set_node_attributes(G, descr_dict, 'description')
    vprint(f"Graph: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges")
    return G

def parse_drugs(drug_list):
    parsed = []
    for d in drug_list:
        if " AND " in d:
            parsed.extend([x.strip() for x in d.split(" AND ")])
        if "; " in d:
            parsed.extend([x.strip() for x in d.split("; ")])
        parsed.append(d.strip())
    return list(dict.fromkeys(parsed))

def remove_parentheses(drugs):
    cleaned = []
    for d in drugs:
        if re.search(r'\s\([^()]*\)', d):
            d = re.sub(r'\s\([^()]*\)', '', d)
            cleaned.append(d)
    return list(dict.fromkeys(cleaned))

def pipeline(nodes_csv, edges_csv, drugs_csv, export_csv):
    nodes, type_dict, descr_dict = load_nodes(nodes_csv)
    edges = load_edges(edges_csv)
    G = build_graph(nodes, edges, type_dict, descr_dict)

    drug_sub = [k for k,v in descr_dict.items() if v.upper()=="DRUG"]
    parsed = parse_drugs(drug_sub)
    combined = list(dict.fromkeys(drug_sub + parsed))

    no_par = remove_parentheses(combined)
    final = list(dict.fromkeys(combined + no_par))

    df = pd.DataFrame(final, columns=["drug"])
    df.to_csv(export_csv, index=False)
    vprint(f"Exported {len(final)} entries to {export_csv}")

if __name__=="__main__":
    import argparse
    p = argparse.ArgumentParser()
    p.add_argument("--nodes")
    p.add_argument("--edges")
    p.add_argument("--drugs")
    p.add_argument("--out", default="final_drugs.csv")
    p.add_argument("--verbose", action="store_true")
    args = p.parse_args()

    VERBOSE = args.verbose
    pipeline(args.nodes, args.edges, args.drugs, args.out)
