#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Create full project structure for covid_networks pipeline.
Generates ALL folders and ALL files with their full content.
"""

from pathlib import Path


# ----------------------------------------------------------------------------------------
# Utility
# ----------------------------------------------------------------------------------------

def write(path: Path, content: str):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content, encoding="utf-8")


# ----------------------------------------------------------------------------------------
# Project root
# ----------------------------------------------------------------------------------------

ROOT = Path("covid_networks")

# ----------------------------------------------------------------------------------------
# CONTENTS OF FILES
# ----------------------------------------------------------------------------------------

# -------------------------- src/config.py --------------------------
config_py = r'''
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]

DATA_RAW = PROJECT_ROOT / "data" / "raw" / "Virus_host_interactomes_thresh25" / "thresh_0.25"
DATA_INTERMEDIATE = PROJECT_ROOT / "data" / "intermediate"
DATA_RESULTS = PROJECT_ROOT / "data" / "results"
R_SCRIPTS_DIR = PROJECT_ROOT / "r"

INTERACTORS_DIR = DATA_INTERMEDIATE / "interactors"
GENE_LISTS_DIR = DATA_INTERMEDIATE / "gene_lists"
ENRICHMENT_DIR = DATA_INTERMEDIATE / "enrichment"

TABLES_DIR = DATA_RESULTS / "tables"
FIGURES_DIR = DATA_RESULTS / "figures"

for p in [
    DATA_RAW,
    DATA_INTERMEDIATE,
    DATA_RESULTS,
    INTERACTORS_DIR,
    GENE_LISTS_DIR,
    ENRICHMENT_DIR,
    TABLES_DIR,
    FIGURES_DIR,
]:
    p.mkdir(parents=True, exist_ok=True)
'''

# -------------------------- src/logging_utils.py --------------------------
logging_utils_py = r'''
import logging
from pathlib import Path

def setup_logging(log_file: Path | None = None) -> None:
    log_format = "%(asctime)s - %(levelname)s - %(name)s - %(message)s"
    logging.basicConfig(
        level=logging.INFO,
        format=log_format,
        handlers=[logging.StreamHandler()]
    )
    if log_file is not None:
        file_handler = logging.FileHandler(log_file)
        file_handler.setFormatter(logging.Formatter(log_format))
        logging.getLogger().addHandler(file_handler)
'''

# -------------------------- src/io_utils.py --------------------------
io_utils_py = r'''
from pathlib import Path
from typing import Dict, List, Tuple
import csv

def list_virus_dirs(rootdir: Path) -> List[Path]:
    return [
        d for d in rootdir.iterdir()
        if d.is_dir()
    ]

def load_nodes_edges(
    virus_dir: Path,
    nodes_filename: str = "nodes.csv",
    edges_filename: str = "edges.csv",
    delimiter_nodes: str = ",",
    delimiter_edges: str = ",",
    has_header: bool = False,
):
    nodes_path = virus_dir / nodes_filename
    edges_path = virus_dir / edges_filename

    descr_dict_nodes = {}
    host_nodes = {}
    virus_nodes = {}

    with nodes_path.open() as f:
        reader = csv.reader(f, delimiter=delimiter_nodes)
        if has_header:
            next(reader)
        for row in reader:
            if not row:
                continue
            node_id = row[0]
            node_name = row[1]
            node_type = row[2].strip()
            descr_dict_nodes[node_id] = node_name
            if node_type in {"1", "1.0"}:
                host_nodes[node_id] = node_name
            elif node_type in {"0", "0.0"}:
                virus_nodes[node_id] = node_name

    edges = []
    with edges_path.open() as f:
        reader = csv.reader(f, delimiter=delimiter_edges)
        if has_header:
            next(reader)
        for row in reader:
            if len(row) >= 2:
                edges.append((row[0], row[1]))

    return descr_dict_nodes, host_nodes, virus_nodes, edges
'''

# -------------------------- src/network.py --------------------------
network_py = r'''
import networkx as nx
from typing import Iterable, Dict, Set

def build_graph(node_ids: Iterable[str], edges: Iterable[tuple[str, str]]) -> nx.Graph:
    G = nx.Graph()
    G.add_nodes_from(node_ids)
    G.add_edges_from(edges)
    return G

def compute_order_interactors(
    G: nx.Graph,
    virus_nodes: Iterable[str],
    host_nodes: set[str],
    max_order: int,
) -> Dict[int, Set[str]]:
    virus_nodes = set(virus_nodes)
    visited = set(virus_nodes)
    frontier = set(virus_nodes)
    results = {}

    for order in range(1, max_order + 1):
        next_frontier = set()
        for u in frontier:
            for v in G.neighbors(u):
                if v not in visited:
                    visited.add(v)
                    next_frontier.add(v)
        results[order] = {n for n in next_frontier if n in host_nodes}
        frontier = next_frontier
    return results
'''

# -------------------------- src/interactors.py --------------------------
interactors_py = r'''
from pathlib import Path
import logging
from .config import DATA_RAW, INTERACTORS_DIR
from .io_utils import list_virus_dirs, load_nodes_edges
from .network import build_graph, compute_order_interactors

logger = logging.getLogger(__name__)

def compute_all_orders_for_all_viruses(
    max_order: int = 2,
    nodes_filename: str = "nodes.csv",
    edges_filename: str = "edges.csv",
) -> None:

    virus_dirs = list_virus_dirs(DATA_RAW)
    logger.info("Found %d virus folders", len(virus_dirs))

    for virus_dir in virus_dirs:
        virus_name = virus_dir.name
        logger.info("Processing virus %s", virus_name)

        descr, host_nodes, virus_nodes, edges = load_nodes_edges(
            virus_dir,
            nodes_filename=nodes_filename,
            edges_filename=edges_filename,
            delimiter_nodes=",",
            delimiter_edges=",",
            has_header=False,
        )

        G = build_graph(descr.keys(), edges)
        logger.info("Graph %s: %d nodes, %d edges", virus_name, G.number_of_nodes(), G.number_of_edges())

        order_to_nodes = compute_order_interactors(
            G,
            virus_nodes=virus_nodes.keys(),
            host_nodes=set(host_nodes.keys()),
            max_order=max_order,
        )

        out_dir = INTERACTORS_DIR / virus_name
        out_dir.mkdir(parents=True, exist_ok=True)

        write_file(order_to_nodes, virus_name, out_dir, descr)


def write_file(order_to_nodes, virus_name, out_dir, descr):
    mapping = {
        1: "direct_interactors.txt",
        2: "only_second_range_interactors.txt",
        3: "only_Third_range_interactors.txt",
        4: "only_Fourth_range_interactors.txt",
        5: "only_Fifth_range_interactors.txt",
    }

    for order, nodes in order_to_nodes.items():
        fname = mapping.get(order, f"only_{order}_range_interactors.txt")
        path = out_dir / fname
        with path.open("w", encoding="utf-8") as f:
            f.write(f"{virus_name}\n")
            for pid in sorted(nodes):
                gene = descr.get(pid, "NA")
                f.write(f"{pid}\t{gene}\n")
        logger.info("Wrote %d entries to %s", len(nodes), path)
'''

# -------------------------- src/tables.py --------------------------
tables_py = r'''
from pathlib import Path
import logging
import pandas as pd
from .config import INTERACTORS_DIR, GENE_LISTS_DIR, TABLES_DIR
from .io_utils import list_virus_dirs

logger = logging.getLogger(__name__)

def _read_genes(path: Path):
    genes = []
    with path.open() as f:
        lines = f.readlines()
    for line in lines:
        if line.startswith("#") or not line.strip():
            continue
        if "\t" in line:
            parts = line.strip().split("\t")
            if len(parts) > 1:
                genes.append(parts[1])
    return genes

def build_gene_lists_from_interactors(range_mode="direct_and_second"):
    pattern = {
        "direct": ["direct_interactors.txt"],
        "direct_and_second": ["direct_interactors.txt", "only_second_range_interactors.txt"],
        "only_second": ["only_second_range_interactors.txt"],
        "only_third": ["only_Third_range_interactors.txt"],
        "orders_1_to_3": [
            "direct_interactors.txt",
            "only_second_range_interactors.txt",
            "only_Third_range_interactors.txt",
        ],
    }

    files = pattern[range_mode]
    virus_to_genes = {}

    for virus_dir in list_virus_dirs(INTERACTORS_DIR):
        virus = virus_dir.name
        genes = set()
        for fname in files:
            path = virus_dir / fname
            if path.exists():
                genes.update(_read_genes(path))
        if genes:
            virus_to_genes[virus] = sorted(genes)

    out = GENE_LISTS_DIR / f"genes_list_{range_mode}.txt"
    with out.open("w") as f:
        for virus, genes in sorted(virus_to_genes.items()):
            f.write(virus + "\t" + "\t".join(genes) + "\n")

    logger.info("Wrote %s", out)
    return out

def build_gene_virus_table(range_mode="direct_and_second"):
    pattern = {
        "direct": ["direct_interactors.txt"],
        "direct_and_second": ["direct_interactors.txt", "only_second_range_interactors.txt"],
        "only_second": ["only_second_range_interactors.txt"],
        "only_third": ["only_Third_range_interactors.txt"],
    }

    files = pattern[range_mode]
    virus_to_genes = {}
    all_genes = set()

    for virus_dir in list_virus_dirs(INTERACTORS_DIR):
        virus = virus_dir.name
        genes = set()
        for fname in files:
            path = virus_dir / fname
            if path.exists():
                genes.update(_read_genes(path))
        if genes:
            virus_to_genes[virus] = genes
            all_genes.update(genes)

    df = pd.DataFrame(
        0,
        index=sorted(all_genes),
        columns=sorted(virus_to_genes.keys()),
        dtype=int,
    )
    for virus, genes in virus_to_genes.items():
        df.loc[list(genes), virus] = 1

    out = TABLES_DIR / f"gene_virus_table_{range_mode}.tsv"
    df.to_csv(out, sep="\t")
    logger.info("Wrote %s", out)
    return out
'''

# -------------------------- src/enrichment_wrapper.py --------------------------
enrichment_wrapper_py = r'''
from pathlib import Path
import subprocess
import logging
from .config import R_SCRIPTS_DIR, ENRICHMENT_DIR

logger = logging.getLogger(__name__)

def run_r_enrichment(
    gene_list_file: Path,
    script_name: str = "enrich_reactome_compareCluster.R",
    output_prefix: str = "enrichPathway",
):
    script = R_SCRIPTS_DIR / script_name
    out = ENRICHMENT_DIR / f"{output_prefix}.tsv"

    cmd = ["Rscript", str(script), str(gene_list_file), str(out)]
    logger.info("Running: %s", " ".join(cmd))
    subprocess.run(cmd, check=True)

    return out
'''

# -------------------------- r/enrich_reactome_compareCluster.R --------------------------
enrich_reactome_R = r'''
#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("Usage: script <genes_list> <out.tsv>")

genes_list <- args[1]
out_tsv <- args[2]

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(ReactomePA)
  library(reactome.db)
  library(readr)
})

gcRaw <- read.table(
  genes_list,
  header = FALSE,
  sep = "\t",
  col.names = paste0("V", seq_len(max(count.fields(genes_list, sep = "\t")))),
  fill = TRUE,
  stringsAsFactors = FALSE
)

df <- data.frame(t(gcRaw[-1]))
colnames(df) <- gcRaw[, 1]

# Ici tu ajoutes ton mapping gene->entrez si besoin.

ck <- compareCluster(
  geneCluster = df,
  fun = "enrichPathway",
  pvalueCutoff = 0.005
)

write_tsv(as.data.frame(ck, stringsAsFactors = FALSE), out_tsv)
'''

# -------------------------- README.md --------------------------
readme_md = r'''
# SARS-CoV-2 – Host Interactome Pipeline (Strict Reproducible Version)

This repository contains a fully structured, reproducible pipeline for computing
viral–host interactors and pathway enrichments using Python + R.

## Steps

### 1. Compute interactors (NetworkX)
    python analysis.py compute_interactors --max-order 4

### 2. Build gene lists for enrichment
    python analysis.py gene_lists --range-mode direct_and_second

### 3. Enrichment (Reactome via Rscript)
    python analysis.py enrich_reactome --range-mode direct_and_second

### 4. Build gene–virus matrices
    python analysis.py gene_virus_table --range-mode direct_and_second

All results are written to:
- data/intermediate/*
- data/results/*
'''

