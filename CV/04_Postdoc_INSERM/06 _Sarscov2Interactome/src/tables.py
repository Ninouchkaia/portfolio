from __future__ import annotations
from pathlib import Path
from typing import Dict, List
import logging

import pandas as pd

from .config import INTERACTORS_DIR, GENE_LISTS_DIR, TABLES_DIR
from .io_utils import list_virus_dirs

logger = logging.getLogger(__name__)


def _read_interactor_file(path: Path) -> List[str]:
    """
    Lit un fichier d'interacteurs au format :
      1: nom du virus
      2-5: éventuellement metadata/commentaires
      6+: "<protein_id>\t<gene_symbol>"
    Retourne la liste de gene_symbol.
    """
    genes: List[str] = []
    with path.open() as f:
        lines = f.readlines()
    for line in lines[1:]:
        if line.startswith("#") or not line.strip():
            continue
        parts = line.strip().split("\t")
        if len(parts) < 2:
            continue
        gene = parts[1]
        genes.append(gene)
    return genes


def build_gene_lists_from_interactors(
    range_mode: str = "direct_and_second",
) -> Path:
    """
    Construit un fichier genes_list_*.txt au format attendu par tes scripts R.
    range_mode peut être :
      - "direct"                  → direct_interactors.txt
      - "only_second"             → only_second_range_interactors.txt
      - "direct_and_second"       → direct_interactors + only_second_range
      - "orders_1_to_3"           → union des ordres 1,2,3
      - etc.
    """
    pattern_by_mode = {
        "direct": ["direct_interactors.txt"],
        "only_second": ["only_second_range_interactors.txt"],
        "only_third": ["only_Third_range_interactors.txt"],
        "only_fourth": ["only_Fourth_range_interactors.txt"],
        "orders_1_to_3": [
            "direct_interactors.txt",
            "only_second_range_interactors.txt",
            "only_Third_range_interactors.txt",
        ],
        "direct_and_second": [
            "direct_interactors.txt",
            "only_second_range_interactors.txt",
        ],
    }

    if range_mode not in pattern_by_mode:
        raise ValueError(f"Unknown range_mode: {range_mode}")

    logger.info("Building gene list for range_mode=%s", range_mode)

    gene_lists: Dict[str, List[str]] = {}

    for virus_dir in list_virus_dirs(INTERACTORS_DIR):
        virus_name = virus_dir.name
        genes_set: set[str] = set()
        for filename in pattern_by_mode[range_mode]:
            fpath = virus_dir / filename
            if not fpath.exists():
                continue
            genes = _read_interactor_file(fpath)
            genes_set.update(genes)
        gene_lists[virus_name] = sorted(genes_set)

    out_path = GENE_LISTS_DIR / f"genes_list_{range_mode}.txt"
    with out_path.open("w", encoding="utf-8") as f:
        for virus, genes in sorted(gene_lists.items()):
            if not genes:
                continue
            line = "\t".join([virus] + genes)
            f.write(line + "\n")

    logger.info("Wrote gene list file: %s", out_path)
    return out_path


def build_gene_virus_table(range_mode: str = "direct_and_second") -> Path:
    """
    Construit une matrice binaire (genes x virus) pour un range donné,
    à partir des fichiers dans INTERACTORS_DIR.
    """
    pattern_by_mode = {
        "direct": ["direct_interactors.txt"],
        "direct_and_second": [
            "direct_interactors.txt",
            "only_second_range_interactors.txt",
        ],
        "only_second": ["only_second_range_interactors.txt"],
        "only_third": ["only_Third_range_interactors.txt"],
    }

    if range_mode not in pattern_by_mode:
        raise ValueError(f"Unknown range_mode: {range_mode}")

    virus_to_genes: Dict[str, set[str]] = {}
    all_genes: set[str] = set()

    for virus_dir in list_virus_dirs(INTERACTORS_DIR):
        virus_name = virus_dir.name
        genes_set: set[str] = set()
        for filename in pattern_by_mode[range_mode]:
            fpath = virus_dir / filename
            if fpath.exists():
                genes_set.update(_read_interactor_file(fpath))
        if genes_set:
            virus_to_genes[virus_name] = genes_set
            all_genes.update(genes_set)

    df = pd.DataFrame(
        0,
        index=sorted(all_genes),
        columns=sorted(virus_to_genes.keys()),
        dtype=int,
    )
    for virus, genes in virus_to_genes.items():
        df.loc[list(genes), virus] = 1

    out_path = TABLES_DIR / f"gene_virus_table_{range_mode}.tsv"
    df.to_csv(out_path, sep="\t")
    logger.info("Wrote gene-virus table: %s", out_path)
    return out_path
