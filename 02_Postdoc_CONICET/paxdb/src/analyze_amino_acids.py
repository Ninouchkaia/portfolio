# paxdb/scripts/analyze_amino_acids.py

from __future__ import annotations

"""
Main entry point for the PaxDB amino acid usage pipeline.

Typical usage:

    python scripts/analyze_amino_acids.py \
        --species-metadata data/metadata/species.tsv \
        --aa-costs data/metadata/amino_acid_costs.tsv \
        --outdir results

This script roughly corresponds to the workflow that used to be
spread across multiple files like:
- amino_acid_count*.py
- newdef_protein*.py
- aa_relationship.py
- protein_relationship.py
but in a cleaner, modular form.
"""

import argparse
import logging
import sys
from pathlib import Path
from typing import List

import pandas as pd

# Make src/ importable when running this script directly
THIS_DIR = Path(__file__).resolve().parent
SRC_DIR = THIS_DIR.parent / "src"
sys.path.insert(0, str(SRC_DIR))

from utils import setup_logging, resolve_path  # type: ignore
from fasta_parser import load_fasta_as_dict  # type: ignore
from abundance_loader import (  # type: ignore
    load_abundance_table,
    load_mapping_table,
    map_abundances_to_fasta_ids,
)
from protein import Protein  # type: ignore
from aa_metrics import (  # type: ignore
    compute_proteome_aa_usage,
    load_amino_acid_costs,
    STANDARD_AA,
)
from relationships import (  # type: ignore
    correlate_usage_with_cost,
    results_to_dataframe,
)

logger = logging.getLogger(__name__)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Analyze amino acid usage under metabolic constraints (PaxDB)."
    )
    parser.add_argument(
        "--species-metadata",
        type=str,
        required=True,
        help="TSV file describing species and paths to FASTA/abundance/mapping.",
    )
    parser.add_argument(
        "--aa-costs",
        type=str,
        required=True,
        help="TSV file with amino acid costs (columns: aa, cost_atp).",
    )
    parser.add_argument(
        "--outdir",
        type=str,
        default="results",
        help="Output directory for results (default: results).",
    )
    parser.add_argument(
        "--logdir",
        type=str,
        default="results/logs",
        help="Directory for log file (default: results/logs).",
    )
    return parser.parse_args()


def load_species_metadata(path: Path) -> pd.DataFrame:
    """
    Load species metadata describing where to find inputs.

    Expected columns:
    - species_id
    - fasta_path          (relative to this file's directory)
    - abundance_path
    - mapping_path        (can be empty)
    """
    df = pd.read_csv(path, sep="\t", comment="#", dtype=str)
    required = ["species_id", "fasta_path", "abundance_path", "mapping_path"]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"Missing columns in species metadata: {missing}")
    return df


def build_protein_objects_for_species(
    species_id: str,
    fasta_path: Path,
    abundance_path: Path,
    mapping_path: Path | None,
) -> List[Protein]:
    """
    For one species, build Protein objects with sequence and abundance.

    Steps:
    1. Load FASTA (sequence_id → sequence).
    2. Load abundance table.
    3. Optional: load mapping table and map to FASTA ids.
    4. Build Protein instances for each mapped entry.
    """
    logger.info("=== Species %s ===", species_id)
    logger.info("FASTA: %s", fasta_path)
    logger.info("Abundance: %s", abundance_path)
    if mapping_path is not None:
        logger.info("Mapping: %s", mapping_path)

    seq_dict = load_fasta_as_dict(fasta_path)
    logger.info("Loaded %d protein sequences from FASTA.", len(seq_dict))

    abund_df = load_abundance_table(abundance_path)

    mapping_df = (
        load_mapping_table(mapping_path) if mapping_path is not None else None
    )

    mapped_abund = map_abundances_to_fasta_ids(
        abund_df=abund_df,
        fasta_ids=set(seq_dict.keys()),
        mapping_df=mapping_df,
    )

    proteins: List[Protein] = []
    for seq_id, abundance in zip(mapped_abund["seq_id"], mapped_abund["abundance"]):
        seq = seq_dict.get(seq_id)
        if not seq:
            continue
        proteins.append(Protein(seq_id=seq_id, sequence=seq, abundance=float(abundance)))

    logger.info("Built %d Protein objects for species %s", len(proteins), species_id)
    return proteins


def main() -> None:
    args = parse_args()

    outdir = Path(args.outdir).resolve()
    logdir = Path(args.logdir).resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    setup_logging(logdir)
    logger.info("Starting PaxDB amino acid analysis pipeline.")
    logger.info("Output directory: %s", outdir)

    species_meta_path = Path(args.species_metadata).resolve()
    aa_costs_path = Path(args.aa_costs).resolve()

    species_df = load_species_metadata(species_meta_path)
    aa_costs = load_amino_acid_costs(aa_costs_path)

    # DataFrames to accumulate AA usage across species
    usage_unweighted_rows = []
    usage_weighted_rows = []
    relationship_results = []

    meta_base = species_meta_path.parent

    for _, row in species_df.iterrows():
        species_id = row["species_id"]

        fasta_path = resolve_path(meta_base, row["fasta_path"])
        abundance_path = resolve_path(meta_base, row["abundance_path"])
        mapping_path = row.get("mapping_path") or ""
        mapping_path_resolved = (
            resolve_path(meta_base, mapping_path) if mapping_path.strip() else None
        )

        proteins = build_protein_objects_for_species(
            species_id=species_id,
            fasta_path=fasta_path,
            abundance_path=abundance_path,
            mapping_path=mapping_path_resolved,
        )

        if not proteins:
            logger.warning(
                "No Protein objects generated for species %s; skipping.", species_id
            )
            continue

        # Amino acid usage (unweighted and abundance-weighted)
        freqs_unweighted = compute_proteome_aa_usage(
            proteins, weighted=False, allowed_aas=STANDARD_AA
        )
        freqs_weighted = compute_proteome_aa_usage(
            proteins, weighted=True, allowed_aas=STANDARD_AA
        )

        row_unweighted = {"species_id": species_id}
        row_unweighted.update(freqs_unweighted)
        usage_unweighted_rows.append(row_unweighted)

        row_weighted = {"species_id": species_id}
        row_weighted.update(freqs_weighted)
        usage_weighted_rows.append(row_weighted)

        # Relationship between usage and cost
        rel_result = correlate_usage_with_cost(
            species_id=species_id,
            freqs=freqs_weighted,
            costs=aa_costs,
        )
        relationship_results.append(rel_result)

    # Save tables
    if usage_unweighted_rows:
        df_unw = pd.DataFrame(usage_unweighted_rows).set_index("species_id")
        df_unw.to_csv(outdir / "amino_acid_usage_unweighted.tsv", sep="\t")
        logger.info(
            "Saved unweighted AA usage table for %d species.",
            df_unw.shape[0],
        )

    if usage_weighted_rows:
        df_w = pd.DataFrame(usage_weighted_rows).set_index("species_id")
        df_w.to_csv(outdir / "amino_acid_usage_weighted.tsv", sep="\t")
        logger.info(
            "Saved abundance-weighted AA usage table for %d species.",
            df_w.shape[0],
        )

    if relationship_results:
        df_rel = results_to_dataframe(relationship_results).set_index("species_id")
        df_rel.to_csv(outdir / "amino_acid_cost_correlations.tsv", sep="\t")
        logger.info(
            "Saved amino-acid usage vs cost correlation table for %d species.",
            df_rel.shape[0],
        )

    logger.info("Pipeline finished.")


if __name__ == "__main__":
    main()
