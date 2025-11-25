from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path
import pandas as pd

THIS_DIR = Path(__file__).resolve().parent
SRC_DIR = THIS_DIR
sys.path.insert(0, str(SRC_DIR))

from paxdb.src.utils import setup_logging, resolve_path
from paxdb.src.fasta_parser import load_fasta_as_dict
from paxdb.src.abundance_loader import load_abundance_table, load_species_metadata
from paxdb.src.protein import Protein
from paxdb.src.aa_metrics import compute_proteome_aa_usage, load_amino_acid_costs, STANDARD_AA
from paxdb.src.relationships import correlate_usage_with_cost, results_to_dataframe

logger = logging.getLogger(__name__)


def load_mapping_file(species_id: str):
    """
    Auto-loads mapping file:
      paxdb/data/metadata/<species_id>_string_to_uniprot.tsv
    """
    path = resolve_path(f"paxdb/data/metadata/{species_id}_string_to_uniprot.tsv")
    try:
        df = pd.read_csv(path, sep="\t")
        return dict(zip(df["string_id"], df["fasta_id"]))
    except Exception as e:
        logger.warning(f"No mapping file found for {species_id}: {e}")
        return {}


def load_proteome_for_species(species_id: str,
                              fasta_path: Path,
                              abundance_path: Path,
                              id_col="string_id"):
    """
    Loads FASTA + Abundance + mapping
    and returns list of Protein objects with abundance.
    """

    seq_dict = load_fasta_as_dict(fasta_path)
    mapping = load_mapping_file(species_id)
    abund = load_abundance_table(abundance_path,
                                 id_col=id_col,
                                 abundance_col="abundance")

    # Apply mapping
    if mapping:
        abund[id_col] = abund[id_col].map(mapping).fillna("")

    proteins = []
    missing = 0

    for _, row in abund.iterrows():
        fid = row[id_col]
        if not fid or fid not in seq_dict:
            missing += 1
            continue
        proteins.append(Protein(fid, seq_dict[fid], float(row["abundance"])))

    if missing > 0:
        logger.warning(f"Missing {missing} proteins matching FASTA for {species_id}")

    return proteins


def parse_args():
    p = argparse.ArgumentParser(description="PaxDB AA analysis pipeline.")
    p.add_argument("--species-metadata", type=Path, required=True)
    p.add_argument("--aa-costs", type=Path, required=True)
    p.add_argument("--outdir", type=Path, default=Path("results"))
    p.add_argument("--logdir", type=Path, default=Path("paxdb/logs"))
    p.add_argument("--species", nargs="*", default=None)
    return p.parse_args()


def main():
    args = parse_args()
    setup_logging(args.logdir)

    logger.info("Starting PaxDB amino acid usage analysis.")
    logger.info(f"Species metadata: {args.species_metadata}")
    logger.info(f"AA costs: {args.aa_costs}")
    logger.info(f"Output dir: {args.outdir}")

    species_df = load_species_metadata(args.species_metadata)

    if args.species:
        species_df = species_df[species_df["species_id"].isin(args.species)]

    species_ids = list(species_df["species_id"])
    aa_costs = load_amino_acid_costs(args.aa_costs)

    results = []

    for _, row in species_df.iterrows():
        sid = row["species_id"]
        logger.info(f"Processing species {sid}")

        fasta_path = resolve_path(row["fasta_path"])
        abundance_path = resolve_path(row["abundance_path"])

        proteome = load_proteome_for_species(
            sid,
            fasta_path,
            abundance_path,
            id_col="string_id",
        )

        if not proteome:
            logger.warning(f"No proteins loaded for species {sid}; skipping.")
            continue

        usage = compute_proteome_aa_usage(proteome, STANDARD_AA)
        (args.outdir / sid).mkdir(parents=True, exist_ok=True)
        usage.to_csv(args.outdir / sid / "aa_usage.tsv", sep="\t")

        r = correlate_usage_with_cost(sid, usage, aa_costs)
        results.append(r)

    if results:
        df = results_to_dataframe(results)
        df.to_csv(args.outdir / "aa_cost_relationships.tsv", sep="\t", index=False)
    else:
        logger.warning("No relationship results generated.")

    logger.info("PaxDB amino acid usage pipeline finished.")


if __name__ == "__main__":
    main()
