#!/usr/bin/env python
from __future__ import annotations

import argparse
import logging

from src.logging_utils import setup_logging
from src.interactors import compute_all_orders_for_all_viruses
from src.tables import build_gene_lists_from_interactors, build_gene_virus_table
from src.enrichment_wrapper import run_r_enrichment


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Systemic effects of SARS-CoV-2 on host cellular functions – pipeline"
    )

    subparsers = parser.add_subparsers(dest="command", required=True)

    # 1. network → interactors
    p_net = subparsers.add_parser("compute_interactors", help="Compute order-k interactors")
    p_net.add_argument("--max-order", type=int, default=2)

    # 2. gene lists for R
    p_genelist = subparsers.add_parser("gene_lists", help="Build gene lists for enrichment")
    p_genelist.add_argument("--range-mode", type=str, default="direct_and_second")

    # 3. gene-virus table
    p_table = subparsers.add_parser("gene_virus_table", help="Build gene-virus table")
    p_table.add_argument("--range-mode", type=str, default="direct_and_second")

    # 4. enrichment via R
    p_enrich = subparsers.add_parser("enrich_reactome", help="Run Reactome enrichment in R")
    p_enrich.add_argument("--range-mode", type=str, default="direct_and_second")
    p_enrich.add_argument("--script-name", type=str, default="enrich_reactome_compareCluster.R")

    return parser.parse_args()


def main() -> None:
    setup_logging()
    logger = logging.getLogger("covid_pipeline")

    args = parse_args()

    if args.command == "compute_interactors":
        logger.info("Step: compute_interactors (max_order=%d)", args.max_order)
        compute_all_orders_for_all_viruses(max_order=args.max_order)

    elif args.command == "gene_lists":
        logger.info("Step: gene_lists (range_mode=%s)", args.range_mode)
        build_gene_lists_from_interactors(range_mode=args.range_mode)

    elif args.command == "gene_virus_table":
        logger.info("Step: gene_virus_table (range_mode=%s)", args.range_mode)
        build_gene_virus_table(range_mode=args.range_mode)

    elif args.command == "enrich_reactome":
        logger.info(
            "Step: enrich_reactome (range_mode=%s, script=%s)",
            args.range_mode, args.script_name
        )
        gene_list_file = build_gene_lists_from_interactors(range_mode=args.range_mode)
        run_r_enrichment(gene_list_file, script_name=args.script_name)

    else:
        raise RuntimeError(f"Unknown command: {args.command}")


if __name__ == "__main__":
    main()
