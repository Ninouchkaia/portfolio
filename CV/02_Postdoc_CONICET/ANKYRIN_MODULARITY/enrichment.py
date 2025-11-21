# enrichment.py
import logging
import math
from pathlib import Path
from typing import Dict

from .config import INTERMEDIATE_DIR, RESULTS_DIR
from .io import ensure_dir, read_simple_two_column_int_table

logger = logging.getLogger(__name__)


def compute_domain_enrichment_from_counts(
    domain_counts_file: Path,
    pfam_counts_file: Path,
    num_proteins_subfamily: int,
    num_proteins_swissprot: int,
    output_file: Path,
) -> Path:
    """
    Version généralisée du bloc qui écrit Pfam_domains_in_BD_2038.txt :contentReference[oaicite:8]{index=8}

    domain_counts_file : table "domain_name <tab> # in subfamily"
    pfam_counts_file   : table "domain_name <tab> # in Uniprot"
    """
    if output_file.exists():
        logger.info("Enrichment file already exists: %s", output_file)
        return output_file

    ensure_dir(output_file.parent)

    sub_counts = read_simple_two_column_int_table(domain_counts_file)
    pfam_counts = read_simple_two_column_int_table(pfam_counts_file)

    with output_file.open("w") as out:
        out.write("domain_name\tcount_subfamily\tcount_uniprot\texpected_count\tlog_obs_exp\n")
        for domain_name, sub_count in sorted(sub_counts.items()):
            if domain_name not in pfam_counts:
                logger.debug("Domain %s not found in pfam_counts, skipping", domain_name)
                continue
            unp_count = pfam_counts[domain_name]
            exp_count = unp_count * num_proteins_subfamily / num_proteins_swissprot
            if exp_count == 0:
                continue
            obs_exp = sub_count / exp_count
            log_obs_exp = math.log(obs_exp)
            out.write(
                f"{domain_name}\t{sub_count}\t{unp_count}\t{exp_count:.6f}\t{log_obs_exp:.6f}\n"
            )

    logger.info("Wrote domain enrichment table to %s", output_file)
    return output_file


def _load_conservation_percentages(path: Path) -> Dict[str, float]:
    """
    Lit domain_conservation_percentages_*.txt.
    """
    d: Dict[str, float] = {}
    with path.open("r") as fh:
        header = next(fh, None)
        for raw in fh:
            cols = raw.rstrip("\n").split()
            if len(cols) < 2:
                continue
            try:
                d[cols[0]] = float(cols[1].replace(",", "."))
            except ValueError:
                continue
    return d


def color_enrichment_by_conservation(
    enrichment_zscore_file: Path,
    conservation_percent_file: Path,
    num_homologs: int,
    output_prefix: str,
) -> Path:
    """
    Version factorisée de domain_enrich_conserv_2015.py et _BD_2015.py 

    enrichment_zscore_file :
        - Pfam_domains_in_Ank1234_Zscores.txt
        - ou Pfam_domains_in_BD_2038_Zscores.txt
      colonnes (au moins) : domain_name, ..., log(obs/exp), Zscore

    conservation_percent_file :
        - domain_conservation_percentages_*.txt

    output_prefix :
        - ex: "Ank1234" ou "BD_2038"
    """
    ensure_dir(RESULTS_DIR)
    output_file = RESULTS_DIR / f"Pfam_domains_in_{output_prefix}_MaxHomologs_{num_homologs}_Zscores_color.txt"

    if output_file.exists():
        logger.info("Colored enrichment file already exists: %s", output_file)
        return output_file

    conserv = _load_conservation_percentages(conservation_percent_file)
    if not conserv:
        logger.warning("Empty conservation file: %s", conservation_percent_file)

    avg_conserv = sum(conserv.values()) / len(conserv) if conserv else 0.0
    logger.info(
        "Average domain conservation (N=%d) = %.4f from %s",
        len(conserv),
        avg_conserv,
        conservation_percent_file,
    )

    with enrichment_zscore_file.open("r") as fh, output_file.open("w") as out:
        header = next(fh, None)
        out.write(
            f"domain_name\tdomain_enrichment_log(obs/exp)\tZscore_enrichment\t"
            f"domain_conservation_over_{num_homologs}Homologs\tcolor\n"
        )

        for raw in fh:
            cols = raw.rstrip("\n").split("\t")
            if len(cols) < 6:
                continue
            domain_name = cols[0]
            try:
                domain_enrichment = float(cols[4].replace(",", "."))
                domain_zscore = float(cols[5].replace(",", "."))
            except ValueError:
                continue

            if domain_name not in conserv:
                logger.debug("Domain %s not in conservation dict", domain_name)
                continue

            domain_conservation = conserv[domain_name]
            color = '"red"' if (domain_zscore > 1 and domain_conservation > avg_conserv) else '"black"'
            out.write(
                f"{domain_name}\t{domain_enrichment:.6f}\t{domain_zscore:.6f}\t"
                f"{domain_conservation:.6f}\t{color}\n"
            )

    logger.info("Wrote colored enrichment table to %s", output_file)
    return output_file
