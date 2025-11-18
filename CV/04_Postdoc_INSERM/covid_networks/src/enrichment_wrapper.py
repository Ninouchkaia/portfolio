from __future__ import annotations
from pathlib import Path
import subprocess
import logging

from .config import R_SCRIPTS_DIR, ENRICHMENT_DIR

logger = logging.getLogger(__name__)


def run_r_enrichment(
    gene_list_file: Path,
    script_name: str = "enrich_reactome_compareCluster.R",
    output_prefix: str = "enrichPathway",
) -> Path:
    """
    Appelle un script R qui lit un fichier genes_list_*.txt,
    fait les mappings + compareCluster, et écrit un TSV dans ENRICHMENT_DIR.
    """
    script_path = R_SCRIPTS_DIR / script_name
    out_tsv = ENRICHMENT_DIR / f"{output_prefix}.tsv"

    cmd = [
        "Rscript",
        str(script_path),
        str(gene_list_file),
        str(out_tsv),
    ]
    logger.info("Running Rscript: %s", " ".join(cmd))
    subprocess.run(cmd, check=True)

    logger.info("R enrichment results written to %s", out_tsv)
    return out_tsv
