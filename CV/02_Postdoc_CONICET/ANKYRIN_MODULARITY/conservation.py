# conservation.py
import logging
from pathlib import Path
from collections import Counter, defaultdict
from typing import Dict, List, Tuple

from .config import INTERMEDIATE_DIR
from .io import ensure_dir, read_uniprot_ids_from_fasta

logger = logging.getLogger(__name__)


def _load_domain_occurrences(path: Path) -> List[Tuple[str, int, int]]:
    """
    Lit un fichier PFAM simplifié:
    uniprot_id \t start \t end \t domain_id
    (ce format correspond aux pfam_domains_in_*_query.txt/homolog.txt utilisés en legacy).
    """
    occurrences: List[Tuple[str, int, int]] = []
    with path.open("r") as fh:
        for raw in fh:
            cols = raw.rstrip("\n").split("\t")
            if len(cols) < 4:
                continue
            unp = cols[0]
            start = int(cols[1])
            end = int(cols[2])
            domain = cols[3]
            occurrences.append((domain, start, end))
    return occurrences


def compute_domain_conservation_for_fasta(
    fasta_path: Path,
    num_homologs: int,
    base_dir: Path,
) -> Tuple[Path, Path]:
    """
    Equivalent structuré à check_domain_conservation_BD_2015.py pour un multifasta donné.
    - lit les uniprot IDs dans fasta_path
    - pour chaque uniprot, lit les fichiers query/homolog
    - calcule la conservation de chaque domaine
    - écrit:
      * conservation_state_...txt
      * domain_conservation_percentages_...txt

    On supppose des chemins du type:
    {basename}_MaxHomologs_MAX{N}_queries_pfam_output/pfam_domains_in_{UNP}_MaxHomologs_{N}_query.txt
    et mêmes patterns pour homologs. :contentReference[oaicite:5]{index=5}
    """
    basename = fasta_path.stem  # ex: binding_partners_2038
    out_state = INTERMEDIATE_DIR / f"conservation_state_{basename}_{num_homologs}.txt"
    out_percent = INTERMEDIATE_DIR / f"domain_conservation_percentages_{basename}_{num_homologs}.txt"

    if out_percent.exists():
        logger.info("Conservation already computed for %s (N=%d)", basename, num_homologs)
        return out_state, out_percent

    ensure_dir(INTERMEDIATE_DIR)

    uniprot_ids = read_uniprot_ids_from_fasta(fasta_path)
    logger.info("Computing conservation for %d proteins in %s (N=%d)", len(uniprot_ids), fasta_path, num_homologs)

    # (domaine → [nb_conservé, nb_non_conservé])
    conservation_counts: Dict[str, List[int]] = defaultdict(lambda: [0, 0])

    with out_state.open("w") as state_out:
        state_out.write("uniprot_id\tdomain_id\tis_conserved\n")

        for unp in uniprot_ids:
            query_file = base_dir / f"{basename}_MaxHomologs_MAX{num_homologs}_queries_pfam_output" / \
                f"pfam_domains_in_{unp}_MaxHomologs_{num_homologs}_query.txt"
            homol_file = base_dir / f"{basename}_MaxHomologs_MAX{num_homologs}_homologs_pfam_output" / \
                f"pfam_domains_in_{unp}_MaxHomologs_{num_homologs}_homolog.txt"

            if not query_file.exists() or not homol_file.exists():
                logger.warning("Missing PFAM files for %s: %s or %s", unp, query_file, homol_file)
                continue

            # Collecte des domaines chez le query
            query_domains = _load_domain_occurrences(query_file)
            # Collecte des domaines homologues
            homolog_domains = _load_domain_occurrences(homol_file)
            homolog_counter = Counter([d for d, _, _ in homolog_domains])

            for domain_id, _, _ in query_domains:
                conserved = domain_id in homolog_counter
                state_out.write(f"{unp}\t{domain_id}\t{int(conserved)}\n")
                if conserved:
                    conservation_counts[domain_id][0] += 1
                else:
                    conservation_counts[domain_id][1] += 1

    # Calcul des pourcentages
    with out_percent.open("w") as perc_out:
        perc_out.write("domain_name\tconservation_percentage\n")
        for domain, (n_conserved, n_not) in sorted(conservation_counts.items()):
            total = n_conserved + n_not
            if total == 0:
                continue
            perc = n_conserved / total
            perc_out.write(f"{domain}\t{perc:.6f}\n")

    logger.info("Wrote conservation state to %s and percentages to %s", out_state, out_percent)
    return out_state, out_percent
