# elm.py
import logging
from pathlib import Path
from collections import Counter, defaultdict
from typing import Dict, List, Tuple

from .io import ensure_dir

logger = logging.getLogger(__name__)


# ------------------------------------------------
# 1. Lecture des tables ELM (SLiMs) prédites
# ------------------------------------------------
def read_elm_table(path: Path) -> List[Tuple[str, str]]:
    """
    Lit un fichier ELM (SLiM predictions), supposé contenir :
        uniprot_id \t elm_name
    """
    els: List[Tuple[str, str]] = []
    with path.open("r") as fh:
        header = next(fh, None)
        for raw in fh:
            cols = raw.rstrip("\n").split("\t")
            if len(cols) < 2:
                continue
            els.append((cols[0], cols[1]))
    logger.info("Read %d ELM predictions from %s", len(els), path)
    return els


# ------------------------------------------------
# 2. Comptage par ELM
# ------------------------------------------------
def count_elms(elm_tuples: List[Tuple[str, str]]) -> Dict[str, int]:
    """
    Compte le nombre d'occurrences de chaque ELM.
    """
    counter = Counter([elm for _, elm in elm_tuples])
    logger.info("Counted %d distinct ELMs", len(counter))
    return dict(counter)


def write_elm_counts(counts: Dict[str, int], output_file: Path) -> Path:
    ensure_dir(output_file.parent)
    with output_file.open("w") as out:
        out.write("elm_name\tcount\n")
        for elm, cnt in sorted(counts.items()):
            out.write(f"{elm}\t{cnt}\n")
    logger.info("Wrote ELM counts → %s", output_file)
    return output_file


# ------------------------------------------------
# 3. Enrichissement ELM (log(obs/exp))
# ------------------------------------------------
def compute_elm_enrichment(
    elm_counts_subfamily: Dict[str, int],
    elm_counts_background: Dict[str, int],
    num_subfamily_proteins: int,
    num_background_proteins: int,
    output_file: Path,
) -> Path:
    """
    Equivalent de domain enrichment, mais pour les ELMs.
    """
    ensure_dir(output_file.parent)

    with output_file.open("w") as out:
        out.write("elm_name\tcount_subfamily\tcount_background\texpected\tlog_obs_exp\n")

        for elm, sub_count in sorted(elm_counts_subfamily.items()):
            bg_count = elm_counts_background.get(elm, 0)
            expected = bg_count * num_subfamily_proteins / num_background_proteins

            if expected > 0 and sub_count > 0:
                log_obs_exp = math.log(sub_count / expected)
            else:
                log_obs_exp = 0.0

            out.write(
                f"{elm}\t{sub_count}\t{bg_count}\t{expected:.4f}\t{log_obs_exp:.6f}\n"
            )

    logger.info("Computed ELM enrichment → %s", output_file)
    return output_file


# ------------------------------------------------
# 4. Co-occurrence ELM ↔ Domain dans interacting pairs
# ------------------------------------------------
def compute_elm_domain_cooccurrences(
    elm_table: List[Tuple[str, str]],
    domain_table: List[Tuple[str, str]],
    output_file: Path,
) -> Path:
    """
    Approche inspirée de ton draft (figure 2) :
    - ELM dans un partenaire
    - domaine dans son interacteur
    - comptage des couples (ELM, domain)

    elm_table: [(unp, elm)]
    domain_table: [(unp, domain)]
    """
    ensure_dir(output_file.parent)

    elms_by_unp = defaultdict(list)
    for unp, elm in elm_table:
        elms_by_unp[unp].append(elm)

    domains_by_unp = defaultdict(list)
    for unp, dom in domain_table:
        domains_by_unp[unp].append(dom)

    pairs_count = Counter()

    all_unp = set(elms_by_unp.keys()) | set(domains_by_unp.keys())

    for unp in all_unp:
        for elm in elms_by_unp.get(unp, []):
            for dom in domains_by_unp.get(unp, []):
                pairs_count[(elm, dom)] += 1

    with output_file.open("w") as out:
        out.write("elm\tdomain\tcount\n")
        for (elm, dom), cnt in pairs_count.items():
            out.write(f"{elm}\t{dom}\t{cnt}\n")

    logger.info("Computed %d ELM-domain cooccurrences → %s", len(pairs_count), output_file)
    return output_file
