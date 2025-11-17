# clans.py
import logging
from pathlib import Path
from typing import Dict

from .config import PFAM_CLANS_FILE, PFAM_FREQUENCIES_FILE, PFAM_CLAN_FREQ_FILE
from .io import ensure_dir

logger = logging.getLogger(__name__)


def load_family_to_clan_mapping(clans_file: Path = PFAM_CLANS_FILE) -> Dict[str, str]:
    """
    Construit un dict PFAM-family → clan à partir de Pfam-A.clans.txt.
    Reste fidèle à ton parsing (colonnes 4 et 2). :contentReference[oaicite:1]{index=1}
    """
    mapping: Dict[str, str] = {}
    with clans_file.open("r") as fh:
        header = next(fh, None)
        for raw in fh:
            cols = raw.rstrip("\n").split("\t")
            if len(cols) < 4:
                continue
            clan = cols[1]
            family = cols[3]
            if clan != r"\N":
                mapping[family] = clan
    logger.info("Loaded %d PFAM family→clan mappings from %s", len(mapping), clans_file)
    return mapping


def aggregate_pfam_frequencies_by_clan(
    pfam_freq_file: Path = PFAM_FREQUENCIES_FILE,
    clans_file: Path = PFAM_CLANS_FILE,
    output_file: Path = PFAM_CLAN_FREQ_FILE,
) -> Path:
    """
    Reprend le comportement de adapt_frequencies_pfam_to_clans.py :
    - lit FrequenciesPfam_mapped.txt
    - regroupe par clan
    - écrit FrequenciesPfam_mapped_CLANS.txt :contentReference[oaicite:2]{index=2}
    """
    family_to_clan = load_family_to_clan_mapping(clans_file)
    clan_counts: Dict[str, int] = {}

    with pfam_freq_file.open("r") as fh:
        header = next(fh, None)
        for raw in fh:
            cols = raw.rstrip("\n").split("\t")
            if len(cols) < 2:
                continue
            domain = cols[1]
            try:
                count = int(cols[-1])
            except ValueError:
                logger.warning("Skipping line with non-int count: %s", raw.strip())
                continue

            clan = family_to_clan.get(domain, domain)
            clan_counts[clan] = clan_counts.get(clan, 0) + count

    ensure_dir(output_file.parent)
    with output_file.open("w") as out:
        out.write("Pfam_clan\tcount_in_swissprot\n")
        for clan, total in sorted(clan_counts.items()):
            out.write(f"{clan}\t{total}\n")

    logger.info("Wrote %d clan counts to %s", len(clan_counts), output_file)
    return output_file


def annotate_domains_with_clans(
    input_file: Path,
    clans_file: Path = PFAM_CLANS_FILE,
    output_file: Path | None = None,
) -> Path:
    """
    Equivalent à rename_domains_to_clans.py :
    ajoute une colonne CLAN à une table contenant un domaine PFAM en première colonne. :contentReference[oaicite:3]{index=3}
    """
    family_to_clan = load_family_to_clan_mapping(clans_file)
    if output_file is None:
        output_file = input_file.with_name(f"{input_file.stem}_with_clans.txt")

    ensure_dir(output_file.parent)

    with input_file.open("r") as fh, output_file.open("w") as out:
        lines = fh.readlines()
        if not lines:
            return output_file

        header_cols = lines[0].rstrip("\n").split("\t")
        out.write(f"{header_cols[0]}\tCLAN\t" + "\t".join(header_cols[1:]) + "\n")

        for raw in lines[1:]:
            cols = raw.rstrip("\n").split("\t")
            if not cols:
                continue
            domain_name = cols[0]
            clan_name = family_to_clan.get(domain_name, "NULL")
            out.write(f"{domain_name}\t{clan_name}\t" + "\t".join(cols[1:]) + "\n")

    logger.info("Annotated %s → %s with clan column", input_file, output_file)
    return output_file
