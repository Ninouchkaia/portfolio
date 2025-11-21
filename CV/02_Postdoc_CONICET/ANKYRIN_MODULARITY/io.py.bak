# io.py
import logging
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

from Bio import SeqIO  # type: ignore

logger = logging.getLogger(__name__)


def ensure_dir(path: Path) -> None:
    """Créer un répertoire (et parents) s'il n'existe pas."""
    path.mkdir(parents=True, exist_ok=True)


def read_uniprot_ids_from_fasta(fasta_path: Path) -> List[str]:
    """
    Reproduit le pattern record.id[3:9] pour extraire l'UNP ID
    depuis un multifasta UniProt-style.
    """
    ids: List[str] = []
    with fasta_path.open("r") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            uniprot_id = record.id[3:9]
            ids.append(uniprot_id)
    logger.info("Read %d UniProt IDs from %s", len(ids), fasta_path)
    return ids


def read_pfam_table(path: Path) -> List[Tuple[str, int, int, str]]:
    """
    Lire un fichier PFAM tabulé de type:
    start\tend\t...\tdomain_id\t...

    Retourne une liste de tuples (start, end, domain_id, full_line_stripped).
    """
    entries: List[Tuple[str, int, int, str]] = []
    with path.open("r") as fh:
        for raw in fh:
            line = raw.rstrip("\n").split("\t")
            if len(line) < 6:
                continue
            start = int(line[1])
            end = int(line[2])
            domain_id = line[5]
            entries.append((domain_id, start, end, raw.rstrip("\n")))
    logger.debug("Read %d PFAM entries from %s", len(entries), path)
    return entries


def read_simple_two_column_int_table(path: Path) -> Dict[str, int]:
    """
    Lire un fichier tabulé simple: key\tvalue
    (pour counts, fréquences, etc.)
    """
    d: Dict[str, int] = {}
    with path.open("r") as fh:
        header = next(fh, None)
        for raw in fh:
            cols = raw.rstrip("\n").split("\t")
            if len(cols) < 2:
                continue
            key = cols[0]
            try:
                value = int(cols[1])
            except ValueError:
                logger.warning("Skipping line (non-int value): %s", raw.strip())
                continue
            d[key] = value
    logger.info("Read %d entries from %s", len(d), path)
    return d


# --- Extraction FASTA depuis positions domaine/homologue --- #

from Bio.SeqRecord import SeqRecord  # type: ignore
from Bio.Seq import Seq  # type: ignore


# Dans io.py (mise à jour de extract_domain_sequences)

def extract_domain_sequences(
    fasta_path: Path,
    pfam_domain_file: Path,
    output_fasta: Path | None = None,
) -> Path:
    """
    Extraction propre des domaines PFAM.
    Si output_fasta n'est pas donné, on génère un nom propre automatiquement.
    """
    from .config import DOMAIN_FASTA_DIR

    ensure_dir(DOMAIN_FASTA_DIR)

    # Déduire UNIPROT du nom du fichier PFAM
    parts = pfam_domain_file.stem.split("_")
    if len(parts) < 5:
        logger.warning("Nom inattendu pour fichier PFAM: %s", pfam_domain_file)
        uniprot = parts[-1]
    else:
        uniprot = parts[-4]  # conforme à tes scripts

    if output_fasta is None:
        output_fasta = DOMAIN_FASTA_DIR / f"{uniprot}_domains.fasta"

    # Index du FASTA
    fasta_index = SeqIO.to_dict(SeqIO.parse(str(fasta_path), "fasta"))

    entries = read_pfam_table(pfam_domain_file)
    seq_records: List[SeqRecord] = []

    for domain_id, start, end, _ in entries:
        if uniprot not in fasta_index:
            logger.warning("UNIPROT %s absent du FASTA %s", uniprot, fasta_path)
            continue

        subseq = fasta_index[uniprot].seq[start - 1 : end]
        rec_id = f"{uniprot}_{domain_id}_{start}_{end}"
        seq_records.append(SeqRecord(Seq(str(subseq)), id=rec_id, description=""))

    with output_fasta.open("w") as out:
        SeqIO.write(seq_records, out, "fasta")

    logger.info(
        "Extracted %d domain sequences for %s → %s",
        len(seq_records), uniprot, output_fasta,
    )
    return output_fasta



def extract_independent_sequences(
    fasta_path: Path,
    cutoff_start: int,
    cutoff_end: int,
    output_fasta: Path | None = None,
) -> Path:
    """
    Génère un FASTA sans la région [cutoff_start, cutoff_end].
    Nomme automatiquement le fichier dans INDEP_FASTA_DIR.
    """
    from .config import INDEP_FASTA_DIR

    ensure_dir(INDEP_FASTA_DIR)

    fasta_index = SeqIO.to_dict(SeqIO.parse(str(fasta_path), "fasta"))

    seq_records: List[SeqRecord] = []
    for uniprot, rec in fasta_index.items():
        full = rec.seq
        new_seq = full[: cutoff_start - 1] + full[cutoff_end:]

        rec_id = f"{uniprot}_indep_without_{cutoff_start}_{cutoff_end}"
        seq_records.append(SeqRecord(new_seq, id=rec_id, description=""))

    if output_fasta is None:
        output_fasta = INDEP_FASTA_DIR / f"indep_without_{cutoff_start}_{cutoff_end}.fasta"

    with output_fasta.open("w") as fh:
        SeqIO.write(seq_records, fh, "fasta")

    logger.info(
        "Generated independent sequences without [%d:%d] → %s",
        cutoff_start, cutoff_end, output_fasta,
    )
    return output_fasta




def find_pfam_files(uniprot_id: str, num_homologs: int) -> Tuple[Path, Path]:
    """
    Trouve automatiquement les fichiers PFAM query/homolog à partir du pattern :

        pfam_domains_in_{UNIPROT}_MaxHomologs_{N}_{query|homolog}.txt

    dans PFAM_QUERY_DIR et PFAM_HOMOLOG_DIR.

    Retourne (query_file, homolog_file).
    """
    from .config import PFAM_QUERY_DIR, PFAM_HOMOLOG_DIR, DEFAULT_FILENAME_PATTERN

    query_pattern = DEFAULT_FILENAME_PATTERN.format(
        unp=uniprot_id, N=num_homologs, kind="query"
    )
    homolog_pattern = DEFAULT_FILENAME_PATTERN.format(
        unp=uniprot_id, N=num_homologs, kind="homolog"
    )

    query_file = PFAM_QUERY_DIR / query_pattern
    homolog_file = PFAM_HOMOLOG_DIR / homolog_pattern

    if not query_file.exists():
        logger.warning("PFAM query file not found: %s", query_file)
    if not homolog_file.exists():
        logger.warning("PFAM homolog file not found: %s", homolog_file)

    return query_file, homolog_file


def concatenate_pfam_outputs_by_threshold(
    input_dir: Path,
    output_file: Path,
    min_occurrences: int = 1,
) -> Path:
    """
    Concatène tous les fichiers PFAM dans un dossier (ex: pfam_queries/)
    et applique un seuil minimal d'occurrences par domaine.

    Reproduit la logique de concaten_pfam_outputs_BD_2015.py.

    Format attendu des fichiers PFAM :
        uniprot_id \t start \t end \t domain_id

    Sortie : fichier tabulé domaine \t total_occurrences
    """
    ensure_dir(output_file.parent)

    domain_counts: Dict[str, int] = {}

    for pfam_file in sorted(input_dir.glob("*.txt")):
        with pfam_file.open("r") as fh:
            for raw in fh:
                cols = raw.rstrip("\n").split("\t")
                if len(cols) < 4:
                    continue
                domain_id = cols[3]
                domain_counts[domain_id] = domain_counts.get(domain_id, 0) + 1

    with output_file.open("w") as out:
        out.write("domain_name\toccurrences\n")
        for domain, count in sorted(domain_counts.items()):
            if count >= min_occurrences:
                out.write(f"{domain}\t{count}\n")

    logger.info(
        "Concatenated %d PFAM domain entries → %s (threshold=%d)",
        sum(domain_counts.values()), output_file, min_occurrences
    )
    return output_file
