# config.py
from pathlib import Path

# Racine du projet
PROJECT_ROOT = Path(__file__).resolve().parent

# ------------------------
# STRUCTURE DES DONNÉES
# ------------------------
DATA_DIR = PROJECT_ROOT / "data"
RAW_DIR = DATA_DIR / "raw"
INTERMEDIATE_DIR = DATA_DIR / "intermediate"
RESULTS_DIR = DATA_DIR / "results"

# Sous-dossiers PFAM (conformes à ton projet original)
PFAM_QUERY_DIR = INTERMEDIATE_DIR / "pfam_queries"
PFAM_HOMOLOG_DIR = INTERMEDIATE_DIR / "pfam_homologs"

# Sous-dossiers de résultats FASTA
SEQUENCES_DIR = RESULTS_DIR / "sequences"
DOMAIN_FASTA_DIR = SEQUENCES_DIR / "domains"
INDEP_FASTA_DIR = SEQUENCES_DIR / "independent"

# ------------------------
# FICHIERS D'ENTRÉE
# ------------------------
PFAM_CLANS_FILE = RAW_DIR / "Pfam-A.clans.txt"
PFAM_FREQUENCIES_FILE = RAW_DIR / "FrequenciesPfam_mapped.txt"

# Multifasta principaux
ANK_FASTA = RAW_DIR / "Ank_uniref50_string_mapped_522.fasta"
BD_FASTA = RAW_DIR / "binding_partners_2038.fasta"

# ------------------------
# FICHIERS DE SORTIE
# ------------------------
PFAM_CLAN_FREQ_FILE = INTERMEDIATE_DIR / "FrequenciesPfam_mapped_CLANS.txt"

# Paramètres par défaut
DEFAULT_MAX_HOMOLOGS = 400
DEFAULT_FILENAME_PATTERN = (
    "pfam_domains_in_{unp}_MaxHomologs_{N}_{kind}.txt"
)
