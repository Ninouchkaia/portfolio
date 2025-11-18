from pathlib import Path

# Racine du projet (fichier analysis.py à la racine)
PROJECT_ROOT = Path(__file__).resolve().parents[1]

DATA_RAW = PROJECT_ROOT / "data" / "raw" / "Virus_host_interactomes_thresh25" / "thresh_0.25"
DATA_INTERMEDIATE = PROJECT_ROOT / "data" / "intermediate"
DATA_RESULTS = PROJECT_ROOT / "data" / "results"
R_SCRIPTS_DIR = PROJECT_ROOT / "r"

INTERACTORS_DIR = DATA_INTERMEDIATE / "interactors"
GENE_LISTS_DIR = DATA_INTERMEDIATE / "gene_lists"
ENRICHMENT_DIR = DATA_INTERMEDIATE / "enrichment"

TABLES_DIR = DATA_RESULTS / "tables"
FIGURES_DIR = DATA_RESULTS / "figures"

# Création des dossiers au besoin
for p in [
    DATA_RAW,
    DATA_INTERMEDIATE,
    DATA_RESULTS,
    INTERACTORS_DIR,
    GENE_LISTS_DIR,
    ENRICHMENT_DIR,
    TABLES_DIR,
    FIGURES_DIR,
]:
    p.mkdir(parents=True, exist_ok=True)
