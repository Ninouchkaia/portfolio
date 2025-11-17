# abm_pipeline/config.py

from pathlib import Path

# Racine du projet (à adapter si nécessaire)
PROJECT_ROOT = Path(__file__).resolve().parent.parent

# Liste patients + % monocytes (reprend ce que tu as dans tes scripts)
PATIENTS_WITH_MONO = [
    "CRE1704-1.1%",
    "LAU1405-2.5%",
    "DES2105-1.25%",
    "ORE1706-0.68%",
    "CAS1802-1.04%",
    "GER160522-0.45%",
    "REI230522-0.95%",
    "PUJ240522-0.34%",
    "LAR300522-0.21%",
    "CAZ310522-3.48%",
]

PATIENT_IDS = [p.split("-")[0] for p in PATIENTS_WITH_MONO]
