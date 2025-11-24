# abm_pipeline/config.py

from pathlib import Path

# ---------------------------------------------------------------------------
# Root directory
# ---------------------------------------------------------------------------

PROJECT_ROOT = Path(__file__).resolve().parents[1]

# ---------------------------------------------------------------------------
# Data folders
# ---------------------------------------------------------------------------

DATA_DIR = PROJECT_ROOT / "data"
EXPERIMENTAL_DATA_DIR = DATA_DIR / "experimental"
PARETO_DIR = DATA_DIR / "pareto"

# ---------------------------------------------------------------------------
# Results folders
# ---------------------------------------------------------------------------

RESULTS_DIR = PROJECT_ROOT / "results"

BEHAVIORSAPCE_GENERAL_DIR = RESULTS_DIR / "behaviorspace" / "general_model"
BEHAVIORSAPCE_PATIENT_DIR = RESULTS_DIR / "behaviorspace" / "patient_specific_models"

VALIDATION_RESULTS_DIR = RESULTS_DIR / "validation"
ADVANCED_RESULTS_DIR = RESULTS_DIR / "advanced_analysis"
SENSITIVITY_RESULTS_DIR = RESULTS_DIR / "sensitivity"

# ---------------------------------------------------------------------------
# Experimental data
# ---------------------------------------------------------------------------

EXPERIMENTAL_FILES = {
    "general": EXPERIMENTAL_DATA_DIR / "expData.csv"
}

EXPERIMENTAL_COLUMNS = {
    "time": "day",
    "viability": "viability",
    "concentration": "concentration",
}

# ---------------------------------------------------------------------------
# Patient information
# ---------------------------------------------------------------------------

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

# ---------------------------------------------------------------------------
# Utilities
# ---------------------------------------------------------------------------

def ensure_dirs_exist():
    for d in [
        DATA_DIR,
        EXPERIMENTAL_DATA_DIR,
        PARETO_DIR,
        RESULTS_DIR,
        BEHAVIORSAPCE_GENERAL_DIR,
        BEHAVIORSAPCE_PATIENT_DIR,
        VALIDATION_RESULTS_DIR,
        ADVANCED_RESULTS_DIR,
        SENSITIVITY_RESULTS_DIR,
    ]:
        d.mkdir(parents=True, exist_ok=True)
````

---

# 🟦 **3. Fichier `sensitivity_config.py` avec listes organisées**

```python
# abm_pipeline/sensitivity_config.py

from typing import List

# ---------------------------------------------------------------------------
# Prediction experiments
# ---------------------------------------------------------------------------

EXP_PREDICTION: List[str] = [
    "prediction-day10-varying-mono-init",
    "prediction-day8-varying-mono-init",
    "prediction-day9-varying-mono-init",
    "prediction-day12-varying-mono-init",
    "prediction-patient1-1.1",
    "prediction-patient2-2.5",
    "prediction-patient3-1.25",
]

# ---------------------------------------------------------------------------
# Stochasticity experiments
# ---------------------------------------------------------------------------

EXP_STOCHASTICITY: List[str] = [
    "Explo-stochasticity-100-simulations",
]

# ---------------------------------------------------------------------------
# Sensitivity experiments
# ---------------------------------------------------------------------------

EXP_SENSITIVITY: List[str] = [
    "PerturbM2PhagoEff2",
    "PerturbMonoPhagoEff",
    "PerturbNLCPhagoEff",
    "PerturbM2KillEff",
    "PerturbNLCProtection",
    "PerturbApoMov",
    "PerturbNeedSigMov",
    "PerturbDiffStd",
    "PerturbDiffMean",
    "PerturbLayers",
    "PerturbCLLSensDist",
    "PerturbMonoSensDist",
    "PerturbNLCSensDist",
    "PerturbMacroSensDist",
    "PerturbNLCThreshold",
    "PerturbSigInit",
    "PerturbSigInitStd",
]

# ---------------------------------------------------------------------------
# Global list (used by CLI)
# ---------------------------------------------------------------------------

EXP_LIST: List[str] = (
    EXP_PREDICTION
    + EXP_STOCHASTICITY
    + EXP_SENSITIVITY
)