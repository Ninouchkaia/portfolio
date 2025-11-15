# abm_pipeline/parameter_exploration/utils.py

import logging
from typing import List, Dict

from abm_pipeline.config import PATIENTS_WITH_MONO, PATIENT_IDS

# --- logging / debug -------------------------------------------------------

logging.basicConfig(
    level=logging.INFO,
    format="%(levelname)s | %(name)s | %(message)s",
)

logger = logging.getLogger("abm_pipeline")


def debug(msg: str):
    """Wrapper simple pour les messages de debug."""
    logger.debug(msg)


# --- patients --------------------------------------------------------------

def get_patient_ids() -> List[str]:
    """['CRE1704', 'LAU1405', ...]."""
    return PATIENT_IDS.copy()


def get_patients_with_mono() -> List[str]:
    """['CRE1704-1.1%', ...]."""
    return PATIENTS_WITH_MONO.copy()


def build_patient_mono_dict() -> Dict[str, float]:
    """
    {'CRE1704': 1.1, ...}
    Utile pour les scripts qui ont besoin du % mono.
    """
    patient_dict = {}
    for patient in PATIENTS_WITH_MONO:
        name = patient.split("-")[0]
        mono = float(patient.split("-")[1].rstrip("%"))
        patient_dict[name] = mono
    return patient_dict
