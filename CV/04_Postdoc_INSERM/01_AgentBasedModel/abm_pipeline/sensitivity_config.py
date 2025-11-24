# abm_pipeline/sensitivity_config.py

"""
Organisation complète des expériences BehaviorSpace extraites
automatiquement du fichier ABM_NLC_CLL.nlogo.

Les expériences sont regroupées en trois catégories :
- PRÉDICTION (prediction-dayX / prediction-patientX)
- STOCHASTICITÉ (exploration stochastique)
- SENSIBILITÉ (PerturbX...)

Une liste globale EXP_LIST est fournie pour compatibilité avec
la CLI et les modules existants.
"""

from typing import List

# ---------------------------------------------------------------------------
# 1. PRÉDICTION — simulations destinées à évaluer le modèle sur des
#    conditions temporelles ou des paramètres patients.
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
# 2. STOCHASTICITÉ — expériences dédiées à l’évaluation du bruit interne du modèle.
# ---------------------------------------------------------------------------

EXP_STOCHASTICITY: List[str] = [
    "Explo-stochasticity-100-simulations",
]

# ---------------------------------------------------------------------------
# 3. SENSIBILITÉ — variations unifactorielle
#    utilisées par la pipeline Python et par OpenMOLE.
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
# 4. LISTE GLOBALE — utilisée par la CLI pour générer les scripts de sensibilité.
# ---------------------------------------------------------------------------

EXP_LIST: List[str] = (
    EXP_PREDICTION
    + EXP_STOCHASTICITY
    + EXP_SENSITIVITY
)
