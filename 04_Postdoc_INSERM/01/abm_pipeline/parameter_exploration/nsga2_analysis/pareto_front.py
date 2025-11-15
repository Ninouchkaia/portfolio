# abm_pipeline/parameter_exploration/nsga2_analysis/pareto_front.py

import math
from pathlib import Path
from typing import List, Tuple

import matplotlib.pyplot as plt

from abm_pipeline.parameter_exploration.utils import logger, debug


def _compute_all_points(lines: List[str]) -> List[List[float]]:
    points: List[List[float]] = []

    for line in lines[1:]:
        line = line.replace("\n", "").split(",")
        apo = int(line[0])
        need_sig = int(line[1])
        layers = int(line[2])
        alpha = int(line[3])
        mono_phago = int(line[4])
        nlc_phago = int(line[5])
        m2_phago = int(line[6])
        m2_kill = int(line[7])
        cll_dist = int(line[8])
        mono_dist = int(line[9])
        nlc_dist = int(line[10])
        macro_dist = int(line[11])
        nlc_threshold = int(line[12])
        signal_init_mean = int(line[13])
        signal_init_std = int(line[14])
        diff_time = int(line[15])
        diff_init_std = int(line[16])
        gamma_life_init = int(line[17])
        alpha_distrib = float(line[18])
        delta_via = float(line[19])
        delta_conc = float(line[20])

        euclid = math.sqrt(delta_via ** 2 + delta_conc ** 2)

        points.append(
            [
                delta_via,
                delta_conc,
                euclid,
                1,  # marqueur Pareto ou non
                apo,
                need_sig,
                layers,
                alpha,
                mono_phago,
                nlc_phago,
                m2_phago,
                m2_kill,
                cll_dist,
                mono_dist,
                nlc_dist,
                macro_dist,
                nlc_threshold,
                signal_init_mean,
                signal_init_std,
                diff_time,
                diff_init_std,
                gamma_life_init,
                alpha_distrib,
            ]
        )
    return points


def _flag_pareto_points(points: List[List[float]]):
    # déduplication
    unique = list({tuple(p) for p in points})
    unique = [list(p) for p in unique]

    unique.sort(key=lambda x: x[2])  # par distance euclidienne

    for i in range(len(unique)):
        x1, y1 = unique[i][0], unique[i][1]
        for j in range(len(unique)):
            x2, y2 = unique[j][0], unique[j][1]
            if (x1, y1) != (x2, y2) and (x2 <= x1) and (y2 <= y1):
                unique[i][3] = 0
                break

    pareto_front = [
        (p[0], p[1], p[4:])
        for p in unique
        if p[3] == 1
    ]
    return pareto_front, unique


def build_pareto_front(input_file: str, output_prefix: str) -> None:
    """
    Construit le Pareto front à partir d'un CSV OpenMOLE.

    - input_file: CSV avec colonnes [apo, needSig, ..., delta_fitness_via, delta_fitness_conc]
    - output_prefix: préfixe pour PNG et TXT ('pareto_ABM_2D_patient' par ex.)
    """
    input_path = Path(input_file)
    logger.info(f"Building pareto front from {input_path}")

    with input_path.open("r", encoding="utf-8") as file_read:
        lines = file_read.readlines()

    points = _compute_all_points(lines)
    pareto_front, all_points = _flag_pareto_points(points)

    # Plot
    x_val = [p[0] for p in all_points]
    y_val = [p[1] for p in all_points]
    x_pareto = [p[0] for p in pareto_front]
    y_pareto = [p[1] for p in pareto_front]

    fig, ax = plt.subplots()
    ax.scatter(x_val, y_val, s=3)
    ax.scatter(x_pareto, y_pareto, s=5)
    ax.set_xlabel(r"$\Delta Viability$ fitness")
    ax.set_ylabel(r"$\Delta Concentration$ fitness")
    ax.set_title("Pareto front")
    ax.grid(True)
    fig.tight_layout()
    png_path = Path(f"{output_prefix}.png")
    fig.savefig(png_path, bbox_inches="tight")

    logger.info(f"len(pareto_front) = {len(pareto_front)}")

    pareto_sorted = sorted(pareto_front, key=lambda x: x[0])
    header = (
        "delta_fitness_via,delta_fitness_conc,apo,needSig,layers,alpha,"
        "monoPhago,NLCPhago,M2Phago,M2Kill,cllDist,MonoDist,nlcDist,"
        "macroDist,nlcThreshold,signalInitMean,signalInitStd,diffTime,"
        "diffInitStd,LifeInitGamma,alphaDistrib\n"
    )

    txt_path = Path(f"{output_prefix}.txt")
    with txt_path.open("w", encoding="utf-8") as file_write:
        file_write.write(header)
        for sets in pareto_sorted:
            line = (",".join(str(x) for x in sets)).replace("[", "").replace("]", "")
            file_write.write(line + "\n")

    logger.info(f"Pareto text written to {txt_path}")
