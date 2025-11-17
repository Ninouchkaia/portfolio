from pathlib import Path
from .io import read_tsv, write_lines

MODEL_TF = [
    "STAT1", "STAT5A", "STAT5B", "NFKB1", "NFKB2",
    "PPARG", "STAT6", "STAT3", "IRF3", "IRF5",
    "IRF4", "KLF4", "HIF1A",
]


def compute_signatures(input_file: Path, output_file: Path):
    rows = read_tsv(input_file)

    M1, M2, NLC = [], [], []
    output_lines = []

    for row in rows:
        TF = row[0]
        m1 = float(row[2])
        m2 = float(row[3])
        nlc = float(row[4])

        if nlc > m1 and nlc > m2:
            NLC.append(TF)
            output_lines.append(f"{TF}\t{m1}\t{m2}\t{nlc}")
        elif m1 > nlc and m1 > m2:
            M1.append(TF)
        elif m2 > m1 and m2 > nlc:
            M2.append(TF)

    write_lines(output_file, output_lines)

    return {
        "M1_signature": set(M1),
        "M2_signature": set(M2),
        "NLC_signature": set(NLC),
    }
