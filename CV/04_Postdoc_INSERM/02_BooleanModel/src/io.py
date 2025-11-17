from pathlib import Path


def read_tsv(filepath: Path):
    rows = []
    with open(filepath, "r", encoding="utf-8") as f:
        for line in f.readlines()[1:]:
            parts = line.strip().split("\t")
            rows.append(parts)
    return rows


def write_lines(filepath: Path, lines):
    with open(filepath, "a+", encoding="utf-8") as f:
        for line in lines:
            f.write(line + "\n")
