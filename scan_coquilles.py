#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Scan automatique des tournures interdites :
- toute référence à "tu", "te", "tes"
- toute tournure possible de collaboration :
  "basé sur ce que tu avais écrit", "adapté depuis", "selon ton fichier", etc.
- toute référence implicite à un dialogue

Sortie : report_coquilles.txt
"""

import os
import re

# mots/tournures que tu veux éliminer :
PATTERNS = [
    r"\btu\b",
    r"\bte\b",
    r"\btes\b",
    r"bas[ée] sur ce que.*tu",
    r"inspir[ée] de ce que.*tu",
    r"selon ce que.*tu",
    r"comme tu .*",
    r"adapt[ée] depuis",
    r"version inspir[ée]",
    r"selon ton",
    r"draft initial.*utili",  # si tu veux le garder, enlève-le
    r"\byou\b",
    r"\byour\b",
    r"\byours\b",
]

EXTENSIONS = {".py", ".md", ".txt", ".rst", ".yml", ".yaml", ".json"}

OUTPUT = "report_coquilles.txt"


def scan_file(path, patterns):
    try:
        with open(path, "r", encoding="utf-8", errors="ignore") as f:
            lines = f.readlines()
    except Exception:
        return []

    findings = []
    for lineno, line in enumerate(lines, start=1):
        for pattern in patterns:
            if re.search(pattern, line, flags=re.IGNORECASE):
                findings.append((lineno, line.strip(), pattern))
    return findings


def main():
    report = []

    for root, _, files in os.walk("."):
        for f in files:
            _, ext = os.path.splitext(f)
            if ext.lower() not in EXTENSIONS:
                continue

            path = os.path.join(root, f)
            findings = scan_file(path, PATTERNS)
            if findings:
                report.append(f"\n=== {path} ===")
                for lineno, text, pattern in findings:
                    report.append(
                        f"  L{lineno}: {text}\n     → match: /{pattern}/"
                    )

    with open(OUTPUT, "w", encoding="utf-8") as out:
        if report:
            out.write("\n".join(report))
        else:
            out.write("Aucune occurrence trouvée.")

    print(f"[DONE] Rapport écrit dans {OUTPUT}")


if __name__ == "__main__":
    main()
