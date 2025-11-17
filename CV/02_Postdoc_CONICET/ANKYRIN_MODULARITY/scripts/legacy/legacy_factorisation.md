Je vais te proposer une **vraie version modulaire moderne** du pipeline Ankyrin / binding partners, en Python 3, en restant au plus près de ce que tu avais codé.

Comme tout refaire *absolument tout* (tous les scripts de `all_python_scripts.txt`) dans un seul message est irréaliste, je me concentre ici sur le **cœur du pipeline “Pfam_local_search/binding_partners”** (conservation + enrichissement + mapping Pfam→clans), celui qui alimente ton poster/ draft.
Je te donne :

* des **modules complets** et utilisables (`config.py`, `io.py`, `clans.py`, `conservation.py`, `enrichment.py`, `pipeline.py`)
* un **fichier de mapping** qui relie tes anciens scripts aux nouvelles fonctions (`legacy_mapping.py`)
* un mini **README + schéma des dépendances**

Tu pourras ensuite me demander d’étendre la même logique aux autres scripts (ELMs, interacting pairs, etc.).

---

## 1. `config.py`

```python
# config.py
from pathlib import Path

# Racine du projet
PROJECT_ROOT = Path(__file__).resolve().parent

DATA_DIR = PROJECT_ROOT / "data"
RAW_DIR = DATA_DIR / "raw"
INTERMEDIATE_DIR = DATA_DIR / "intermediate"
RESULTS_DIR = DATA_DIR / "results"

# Fichiers d'entrée principaux
PFAM_CLANS_FILE = RAW_DIR / "Pfam-A.clans.txt"
PFAM_FREQUENCIES_FILE = RAW_DIR / "FrequenciesPfam_mapped.txt"

ANK_FASTA = RAW_DIR / "Ank_uniref50_string_mapped_522.fasta"
BD_FASTA = RAW_DIR / "binding_partners_2038.fasta"

# Fichiers de sortie "standards"
PFAM_CLAN_FREQ_FILE = INTERMEDIATE_DIR / "FrequenciesPfam_mapped_CLANS.txt"

# Paramètres par défaut
DEFAULT_MAX_HOMOLOGS = 400
```

---

## 2. `io.py`

Fonctions génériques d’I/O, parsing FASTA, et utilitaires.

```python
# io.py
import logging
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

from Bio import SeqIO  # type: ignore

logger = logging.getLogger(__name__)


def ensure_dir(path: Path) -> None:
    """Créer un répertoire (et parents) s'il n'existe pas."""
    path.mkdir(parents=True, exist_ok=True)


def read_uniprot_ids_from_fasta(fasta_path: Path) -> List[str]:
    """
    Reproduit ton pattern record.id[3:9] pour extraire l'UNP ID
    depuis un multifasta UniProt-style.
    """
    ids: List[str] = []
    with fasta_path.open("r") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            uniprot_id = record.id[3:9]
            ids.append(uniprot_id)
    logger.info("Read %d UniProt IDs from %s", len(ids), fasta_path)
    return ids


def read_pfam_table(path: Path) -> List[Tuple[str, int, int, str]]:
    """
    Lire un fichier PFAM tabulé de type:
    start\tend\t...\tdomain_id\t...

    Retourne une liste de tuples (start, end, domain_id, full_line_stripped).
    """
    entries: List[Tuple[str, int, int, str]] = []
    with path.open("r") as fh:
        for raw in fh:
            line = raw.rstrip("\n").split("\t")
            if len(line) < 6:
                continue
            start = int(line[1])
            end = int(line[2])
            domain_id = line[5]
            entries.append((domain_id, start, end, raw.rstrip("\n")))
    logger.debug("Read %d PFAM entries from %s", len(entries), path)
    return entries


def read_simple_two_column_int_table(path: Path) -> Dict[str, int]:
    """
    Lire un fichier tabulé simple: key\tvalue
    (pour counts, fréquences, etc.)
    """
    d: Dict[str, int] = {}
    with path.open("r") as fh:
        header = next(fh, None)
        for raw in fh:
            cols = raw.rstrip("\n").split("\t")
            if len(cols) < 2:
                continue
            key = cols[0]
            try:
                value = int(cols[1])
            except ValueError:
                logger.warning("Skipping line (non-int value): %s", raw.strip())
                continue
            d[key] = value
    logger.info("Read %d entries from %s", len(d), path)
    return d
```

---

## 3. `clans.py`

Regroupe la logique de `adapt_frequencies_pfam_to_clans.py` et `rename_domains_to_clans.py`.

```python
# clans.py
import logging
from pathlib import Path
from typing import Dict

from .config import PFAM_CLANS_FILE, PFAM_FREQUENCIES_FILE, PFAM_CLAN_FREQ_FILE
from .io import ensure_dir

logger = logging.getLogger(__name__)


def load_family_to_clan_mapping(clans_file: Path = PFAM_CLANS_FILE) -> Dict[str, str]:
    """
    Construit un dict PFAM-family → clan à partir de Pfam-A.clans.txt.
    Reste fidèle à ton parsing (colonnes 4 et 2). :contentReference[oaicite:1]{index=1}
    """
    mapping: Dict[str, str] = {}
    with clans_file.open("r") as fh:
        header = next(fh, None)
        for raw in fh:
            cols = raw.rstrip("\n").split("\t")
            if len(cols) < 4:
                continue
            clan = cols[1]
            family = cols[3]
            if clan != r"\N":
                mapping[family] = clan
    logger.info("Loaded %d PFAM family→clan mappings from %s", len(mapping), clans_file)
    return mapping


def aggregate_pfam_frequencies_by_clan(
    pfam_freq_file: Path = PFAM_FREQUENCIES_FILE,
    clans_file: Path = PFAM_CLANS_FILE,
    output_file: Path = PFAM_CLAN_FREQ_FILE,
) -> Path:
    """
    Reprend le comportement de adapt_frequencies_pfam_to_clans.py :
    - lit FrequenciesPfam_mapped.txt
    - regroupe par clan
    - écrit FrequenciesPfam_mapped_CLANS.txt :contentReference[oaicite:2]{index=2}
    """
    family_to_clan = load_family_to_clan_mapping(clans_file)
    clan_counts: Dict[str, int] = {}

    with pfam_freq_file.open("r") as fh:
        header = next(fh, None)
        for raw in fh:
            cols = raw.rstrip("\n").split("\t")
            if len(cols) < 2:
                continue
            domain = cols[1]
            try:
                count = int(cols[-1])
            except ValueError:
                logger.warning("Skipping line with non-int count: %s", raw.strip())
                continue

            clan = family_to_clan.get(domain, domain)
            clan_counts[clan] = clan_counts.get(clan, 0) + count

    ensure_dir(output_file.parent)
    with output_file.open("w") as out:
        out.write("Pfam_clan\tcount_in_swissprot\n")
        for clan, total in sorted(clan_counts.items()):
            out.write(f"{clan}\t{total}\n")

    logger.info("Wrote %d clan counts to %s", len(clan_counts), output_file)
    return output_file


def annotate_domains_with_clans(
    input_file: Path,
    clans_file: Path = PFAM_CLANS_FILE,
    output_file: Path | None = None,
) -> Path:
    """
    Equivalent à rename_domains_to_clans.py :
    ajoute une colonne CLAN à une table contenant un domaine PFAM en première colonne. :contentReference[oaicite:3]{index=3}
    """
    family_to_clan = load_family_to_clan_mapping(clans_file)
    if output_file is None:
        output_file = input_file.with_name(f"{input_file.stem}_with_clans.txt")

    ensure_dir(output_file.parent)

    with input_file.open("r") as fh, output_file.open("w") as out:
        lines = fh.readlines()
        if not lines:
            return output_file

        header_cols = lines[0].rstrip("\n").split("\t")
        out.write(f"{header_cols[0]}\tCLAN\t" + "\t".join(header_cols[1:]) + "\n")

        for raw in lines[1:]:
            cols = raw.rstrip("\n").split("\t")
            if not cols:
                continue
            domain_name = cols[0]
            clan_name = family_to_clan.get(domain_name, "NULL")
            out.write(f"{domain_name}\t{clan_name}\t" + "\t".join(cols[1:]) + "\n")

    logger.info("Annotated %s → %s with clan column", input_file, output_file)
    return output_file
```

---

## 4. `conservation.py`

Réécriture structurée de `check_domain_conservation_BD_2015.py` + calcul des pourcentages. 

```python
# conservation.py
import logging
from pathlib import Path
from collections import Counter, defaultdict
from typing import Dict, List, Tuple

from .config import INTERMEDIATE_DIR
from .io import ensure_dir, read_uniprot_ids_from_fasta

logger = logging.getLogger(__name__)


def _load_domain_occurrences(path: Path) -> List[Tuple[str, int, int]]:
    """
    Lit un fichier PFAM simplifié:
    uniprot_id \t start \t end \t domain_id
    (ce format correspond à tes pfam_domains_in_*_query.txt/homolog.txt).
    """
    occurrences: List[Tuple[str, int, int]] = []
    with path.open("r") as fh:
        for raw in fh:
            cols = raw.rstrip("\n").split("\t")
            if len(cols) < 4:
                continue
            unp = cols[0]
            start = int(cols[1])
            end = int(cols[2])
            domain = cols[3]
            occurrences.append((domain, start, end))
    return occurrences


def compute_domain_conservation_for_fasta(
    fasta_path: Path,
    num_homologs: int,
    base_dir: Path,
) -> Tuple[Path, Path]:
    """
    Equivalent structuré à check_domain_conservation_BD_2015.py pour un multifasta donné.
    - lit les uniprot IDs dans fasta_path
    - pour chaque uniprot, lit les fichiers query/homolog
    - calcule la conservation de chaque domaine
    - écrit:
      * conservation_state_...txt
      * domain_conservation_percentages_...txt

    On supppose des chemins du type:
    {basename}_MaxHomologs_MAX{N}_queries_pfam_output/pfam_domains_in_{UNP}_MaxHomologs_{N}_query.txt
    et mêmes patterns pour homologs (comme dans ton script original). :contentReference[oaicite:5]{index=5}
    """
    basename = fasta_path.stem  # ex: binding_partners_2038
    out_state = INTERMEDIATE_DIR / f"conservation_state_{basename}_{num_homologs}.txt"
    out_percent = INTERMEDIATE_DIR / f"domain_conservation_percentages_{basename}_{num_homologs}.txt"

    if out_percent.exists():
        logger.info("Conservation already computed for %s (N=%d)", basename, num_homologs)
        return out_state, out_percent

    ensure_dir(INTERMEDIATE_DIR)

    uniprot_ids = read_uniprot_ids_from_fasta(fasta_path)
    logger.info("Computing conservation for %d proteins in %s (N=%d)", len(uniprot_ids), fasta_path, num_homologs)

    # (domaine → [nb_conservé, nb_non_conservé])
    conservation_counts: Dict[str, List[int]] = defaultdict(lambda: [0, 0])

    with out_state.open("w") as state_out:
        state_out.write("uniprot_id\tdomain_id\tis_conserved\n")

        for unp in uniprot_ids:
            query_file = base_dir / f"{basename}_MaxHomologs_MAX{num_homologs}_queries_pfam_output" / \
                f"pfam_domains_in_{unp}_MaxHomologs_{num_homologs}_query.txt"
            homol_file = base_dir / f"{basename}_MaxHomologs_MAX{num_homologs}_homologs_pfam_output" / \
                f"pfam_domains_in_{unp}_MaxHomologs_{num_homologs}_homolog.txt"

            if not query_file.exists() or not homol_file.exists():
                logger.warning("Missing PFAM files for %s: %s or %s", unp, query_file, homol_file)
                continue

            # Collecte des domaines chez le query
            query_domains = _load_domain_occurrences(query_file)
            # Collecte des domaines homologues
            homolog_domains = _load_domain_occurrences(homol_file)
            homolog_counter = Counter([d for d, _, _ in homolog_domains])

            for domain_id, _, _ in query_domains:
                conserved = domain_id in homolog_counter
                state_out.write(f"{unp}\t{domain_id}\t{int(conserved)}\n")
                if conserved:
                    conservation_counts[domain_id][0] += 1
                else:
                    conservation_counts[domain_id][1] += 1

    # Calcul des pourcentages
    with out_percent.open("w") as perc_out:
        perc_out.write("domain_name\tconservation_percentage\n")
        for domain, (n_conserved, n_not) in sorted(conservation_counts.items()):
            total = n_conserved + n_not
            if total == 0:
                continue
            perc = n_conserved / total
            perc_out.write(f"{domain}\t{perc:.6f}\n")

    logger.info("Wrote conservation state to %s and percentages to %s", out_state, out_percent)
    return out_state, out_percent
```

> Remarque : on reproduit la logique de ton script (un fichier d’état + un fichier de pourcentages) mais de manière plus nette, avec logs et chemins paramétrables.

---

## 5. `enrichment.py`

Réécriture propre de la combinaison :

* calcul enrichissement log(obs/exp) pour BD (snippet avec `Pfam_domains_in_BD_2038.txt`) 
* colorisation par conservation (`domain_enrich_conserv_2015.py` et `domain_enrich_conserv_BD_2015.py`)

```python
# enrichment.py
import logging
import math
from pathlib import Path
from typing import Dict

from .config import INTERMEDIATE_DIR, RESULTS_DIR
from .io import ensure_dir, read_simple_two_column_int_table

logger = logging.getLogger(__name__)


def compute_domain_enrichment_from_counts(
    domain_counts_file: Path,
    pfam_counts_file: Path,
    num_proteins_subfamily: int,
    num_proteins_swissprot: int,
    output_file: Path,
) -> Path:
    """
    Version généralisée de ton bloc qui écrit Pfam_domains_in_BD_2038.txt :contentReference[oaicite:8]{index=8}

    domain_counts_file : table "domain_name <tab> # in subfamily"
    pfam_counts_file   : table "domain_name <tab> # in Uniprot"
    """
    if output_file.exists():
        logger.info("Enrichment file already exists: %s", output_file)
        return output_file

    ensure_dir(output_file.parent)

    sub_counts = read_simple_two_column_int_table(domain_counts_file)
    pfam_counts = read_simple_two_column_int_table(pfam_counts_file)

    with output_file.open("w") as out:
        out.write("domain_name\tcount_subfamily\tcount_uniprot\texpected_count\tlog_obs_exp\n")
        for domain_name, sub_count in sorted(sub_counts.items()):
            if domain_name not in pfam_counts:
                logger.debug("Domain %s not found in pfam_counts, skipping", domain_name)
                continue
            unp_count = pfam_counts[domain_name]
            exp_count = unp_count * num_proteins_subfamily / num_proteins_swissprot
            if exp_count == 0:
                continue
            obs_exp = sub_count / exp_count
            log_obs_exp = math.log(obs_exp)
            out.write(
                f"{domain_name}\t{sub_count}\t{unp_count}\t{exp_count:.6f}\t{log_obs_exp:.6f}\n"
            )

    logger.info("Wrote domain enrichment table to %s", output_file)
    return output_file


def _load_conservation_percentages(path: Path) -> Dict[str, float]:
    """
    Lit domain_conservation_percentages_*.txt.
    """
    d: Dict[str, float] = {}
    with path.open("r") as fh:
        header = next(fh, None)
        for raw in fh:
            cols = raw.rstrip("\n").split()
            if len(cols) < 2:
                continue
            try:
                d[cols[0]] = float(cols[1].replace(",", "."))
            except ValueError:
                continue
    return d


def color_enrichment_by_conservation(
    enrichment_zscore_file: Path,
    conservation_percent_file: Path,
    num_homologs: int,
    output_prefix: str,
) -> Path:
    """
    Version factorisée de domain_enrich_conserv_2015.py et _BD_2015.py 

    enrichment_zscore_file :
        - Pfam_domains_in_Ank1234_Zscores.txt
        - ou Pfam_domains_in_BD_2038_Zscores.txt
      colonnes (au moins) : domain_name, ..., log(obs/exp), Zscore

    conservation_percent_file :
        - domain_conservation_percentages_*.txt

    output_prefix :
        - ex: "Ank1234" ou "BD_2038"
    """
    ensure_dir(RESULTS_DIR)
    output_file = RESULTS_DIR / f"Pfam_domains_in_{output_prefix}_MaxHomologs_{num_homologs}_Zscores_color.txt"

    if output_file.exists():
        logger.info("Colored enrichment file already exists: %s", output_file)
        return output_file

    conserv = _load_conservation_percentages(conservation_percent_file)
    if not conserv:
        logger.warning("Empty conservation file: %s", conservation_percent_file)

    avg_conserv = sum(conserv.values()) / len(conserv) if conserv else 0.0
    logger.info(
        "Average domain conservation (N=%d) = %.4f from %s",
        len(conserv),
        avg_conserv,
        conservation_percent_file,
    )

    with enrichment_zscore_file.open("r") as fh, output_file.open("w") as out:
        header = next(fh, None)
        out.write(
            f"domain_name\tdomain_enrichment_log(obs/exp)\tZscore_enrichment\t"
            f"domain_conservation_over_{num_homologs}Homologs\tcolor\n"
        )

        for raw in fh:
            cols = raw.rstrip("\n").split("\t")
            if len(cols) < 6:
                continue
            domain_name = cols[0]
            try:
                domain_enrichment = float(cols[4].replace(",", "."))
                domain_zscore = float(cols[5].replace(",", "."))
            except ValueError:
                continue

            if domain_name not in conserv:
                logger.debug("Domain %s not in conservation dict", domain_name)
                continue

            domain_conservation = conserv[domain_name]
            color = '"red"' if (domain_zscore > 1 and domain_conservation > avg_conserv) else '"black"'
            out.write(
                f"{domain_name}\t{domain_enrichment:.6f}\t{domain_zscore:.6f}\t"
                f"{domain_conservation:.6f}\t{color}\n"
            )

    logger.info("Wrote colored enrichment table to %s", output_file)
    return output_file
```

---

## 6. `pipeline.py`

Script exécutable, avec **contrôle manuel** par flags, **logging**, et **reprise via fichiers intermédiaires**.

```python
# pipeline.py
import logging
from pathlib import Path

from .config import (
    PROJECT_ROOT,
    RAW_DIR,
    INTERMEDIATE_DIR,
    RESULTS_DIR,
    PFAM_FREQUENCIES_FILE,
    PFAM_CLANS_FILE,
    PFAM_CLAN_FREQ_FILE,
    ANK_FASTA,
    BD_FASTA,
    DEFAULT_MAX_HOMOLOGS,
)
from .io import ensure_dir
from .clans import aggregate_pfam_frequencies_by_clan, annotate_domains_with_clans
from .conservation import compute_domain_conservation_for_fasta
from .enrichment import color_enrichment_by_conservation

# -------------------------
# Configuration des étapes
# -------------------------
RUN_CLAN_AGGREGATION = True
RUN_CONSERVATION_ANK = False
RUN_CONSERVATION_BD = False
RUN_ENRICHMENT_ANK = False
RUN_ENRICHMENT_BD = False
# -------------------------


def setup_logging() -> None:
    log_file = PROJECT_ROOT / "ankyrin_modularity.log"
    logging.basicConfig(
        filename=str(log_file),
        level=logging.INFO,
        format="%(asctime)s — %(name)s — %(levelname)s — %(message)s",
    )
    # logger console minimal
    console = logging.StreamHandler()
    console.setLevel(logging.INFO)
    formatter = logging.Formatter("%(name)s — %(levelname)s — %(message)s")
    console.setFormatter(formatter)
    logging.getLogger().addHandler(console)


def main() -> None:
    setup_logging()
    logger = logging.getLogger(__name__)
    ensure_dir(RAW_DIR)
    ensure_dir(INTERMEDIATE_DIR)
    ensure_dir(RESULTS_DIR)

    logger.info("Starting Ankyrin modularity pipeline")

    # 1) Agrégation Pfam→clans
    if RUN_CLAN_AGGREGATION:
        logger.info("Step: aggregate PFAM frequencies by clan")
        aggregate_pfam_frequencies_by_clan(
            pfam_freq_file=PFAM_FREQUENCIES_FILE,
            clans_file=PFAM_CLANS_FILE,
            output_file=PFAM_CLAN_FREQ_FILE,
        )

    # 2) Conservation Ankyrin (si tu as les pfam_domains_in_* pour ANK)
    if RUN_CONSERVATION_ANK:
        logger.info("Step: compute conservation for Ankyrin family")
        compute_domain_conservation_for_fasta(
            fasta_path=ANK_FASTA,
            num_homologs=DEFAULT_MAX_HOMOLOGS,
            base_dir=INTERMEDIATE_DIR,  # à adapter si nécessaire
        )

    # 3) Conservation Binding Partners
    if RUN_CONSERVATION_BD:
        logger.info("Step: compute conservation for Binding Partners")
        compute_domain_conservation_for_fasta(
            fasta_path=BD_FASTA,
            num_homologs=DEFAULT_MAX_HOMOLOGS,
            base_dir=INTERMEDIATE_DIR,
        )

    # 4) Enrichissement + couleur pour Ank
    if RUN_ENRICHMENT_ANK:
        logger.info("Step: enrichment + conservation coloring for Ankyrin")
        conservation_file = INTERMEDIATE_DIR / f"domain_conservation_percentages_{ANK_FASTA.stem}_{DEFAULT_MAX_HOMOLOGS}.txt"
        zscore_file = RAW_DIR / "Pfam_domains_in_Ank1234_Zscores.txt"  # même nom que ton fichier original :contentReference[oaicite:10]{index=10}
        color_enrichment_by_conservation(
            enrichment_zscore_file=zscore_file,
            conservation_percent_file=conservation_file,
            num_homologs=DEFAULT_MAX_HOMOLOGS,
            output_prefix="Ank1234",
        )

    # 5) Enrichissement + couleur pour BD
    if RUN_ENRICHMENT_BD:
        logger.info("Step: enrichment + conservation coloring for Binding Partners")
        conservation_file = INTERMEDIATE_DIR / f"domain_conservation_percentages_{BD_FASTA.stem}_{DEFAULT_MAX_HOMOLOGS}.txt"
        zscore_file = RAW_DIR / "Pfam_domains_in_BD_2038_Zscores.txt"  # même nom que ton fichier original :contentReference[oaicite:11]{index=11}
        color_enrichment_by_conservation(
            enrichment_zscore_file=zscore_file,
            conservation_percent_file=conservation_file,
            num_homologs=DEFAULT_MAX_HOMOLOGS,
            output_prefix="binding_partners_2038.fasta",
        )

    logger.info("Pipeline finished.")


if __name__ == "__main__":
    main()
```

Tu peux contrôler l’exécution des étapes uniquement en changeant les booleans au début du fichier, au lieu de commenter/décommenter des blocs de code partout.

---

## 7. `legacy_mapping.py`

Mapping entre tes anciens scripts et les nouvelles fonctions (uniformisation demandée).

```python
# legacy_mapping.py
"""
Mapping entre les anciens scripts "one-shot" et les fonctions du pipeline modulaire.
"""

LEGACY_SCRIPT_TO_FUNCTION = {
    # CLANS / PFAM
    "adapt_frequencies_pfam_to_clans.py": "clans.aggregate_pfam_frequencies_by_clan",
    "rename_domains_to_clans.py": "clans.annotate_domains_with_clans",

    # CONSERVATION
    "check_domain_conservation_BD_2015.py": "conservation.compute_domain_conservation_for_fasta",
    # (pour Ankyrin : même fonction mais avec ANK_FASTA)

    # ENRICHISSEMENT + COULEUR
    "domain_enrich_conserv_2015.py": "enrichment.color_enrichment_by_conservation (version Ank)",
    "domain_enrich_conserv_BD_2015.py": "enrichment.color_enrichment_by_conservation (version BD)",

    # COMPTEURS / PARSING (patterns, à implémenter si besoin)
    "parse_and_count_pfam_BD_bis.py": "io / parsing utilitaire pour counts BD",
    "count_pfam_domains_BD.py": "enrichment.compute_domain_enrichment_from_counts",
    "concaten_pfam_outputs_BD_2015.py": "io.concatenate_pfam_outputs_by_threshold (à factoriser)",
    "make_fasta_BD.py": "io / extraction de séquences domaine/homologue (à factoriser)",
    "make_indep_fasta_BD_XXX.py": "io / extraction de séquences indépendantes (à factoriser)",

    # CHECK / UTILITAIRES
    "check_unp_id.py": "io.read_uniprot_ids_from_fasta + comparaison de listes",
    "general_naming.py": "comparaison de listes de noms (peut devenir une fonction utils.compare_name_lists)",
}
```

---

## 8. Mini README + schéma de dépendances

Tu peux mettre ceci dans `README.md` à la racine du projet ankyrin_modularity :

````markdown
# Ankyrin modularity – Domain & Conservation pipeline (refactor)

Ce dossier contient une version modulaire et moderne du pipeline utilisé pour :

- l’enrichissement en domaines Pfam / clans dans les Ankyrins et leurs partenaires,
- la mesure de conservation des domaines dans les homologues,
- la combinaison enrichissement + conservation (log(obs/exp) + Z-scores + couleur),
- en cohérence avec le poster ISMB 2014 et le draft d’article.

## Structure

```text
ankyrin_modularity/
├── config.py
├── io.py
├── clans.py
├── conservation.py
├── enrichment.py
├── pipeline.py
├── legacy_mapping.py
└── data/
    ├── raw/
    │   ├── Ank_uniref50_string_mapped_522.fasta
    │   ├── binding_partners_2038.fasta
    │   ├── Pfam-A.clans.txt
    │   ├── FrequenciesPfam_mapped.txt
    │   ├── Pfam_domains_in_Ank1234_Zscores.txt
    │   └── Pfam_domains_in_BD_2038_Zscores.txt
    ├── intermediate/
    └── results/
````

## Dépendances entre modules

```text
          +-----------+
          | config.py |
          +-----+-----+
                |
      +---------+----------+
      |                    |
   +--v--+            +----v-----+
   |io.py|            |clans.py  |
   +--+--+            +----+-----+
      |                    |
      |                    |
   +--v-----------+   +----v----------+
   |conservation.py|   |enrichment.py |
   +-------+------+   +------+--------+
           \             /
            \           /
             +----v----+
             |pipeline |
             +---------+
```

* `io.py` fournit les fonctions génériques de parsing (FASTA, tables).
* `clans.py` gère le mapping domaine→clan et l’agrégation des fréquences.
* `conservation.py` calcule la conservation des domaines sur les homologues.
* `enrichment.py` combine enrichissement & conservation pour générer les fichiers finalisés (avec couleur).
* `pipeline.py` orchestre tout ça avec des flags, du logging, et la reprise sur fichiers intermédiaires.
* `legacy_mapping.py` documente la correspondance avec les anciens scripts one-shot.

```

---

Si tu veux, on peut ensuite :

- ajouter une fonction `concatenate_pfam_outputs_by_threshold` dans `io.py` qui encapsule proprement `concaten_pfam_outputs_BD_2015.py`,  
- factoriser `make_fasta_BD.py` / `make_indep_fasta_BD_XXX.py` dans des fonctions d’extraction de séquences domaine/homologue,  
- et étendre la même logique modulaire à la partie ELM / interacting_pairs.

---

**Q1 : Veux-tu que la prochaine étape soit de factoriser `make_fasta_BD.py` et `make_indep_fasta_BD_XXX.py` dans `io.py` (extraction de séquences domaines/homologues) ?**  
**Q2 : Tu préfères garder les chemins “hardcodés” vers tes dossiers d’alignements (`/home/nina/...`) dans une section spéciale de `config.py`, ou qu’on les paramètre par CLI (arguments) ?**  
**Q3 : Souhaites-tu que je te propose ensuite une version “portfolio” de ce projet (README + petit diagramme pipeline + extrait de code) prête à être mise sur GitHub ?**
```
