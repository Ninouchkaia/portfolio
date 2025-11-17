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
    PFAM_QUERY_DIR,
    PFAM_HOMOLOG_DIR,
    DOMAIN_FASTA_DIR,
    INDEP_FASTA_DIR,
)

from .io import (
    ensure_dir,
    concatenate_pfam_outputs_by_threshold,
    find_pfam_files,
    extract_domain_sequences,
    extract_independent_sequences,
    read_uniprot_ids_from_fasta,
)

from .clans import aggregate_pfam_frequencies_by_clan
from .conservation import compute_domain_conservation_for_fasta
from .enrichment import color_enrichment_by_conservation

from .elm import (
    read_elm_table,
    count_elms,
    write_elm_counts,
    compute_elm_enrichment,
    compute_elm_domain_cooccurrences,
)

# ---------------------------------------------------------
# ACTIVATION DES ÉTAPES (contrôle manuel)
# ---------------------------------------------------------
RUN_CLAN_AGGREGATION = False
RUN_PFAM_CONCATENATION = False
RUN_DOMAIN_EXTRACTION = False
RUN_INDEP_EXTRACTION = False

RUN_CONSERVATION_ANK = False
RUN_CONSERVATION_BD = False

RUN_ENRICHMENT_ANK = False
RUN_ENRICHMENT_BD = False

RUN_ELM_PROCESSING = False
RUN_ELM_ENRICHMENT = False
RUN_ELM_DOMAIN_ASSOCIATION = False
# ---------------------------------------------------------


def setup_logging() -> None:
    log_file = PROJECT_ROOT / "ankyrin_modularity.log"
    logging.basicConfig(
        filename=str(log_file),
        level=logging.INFO,
        format="%(asctime)s — %(name)s — %(levelname)s — %(message)s",
    )
    console = logging.StreamHandler()
    console.setLevel(logging.INFO)
    formatter = logging.Formatter("%(name)s — %(levelname)s — %(message)s")
    console.setFormatter(formatter)
    logging.getLogger().addHandler(console)


def main() -> None:

    setup_logging()
    logger = logging.getLogger(__name__)

    logger.info("=== Starting Ankyrin modularity pipeline ===")

    # Assure la présence des répertoires
    ensure_dir(RAW_DIR)
    ensure_dir(INTERMEDIATE_DIR)
    ensure_dir(RESULTS_DIR)
    ensure_dir(PFAM_QUERY_DIR)
    ensure_dir(PFAM_HOMOLOG_DIR)
    ensure_dir(DOMAIN_FASTA_DIR)
    ensure_dir(INDEP_FASTA_DIR)

    # ---------------------------------------------------------
    # 1) PFAM → CLANS
    # ---------------------------------------------------------
    if RUN_CLAN_AGGREGATION:
        logger.info("Step: PFAM → Clan aggregation")
        aggregate_pfam_frequencies_by_clan(
            pfam_freq_file=PFAM_FREQUENCIES_FILE,
            clans_file=PFAM_CLANS_FILE,
            output_file=PFAM_CLAN_FREQ_FILE,
        )

    # ---------------------------------------------------------
    # 2) PFAM → Concaténation tous fichiers + seuil
    # ---------------------------------------------------------
    if RUN_PFAM_CONCATENATION:
        logger.info("Step: Concatenate PFAM outputs (query files)")
        concatenated = concatenate_pfam_outputs_by_threshold(
            input_dir=PFAM_QUERY_DIR,
            output_file=INTERMEDIATE_DIR / "PFAM_domains_concatenated.txt",
            min_occurrences=1,
        )
        logger.info("Concatenated PFAM file saved to %s", concatenated)

    # ---------------------------------------------------------
    # 3) EXTRACTION DE DOMAINES (FASTA)
    # ---------------------------------------------------------
    if RUN_DOMAIN_EXTRACTION:
        logger.info("Step: Extract domain FASTA for Binding Partners")

        uniprot_list = read_uniprot_ids_from_fasta(BD_FASTA)

        for unp in uniprot_list:
            query_file, homolog_file = find_pfam_files(unp, DEFAULT_MAX_HOMOLOGS)

            if query_file.exists():
                extract_domain_sequences(
                    fasta_path=BD_FASTA,
                    pfam_domain_file=query_file,
                )

    # ---------------------------------------------------------
    # 4) EXTRACTION DE SÉQUENCES INDÉPENDANTES
    # ---------------------------------------------------------
    if RUN_INDEP_EXTRACTION:
        logger.info("Step: Extract independent sequences")
        extract_independent_sequences(
            fasta_path=BD_FASTA,
            cutoff_start=150,
            cutoff_end=350,
        )

    # ---------------------------------------------------------
    # 5) CONSERVATION Ankyrin
    # ---------------------------------------------------------
    if RUN_CONSERVATION_ANK:
        logger.info("Step: Compute conservation for Ankyrin family")
        compute_domain_conservation_for_fasta(
            fasta_path=ANK_FASTA,
            num_homologs=DEFAULT_MAX_HOMOLOGS,
            base_dir=INTERMEDIATE_DIR,
        )

    # ---------------------------------------------------------
    # 6) CONSERVATION Binding Partners
    # ---------------------------------------------------------
    if RUN_CONSERVATION_BD:
        logger.info("Step: Compute conservation for Binding Partners")
        compute_domain_conservation_for_fasta(
            fasta_path=BD_FASTA,
            num_homologs=DEFAULT_MAX_HOMOLOGS,
            base_dir=INTERMEDIATE_DIR,
        )

    # ---------------------------------------------------------
    # 7) ENRICHISSEMENT + COULEUR Ankyrin
    # ---------------------------------------------------------
    if RUN_ENRICHMENT_ANK:
        logger.info("Step: Enrichment + color for Ankyrins")

        conservation_file = (
            INTERMEDIATE_DIR / f"domain_conservation_percentages_{ANK_FASTA.stem}_{DEFAULT_MAX_HOMOLOGS}.txt"
        )
        zscore_file = RAW_DIR / "Pfam_domains_in_Ank1234_Zscores.txt"

        color_enrichment_by_conservation(
            enrichment_zscore_file=zscore_file,
            conservation_percent_file=conservation_file,
            num_homologs=DEFAULT_MAX_HOMOLOGS,
            output_prefix="Ank1234",
        )

    # ---------------------------------------------------------
    # 8) ENRICHISSEMENT + COULEUR Binding Partners
    # ---------------------------------------------------------
    if RUN_ENRICHMENT_BD:
        logger.info("Step: Enrichment + color for Binding Partners")

        conservation_file = (
            INTERMEDIATE_DIR / f"domain_conservation_percentages_{BD_FASTA.stem}_{DEFAULT_MAX_HOMOLOGS}.txt"
        )
        zscore_file = RAW_DIR / "Pfam_domains_in_BD_2038_Zscores.txt"

        color_enrichment_by_conservation(
            enrichment_zscore_file=zscore_file,
            conservation_percent_file=conservation_file,
            num_homologs=DEFAULT_MAX_HOMOLOGS,
            output_prefix="BD_2038",
        )

    # ---------------------------------------------------------
    # 9) SECTION ELM (SLiMs)
    # ---------------------------------------------------------
    if RUN_ELM_PROCESSING:
        logger.info("Step: Counting ELMs")

        elm_file_ank = RAW_DIR / "ELM_predictions_ankyrins.txt"
        elm_file_bd = RAW_DIR / "ELM_predictions_binding_partners.txt"

        ank_elms = read_elm_table(elm_file_ank)
        bd_elms = read_elm_table(elm_file_bd)

        ank_counts = count_elms(ank_elms)
        bd_counts = count_elms(bd_elms)

        write_elm_counts(ank_counts, RESULTS_DIR / "ELM_counts_ankyrins.txt")
        write_elm_counts(bd_counts, RESULTS_DIR / "ELM_counts_binding_partners.txt")

    if RUN_ELM_ENRICHMENT:
        logger.info("Step: ELM enrichment")

        compute_elm_enrichment(
            elm_counts_subfamily=ank_counts,
            elm_counts_background=bd_counts,
            num_subfamily_proteins=len(ank_counts),
            num_background_proteins=len(bd_counts),
            output_file=RESULTS_DIR / "ELM_enrichment_ank_vs_bd.txt",
        )

    if RUN_ELM_DOMAIN_ASSOCIATION:
        logger.info("Step: ELM-domain cooccurrence")

        # Domain table BD format: uniprot\t domain
        domain_table_file = RESULTS_DIR / "domains_binding_partners_for_elm.txt"
        domain_table = []
        with domain_table_file.open() as fh:
            next(fh)
            for raw in fh:
                cols = raw.rstrip("\n").split("\t")
                if len(cols) >= 2:
                    domain_table.append((cols[0], cols[1]))

        compute_elm_domain_cooccurrences(
            elm_table=bd_elms,
            domain_table=domain_table,
            output_file=RESULTS_DIR / "ELM_domain_cooccurrences.txt",
        )

    # ---------------------------------------------------------
    logger.info("=== Pipeline finished ===")


if __name__ == "__main__":
    main()
