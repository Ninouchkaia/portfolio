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
