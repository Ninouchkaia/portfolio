"""
Amino Acid Usage Analysis Pipeline (PaxDB-based)
Author: Nina Verstraete (refactored)

This module orchestrates the complete workflow used in the
"Evolutionary Constraints on Amino Acid Usage in Proteomes" project
(INQUIMAE–CONICET, Buenos Aires, 2013–2015).

The pipeline follows these major steps:

    1. Load PaxDB proteomes (abundance + FASTA sequences)
    2. Map PaxDB protein IDs to UniProt IDs
    3. Compute unweighted AA frequencies (DS1-like)
    4. Compute abundance-weighted AA frequencies (DS2-like)
    5. Compute model expectations:
         - Genetic code (GC-content dependent)
         - Metabolic cost model
    6. Compute correlations between observed and model predictions
    7. Generate paper-like figures
"""

from __future__ import annotations
import pathlib
from typing import Sequence, Dict

# Internal modules (to be implemented based on your old scripts)
from . import io, mapping, counting, models, plots


# ---------------------------------------------------------------------
# CONFIGURATION OBJECT
# ---------------------------------------------------------------------

class AminoAcidPipelineConfig:
    """
    Central configuration for defining input and output paths.
    """

    def __init__(
        self,
        raw_dir: str | pathlib.Path,
        processed_dir: str | pathlib.Path,
        organisms: Sequence[str],
        abundance_suffix: str = "_abundance.tsv",
        fasta_suffix: str = "_proteome.fasta",
    ):
        self.raw_dir = pathlib.Path(raw_dir)
        self.processed_dir = pathlib.Path(processed_dir)
        self.organisms = list(organisms)
        self.abundance_suffix = abundance_suffix
        self.fasta_suffix = fasta_suffix

        self.processed_dir.mkdir(parents=True, exist_ok=True)

    def abundance_path(self, org: str) -> pathlib.Path:
        return self.raw_dir / f"{org}{self.abundance_suffix}"

    def fasta_path(self, org: str) -> pathlib.Path:
        return self.raw_dir / f"{org}{self.fasta_suffix}"

    def mapped_proteome_path(self, org: str) -> pathlib.Path:
        return self.processed_dir / f"{org}_proteome_mapped.tsv"

    def aa_freq_unweighted_path(self, org: str) -> pathlib.Path:
        return self.processed_dir / f"{org}_aa_freq_unweighted.tsv"

    def aa_freq_weighted_path(self, org: str) -> pathlib.Path:
        return self.processed_dir / f"{org}_aa_freq_weighted.tsv"

    def merged_aa_freq_unweighted(self) -> pathlib.Path:
        return self.processed_dir / "DS1_aa_freq_unweighted.tsv"

    def merged_aa_freq_weighted(self) -> pathlib.Path:
        return self.processed_dir / "DS2_aa_freq_weighted.tsv"


# ---------------------------------------------------------------------
# STEP 1 — PREPARE PROTEOMES
# ---------------------------------------------------------------------

def prepare_proteomes(cfg: AminoAcidPipelineConfig) -> None:
    """
    Builds unified proteome tables for each organism:

        protein_id, uniprot_id, sequence, abundance

    This step merges the logic of:
        - parsing_fasta_file*.py
        - retrieve_uniprot_id.py / map_uniprot_id.py
        - cumul_mapping*.py
        - coverage_summary*.py
    """

    for org in cfg.organisms:
        print(f"[1/4] Preparing proteome for {org}…")

        abund_path = cfg.abundance_path(org)
        fasta_path = cfg.fasta_path(org)

        abund_df = io.load_abundance_table(abund_path)
        seq_dict = io.load_fasta_sequences(fasta_path)

        mapped_df = mapping.map_ids_and_merge(abund_df, seq_dict)

        cov_stats = mapping.compute_coverage_stats(mapped_df, abund_df)
        mapping.save_coverage_stats(
            cov_stats, cfg.processed_dir / f"{org}_coverage.txt"
        )

        io.save_proteome_table(mapped_df, cfg.mapped_proteome_path(org))


# ---------------------------------------------------------------------
# STEP 2 — COMPUTE AA FREQUENCIES (DS1 & DS2)
# ---------------------------------------------------------------------

def compute_amino_acid_frequencies(cfg: AminoAcidPipelineConfig) -> None:
    """
    Computes amino acid frequencies for each organism:

       - DS1-like: unweighted (all proteins equal)
       - DS2-like: weighted by PaxDB abundance

    This merges the logic of:
        amino_acid_count*.py
        amino_acid_count_aux.py
    """

    all_unweighted = []
    all_weighted = []

    for org in cfg.organisms:
        print(f"[2/4] Counting amino acids for {org}…")

        proteome_path = cfg.mapped_proteome_path(org)
        proteome = io.load_proteome_table(proteome_path)

        aa_unw = counting.amino_acid_frequencies_unweighted(proteome, org)
        aa_w = counting.amino_acid_frequencies_weighted(proteome, org)

        io.save_table(aa_unw, cfg.aa_freq_unweighted_path(org))
        io.save_table(aa_w, cfg.aa_freq_weighted_path(org))

        all_unweighted.append(aa_unw)
          all_weighted.append(aa_w)

    merged_unweighted = counting.merge_aa_tables(all_unweighted)
    merged_weighted = counting.merge_aa_tables(all_weighted)

    io.save_table(merged_unweighted, cfg.merged_aa_freq_unweighted())
    io.save_table(merged_weighted, cfg.merged_aa_freq_weighted())


# ---------------------------------------------------------------------
# STEP 3 — MODELS & CORRELATIONS
# ---------------------------------------------------------------------

def compute_models_and_correlations(
    cfg: AminoAcidPipelineConfig,
    gc_table_path: str | pathlib.Path
) -> None:
    """
    Computes:
        - Genetic code model (GC → expected AA frequencies)
        - Metabolic cost model (ATP/time)
        - Correlations with DS1/DS2 observed values

    Corresponds to paper Figures 1–2.
    """

    print("[3/4] Loading DS1 / DS2…")

    ds1 = io.load_table(cfg.merged_aa_freq_unweighted())
    ds2 = io.load_table(cfg.merged_aa_freq_weighted())
    gc_df = io.load_table(gc_table_path)

    print("[3/4] Computing genetic code model expectations…")

    expected_genetic_ds1 = models.genetic_code_model(ds1, gc_df)
    expected_genetic_ds2 = models.genetic_code_model(ds2, gc_df)

    print("[3/4] Computing metabolic cost model expectations…")

    cost_table = models.load_amino_acid_costs()
    expected_metabolic_ds1 = models.metabolic_model(ds1, cost_table)
    expected_metabolic_ds2 = models.metabolic_model(ds2, cost_table)

    correlations = models.compute_correlations(
        ds1, ds2,
        expected_genetic_ds1, expected_genetic_ds2,
        expected_metabolic_ds1, expected_metabolic_ds2,
    )

    io.save_table(correlations, cfg.processed_dir / "model_correlations.tsv")


# ---------------------------------------------------------------------
# STEP 4 — FIGURES
# ---------------------------------------------------------------------

def generate_figures(
    cfg: AminoAcidPipelineConfig,
    figures_dir: str | pathlib.Path
) -> None:

    figures_dir = pathlib.Path(figures_dir)
    figures_dir.mkdir(parents=True, exist_ok=True)

    print("[4/4] Generating paper-like figures…")

    ds1 = io.load_table(cfg.merged_aa_freq_unweighted())
    ds2 = io.load_table(cfg.merged_aa_freq_weighted())

    plots.scatter_observed_vs_genetic(ds1, ds2, out_dir=figures_dir)
    plots.scatter_observed_vs_metabolic(ds1, ds2, out_dir=figures_dir)
    plots.organism_level_correlations(ds1, ds2, out_dir=figures_dir)


# ---------------------------------------------------------------------
# FULL PIPELINE
# ---------------------------------------------------------------------

def run_full_pipeline():
    """
    Full workflow replicating the analysis pipeline behind:

        Krick, Verstraete et al., MBE 2014
        "Amino Acid Metabolism Conflicts with Protein Diversity"
    """

    cfg = AminoAcidPipelineConfig(
        raw_dir="data/raw",
        processed_dir="data/processed",
        organisms=[
            # Example, replace with your actual 17:
            "Yeast",
            "Ecoli",
            "Human",
        ],
    )

    prepare_proteomes(cfg)
    compute_amino_acid_frequencies(cfg)
    compute_models_and_correlations(cfg, gc_table_path="data/gc_content.tsv")
    generate_figures(cfg, figures_dir="results/figures")


if __name__ == "__main__":
    run_full_pipeline()
