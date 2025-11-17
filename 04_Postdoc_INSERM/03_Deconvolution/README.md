# GEM-DeCan — Methodological Contributions

This section describes the methodological work I performed on the GEM-DeCan project, focusing on the reproducibility, validation, and biological coherence of the multi-omics pipeline used to generate and evaluate the GEM-DeCan signatures.

---

## 🔧 1. Pipeline Audit, Standardization and Documentation (Snakemake)

### Comprehensive pipeline audit
- Inspected all Snakemake rules (inputs, outputs, wildcards, DAG).
- Corrected missing or inconsistent dependencies.
- Standardized folder structure and naming conventions across all workflow steps.
- Identified and removed orphan or redundant rules.

### Reproducible environments
- Created and harmonized conda environments for each module (QC, alignment, quantification, deconvolution).
- Locked tool versions for full end-to-end reproducibility (FastQC, TrimGalore, STAR, Salmon, EpiDISH, deconRNAseq, quanTIseq, MCP-counter).
- Ensured cross-platform compatibility.

### Documentation
- Wrote detailed workflow documentation:
  - step-by-step description of each rule,
  - configuration parameters,
  - execution examples,
  - schema of the full pipeline from FASTQ to cell-type proportions.

---

## 🧮 2. Validation of Deconvolution Modules

### EpiDISH / RPC
- Verified correct loading of BPmet and BPmetCan methylation signatures.
- Checked CpG/group matching and format consistency.
- Validated estimates on PBMC mixtures (simulated and real).

### deconRNAseq
- Evaluated BPRNA, BPRNACan, and CCLE_TIL10 signatures.
- Verified mixture reconstruction stability and positivity constraints.
- Tested performance on in silico mixtures.

### MCP-counter & quanTIseq
- Ensured correct mapping between signature genes and sample matrices.
- Verified output structure, reproducibility, and score interpretation.

### Immunotherapy prediction module
- Checked consistency of generated matrices (proportions → model input).
- Verified train/test splits, seeds, and feature alignment.

---

## 🧬 3. Transcription Factor Activity Analysis

### Data preparation
- Generated expression matrices (TPM / logTPM) compatible with TF activity inference.
- Harmonized gene naming conventions and filtered relevant gene sets.

### TF activity inference
- Computed regulatory activity for TFs associated with:
  - M1 macrophages (IRF1, IRF5, STAT1, NFκB),
  - M2 macrophages (STAT6, PPARγ, C/EBPβ),
  - T/NK/B-cell programs (T-bet, BATF, Eomes, NFAT).

### Biological validation
- Confirmed biological coherence between inferred cell-type proportions and TF activation patterns.
- Demonstrated improved macrophage M1/M2 discrimination when using signatures augmented by methylation and Hi-C (BPRNACanProMet, BPRNACan3DProMet).
- Identified and documented discordances for quality control and further model refinement.

---

## 📊 4. End-to-End Quality Control

- Ran complete FASTQ → deconvolution tests on:
  - PBMC (public datasets),
  - TCGA bulk samples,
  - in silico mixtures.
- Verified consistency across deconvolution methods (EpiDISH, deconRNAseq, MCP-counter, quanTIseq).
- Checked purity estimates against reference methods (ABSOLUTE, LUMP, IHC, ESTIMATE).
- Ensured stability of signatures across datasets and species.

---

## ✔ Summary

My work ensured that GEM-DeCan is:
- fully reproducible (Snakemake, conda, version control),
- biologically coherent (TF activity validation),
- technically robust (QC on all modules and datasets),
- suitable for external users and reviewers.

```mermaid
flowchart TD

    %% --------------------
    %% RAW DATA
    %% --------------------
    A[FASTQ\nBulk RNA-seq] --> B[Quality Control\nFastQC / TrimGalore]
    B --> C[Quantification\nSTAR / Salmon TPM]

    %% --------------------
    %% SIGNATURES
    %% --------------------
    C --> D1[BPRNA\nBlueprint RNA signature]
    C --> D2[BPRNACan\nRNA + cancer-enriched genes]
    C --> D3[CCLE_TIL10\nTIL10 extended with CCLE/GTEx]

    C --> E1[BPmet\nBlueprint WGBS methylation]
    C --> E2[BPmetCan\nImmune + normal + cancer CpGs]

    %% Hybrid signatures
    E2 --> F1[BPRNACanProMet\nAdd genes with CpG diff.\nin promoters]
    E2 --> F2[BPRNACan3DProMet\nAdd genes with CpG diff.\nin promoters\n+ 3D-contact regions]

    %% --------------------
    %% DECONVOLUTION METHODS
    %% --------------------
    D1 --> G1[EpiDISH / RPC]
    D2 --> G1
    D3 --> G1

    E1 --> G2[deconRNAseq]
    E2 --> G2

    F1 --> G1
    F2 --> G1

    C --> G3[MCP-counter\n non-proportional counts]
    C --> G4[quanTIseq]

    %% --------------------
    %% OUTPUT: PROPORTIONS
    %% --------------------
    G1 --> H[Cell-type Proportions]
    G2 --> H
    G3 --> H
    G4 --> H

    %% --------------------
    %% QC / VALIDATION
    %% --------------------
    H --> I1[Validation\nPBMC real + simulated]
    H --> I2[TCGA bulk\nTumor Purity vs ABSOLUTE/LUMP/IHC]
    H --> I3[H&E image-based cell density\n Saltz et al.]
    H --> I4[Single-cell reconstructed bulk\n Melanoma - Tirosh et al.]

    %% --------------------
    %% IMMUNOTHERAPY MODELS
    %% --------------------
    H --> J1[ElasticNet models\nHugo dataset]
    H --> J2[ElasticNet models\nGide dataset]
    H --> J3[ElasticNet models\nRiaz dataset]

    J1 --> K[Prediction of response\nto anti-PD1 therapy]
    J2 --> K
    J3 --> K
```