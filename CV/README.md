# **Portfolio - Research, Bioinformatics & Computational Biology**

## Research Focus

- **Bioinformatics & Data Analysis**
  - RNA-seq (bulk & single-cell), deconvolution, and gene regulatory inference  
  - Network medicine and multi-omics data integration  
  - Statistical modeling and visualization of biological data  

- **Modeling & Systems Biology**
  - Agent-based modeling of tumor-immune ecosystems  
  - Multi-scale simulations of cellular dynamics 
  - Network perturbation and system-level functional analysis  

- **Molecular & Structural Biology**
  - Protein evolution and structural modularity  
  - Regulatory mechanisms in transcription and post-transcriptional control  

- **Scientific Communication & Pedagogy**
  - Teaching, mentoring, and outreach in coding
  - Integration of artistic coding practices for public engagement  

---

## List of Research Projects

### [**INSERM/CRCT - Postdoc (2020-2023)**](04_Postdoc_INSERM/)

#### [1. Tumor Ecosystem Modeling](04_Postdoc_INSERM/01_AgentBasedModel/)

This project models the spatial and temporal dynamics of CLL tumor cells and nurse-like cells to understand how microenvironmental interactions shape patient-specific growth trajectories. The goal was to calibrate an agent-based model against experimental co-culture data using multi-objective optimization. The study identifies parameter regimes that reproduce observed CLL-NLC coexistence patterns.

#### [2. Macrophage Polarization in CLL](04_Postdoc_INSERM/02_BooleanModel)

This project investigates how transcription factors drive monocyte-to-macrophage polarization in the CLL microenvironment. A Boolean regulatory model was used to map distinct polarization trajectories and link them to experimental signatures. The analysis identifies regulatory determinants underlying M1/M2/NLC-like phenotypes.

#### [3. Predicting Response to Immunotherapy (GEMDECAN)](04_Postdoc_INSERM/03_RNASeqDeconvolution/)

The project explores how immune composition and gene expression programs predict checkpoint inhibitor response in solid tumors. A multi-cohort RNA-seq workflow quantifies immune infiltration and computes signature-based predictors of clinical outcome. The analysis aims to identify biomarkers associated with therapeutic benefit.

#### [4. LungPredict - Transcription Factor Network Analysis](04_Postdoc_INSERM/04_LungPredict)

This project examines how transcription factor activity patterns stratify lung cancer patients into biologically meaningful subgroups. Multi-omic annotations and PCA reveal TF-driven regulatory programs associated with clinical phenotypes. The analysis highlights regulatory axes shaping tumor heterogeneity.

#### [5. COVID-19 Drug Repurposing via Network Medicine](04_Postdoc_INSERM/05_NetworkMedicine/)

This work models how pharmacological perturbations propagate through SARS-CoV-2 host interaction layers. By comparing observed multilayer network structure to bootstrapped null models, it identifies drugs targeting essential viral-host modules. The objective is to highlight compounds with mechanistically grounded repurposing potential.

#### [6. Systemic Effects of SARS-CoV-2](04_Postdoc_INSERM/06_Sarscov2Interactome)

This project characterizes the systemic impact of SARS-CoV-2 by integrating viral-host interaction datasets and pathway annotations. Enrichment analysis highlights the cellular pathways and tissues most affected by viral proteins. The goal is to map functional vulnerabilities and host responses.

#### [7. Clonal Dynamics and scRNA-seq under Treatment](04_Postdoc_INSERM/07_BarcodesDrugScreening/)

This project investigates how clonal populations respond to targeted and chemotherapeutic drugs by quantifying barcode representation under treatment. Differential abundance analysis reveals drug-specific clonal signatures and functional clusters. Correlation networks map convergent and divergent therapeutic effects.

---

### [**Software Engineering for Airbus Systems (Toulouse, 2017-2020)**](03_Industry_AIRBUS/)
*Areas of interest: Development, integration and maintenance of PLM/PDM systems (APS, CASPARE, PASS V3 migration, Windchill v6→v11).*

---

### [**CONICET/UBA - Postdoc (2013-2015)**](02_Postdoc_CONICET)
*Areas of interest: Structural bioinformatics, protein evolution, and interaction modularity*

#### [1. Amino Acid Usage under Evolutionary Constraints](02_Postdoc_CONICET/paxdb/)

This work explores how metabolic cost and biosynthetic limitations influence amino acid composition across proteomes. By combining abundance-weighted PaxDB datasets with biochemical cost models, it tests evolutionary trade-offs between energy efficiency and proteome diversity. The results provide a quantitative framework linking metabolism to sequence composition.

#### [2. Structure and Dynamics of Ankyrin Repeats](02_Postdoc_CONICET/AnkyrinStructure/)

This project examines how structural variability and conserved features define the stability of ankyrin repeat proteins. A curated structural dataset enables analysis of contacts, repeat architecture, and positional variability. The study provides insight into how repeat proteins achieve robustness despite modular design.

#### [3. Functional Modularity of Ankyrin Proteins and Partners](02_Postdoc_CONICET/ANKYRIN_MODULARITY/)

The goal of this project is to understand how ankyrin domains combine with additional motifs to encode specific interaction functions. Conservation, Pfam/ELM enrichment and co-occurrence analyses reveal modular rules governing partner specificity. The work outlines how repeat and non-repeat elements assemble into functional architectures.

---

### [**CNRS/IBENS - PhD (2008-2012)**](01_PhD_CNRS)
*Areas of interest: Regulation of transcriptional elongation and structure-function analysis*

#### [1. P-TEFb Regulation by HEXIM1 and HIV-1 Tat](01_PhD_CNRS)

This project identifies the Cyclin T1 structural determinants underlying binding to HEXIM1 and HIV-1 Tat. Mutagenesis and biochemical assays map critical hotspots involved in transcriptional elongation regulation. The findings clarify how viral proteins hijack the host machinery.

#### [2. P-TEFb Mobility in Living Cells](01_PhD_CNRS)
This study examines how complex formation modulates P-TEFb nuclear mobility using FRAP/FLIP imaging. Truncation mutants and chemical perturbations reveal diffusion-binding equilibria within the nucleus. The work connects molecular assembly to dynamic regulation of transcription.

#### [3. P-TEFb Conservation Across Metazoa](01_PhD_CNRS)

This project investigates whether CyclinT-HEXIM interactions are conserved across metazoans. Comparative analysis and *C. elegans* validation show that critical interface residues remain functionally preserved. The study highlights ancestral mechanisms controlling transcriptional elongation.

---

## Research Projects Summary

| Project | Summary | Biological Question | Contribution |
|--------|---------|----------------------|---------------|
| [Tumor Ecosystem Modeling](04_Postdoc_INSERM/01_AgentBasedModel/) | Models how CLL-NLC interactions generate patient-specific tumor dynamics through calibrated agent-based simulations. | How do CLL-NLC interactions reproduce patient-specific tumor dynamics? | Full ABM pipeline (NetLogo → OpenMOLE → Python), NSGA-II calibration, patient fitting, validation, advanced stats |
| [Macrophage Polarization in CLL](04_Postdoc_INSERM/02_BooleanModel) | Explores transcription factor programs driving monocyte-to-macrophage polarization in the CLL microenvironment. | Which TFs drive monocyte → macrophage polarization in CLL? | TF activity scoring (VIPER/DoRothEA), scaling/normalisation |
| [Predicting Response to Immunotherapy (GEMDECAN)](04_Postdoc_INSERM/03_RNASeqDeconvolution/) | Identifies immune and transcriptomic predictors of checkpoint inhibitor response across patient cohorts. | Can RNA-seq deconvolution predict responders? | Snakemake workflow: TPM, immune deconvolution, modelling |
| [LungPredict - Transcription Factor Network Analysis](04_Postdoc_INSERM/04_LungPredict) | Reveals transcription factor activity patterns that stratify lung cancer patients into regulatory subgroups. | How do TF programs stratify lung cancer patients? | TF activity inference, heatmaps, multi-omic annotation, PCA |
| [COVID-19 Drug Repurposing via Network Medicine](04_Postdoc_INSERM/05_NetworkMedicine/) | Evaluates how drugs perturb multilayer SARS-CoV-2 host networks to highlight mechanistically relevant candidates. | Which drugs most strongly perturb human multilayer networks implicated in SARS-CoV-2 infection? | Multilayer networks, bootstrap null models, z-score ranking |
| [Systemic Effects of SARS-CoV-2](04_Postdoc_INSERM/06_Sarscov2Interactome) | Maps viral-host protein interactions to identify the pathways and tissues most affected by SARS-CoV-2 infection. | Which pathways and tissues are most impacted? | Viral-host merge, Reactome enrichment, annotation, figures |
| [Clonal Dynamics and scRNA-seq under Treatment](04_Postdoc_INSERM/07_BarcodesDrugScreening/) | Quantifies clonal responses to drug treatments through barcode abundance shifts and correlation-based drug clustering. | How do clonal populations respond to drugs? | QC, DESeq2 input generation, logFC matrices, correlation networks |
| [Amino Acid Usage under Evolutionary Constraints](02_Postdoc_CONICET/paxdb/) | Links proteome amino acid frequencies to metabolic costs to reveal evolutionary constraints on sequence composition. | How do metabolic costs shape amino acid usage? | PaxDB parsing, weighted AA frequencies, cost modelling |
| [Structure and Dynamics of Ankyrin Repeats](02_Postdoc_CONICET/AnkyrinStructure/) | Characterizes structural variability and conserved features shaping the stability of ankyrin repeat proteins. | What features shape ankyrin repeat stability? | Dataset creation, repeat detection, structural analysis |
| [Functional Modularity of Ankyrin Proteins and Partners](02_Postdoc_CONICET/ANKYRIN_MODULARITY/) | Examines how ankyrin domains combine with Pfam/ELM motifs to encode modular interaction specificity. | How do ankyrin modules encode function? | Conservation, Pfam/ELM enrichment, co-occurrence networks |
| [P-TEFb Regulation by HEXIM1 and HIV-1 Tat](01_PhD_CNRS) | Maps CyclinT1 residues controlling HEXIM1 and Tat binding during transcriptional elongation regulation. | Which residues control HEXIM1/Tat interactions? | Mutagenesis, interface mapping, assays |
| [P-TEFb Mobility in Living Cells](01_PhD_CNRS) | Uses FRAP/FLIP to determine how complex assembly modulates P-TEFb nuclear diffusion dynamics. | How does complex assembly affect mobility? | FRAP/FLIP quantification, diffusion modelling |
| [P-TEFb Conservation Across Metazoa](01_PhD_CNRS) | Shows that CyclinT-HEXIM regulatory interactions are conserved across metazoans, including *C. elegans*. | Is CyclinT-HEXIM interaction conserved? | Conservation analysis + *C. elegans* validation |

---

## **Contact**

For collaborations or technical discussions, contact me at verstraete.nina[at]gmail.com.

---











































