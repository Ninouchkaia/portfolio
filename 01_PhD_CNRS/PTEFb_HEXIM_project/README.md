# Structure–Function Analysis of P-TEFb Transcription Regulation by HEXIM1 and HIV-1 Tat

**Affiliation:** CNRS / École Normale Supérieure – Institut de Biologie de l’ENS (IBENS), Département de Génomique Fonctionnelle  
**Period:** 2008–2012  
**Supervision:** Philippe Bensaude, CNRS UMR 8197  
**Thesis Title:** *Régulation du facteur d’élongation transcriptionnelle P-TEFb par les protéines HEXIM1 et TAT du VIH-1*  

---

## 🧭 Context
The P-TEFb complex (Cyclin T1/CDK9) is a central regulator of RNA polymerase II transcriptional elongation.  
Its activity is controlled by the **7SK RNP complex**, where **HEXIM1** acts as an inhibitory subunit, and hijacked by the **HIV-1 Tat** protein to promote viral gene expression.  

This thesis investigated the structural and functional interfaces involved in these interactions to elucidate the molecular mechanisms of transcriptional control.

---

## 🎯 Objectives
- Map the protein–protein interaction surfaces within the P-TEFb complex.  
- Characterize how HEXIM1 and HIV-1 Tat compete for Cyclin T1 binding.  
- Explore the impact of these interactions on transcriptional regulation in mammalian and model systems.  

This research project combined:

* **High-throughput yeast genetics (reverse two-hybrid)**
* **Mutagenesis of Cyclin T1 and HEXIM1**
* **Cell-based assays in human and mouse cells**
* **Fluorescence microscopy (localization, FLIP, single-molecule tracking)**
* **Biochemical assays (co-IP, CTD kinase assays)**
* **RNA level analyses (7SK detection, phylogeny, structural modeling)**
* **Comparative genomics and evolutionary reconstruction of 7SK RNA and its protein partners**

---

## 💡 Contributions
- Constructed mutant libraries of Cyclin T1 and HEXIM1 to identify contact regions with P-TEFb and Tat.  
- Mapped functional domains responsible for transcriptional activation or inhibition.  
- Combined biochemical and imaging approaches to establish a structure–function model of P-TEFb regulation.  
- Contributed to internal publications and collaborative manuscripts within the Bensaude group.

---

Voici une **synthèse unique**, **non redondante**, **compacte mais exhaustive**, qui fusionne *toutes* les sections que tu m’as données.
Elle garde **toutes les informations importantes**, sans répétition, et en un bloc parfaitement clair pour un README, une thèse, un portfolio ou un dossier de candidature.

---

## **Experimental Methods**

### **Yeast Genetics & Mutagenesis**

* Random mutagenesis of Cyclin T1 by error-prone PCR.
* Targeted mutagenesis of HEXIM1 (alanine scanning; conserved motifs).
* In vivo homologous recombination in yeast to rebuild full-length constructs.
* Two-hybrid and reverse two-hybrid assays to screen interaction-deficient mutants.
* 5-FOA counter-selection to isolate Cyclin T1 variants losing HEXIM1 binding.
* Validation of mutant phenotypes using HEXIM1 and CDK9 interaction assays.

---

### **Molecular Cloning & Construct Engineering**

* Design and assembly of mutagenized Cyclin T1 and HEXIM1 alleles.
* Generation of Dendra-tagged Cyclin T1 constructs for live imaging.
* Reconstruction and mapping of interaction surfaces through combinatorial cloning.

---

### **Mammalian Cell Culture & Functional Assays**

* Transient expression in HEK293T and NIH-3T3 cells.
* Establishment of stable cell lines expressing Dendra-CycT1.
* Functional readouts through HIV LTR reporter activation (Gal4-CycT1).
* Tat rescue assays to assess the ability of mutant Cyclin T1 to support viral transcription.

---

### **Protein Interaction & Biochemistry**

* Immunoprecipitation and co-IP of Cyclin T1, CDK9, and HEXIM1 complexes.
* Western blot analysis and quantification of protein–protein interactions.
* CTD kinase assays to measure P-TEFb catalytic activity.
* RNAse-dependent assays to evaluate the contribution of 7SK-HEXIM regulation.

---

### **Fluorescence Imaging & Live-Cell Dynamics**

* Confocal microscopy to assess localization of Cyclin T1 mutants.
* FLIP measurements to quantify nuclear mobility and binding engagement.
* Single-molecule tracking in Dendra-CycT1 stable cell lines.
* Analysis of C-terminal truncations and mobility phenotypes.

---

### **RNA & Genomics Analyses**

* Identification of functional 7SK RNA homologs in *C. elegans*.
* Evolutionary inference of the ancient metazoan origin of the 7SK–HEXIM regulatory system.

---

## **Key Results **

* Identification of a **Cyclin T1 interaction groove** required for HEXIM1 binding, centered on **hotspot Y175**.
* Demonstration that **Y175 mutations disrupt HEXIM1 and HIV Tat binding**, abolishing TAR association and LTR activation.
* Mapping of **critical HEXIM1 residues** (F262, F267, H275) located in an **intrinsically disordered regulatory region**.
* Discovery of **Cyclin T1 mutants** that lose HEXIM1 binding but retain CDK9 binding, producing **hyperactive P-TEFb** insensitive to 7SK repression.
* Definition of how **C-terminal regions of Cyclin T1** influence RNAPII interaction and **nuclear mobility** (FLIP, single-molecule imaging).
* Evolutionary reconstruction showing that the **7SK–HEXIM regulatory system predates vertebrates**, including discovery of a **functional nematode 7SK RNA** and conserved structural modules (M1–M8).

---

## **Publications**

* **N. Verstraete et al. (2014)**
  *A Cyclin T1 point mutation that abolishes P-TEFb binding to HEXIM1 and HIV Tat.* Retrovirology.

* **M. Marz, A. Donath, N. Verstraete et al. (2009)**
  *Evolution of 7SK RNA and its protein partners in Metazoa.* Mol Biol Evol.

---

## **7. Repository Structure**

```
PTEFb_HEXIM_project/
│
├── README.md
├── docs/
│   ├── figures/
│   │   ├── pipeline_mermaid.svg
│   │   ├── interaction_map.svg
│   │   └── 7SK_secondary_structure.png
│   ├── thesis_summary.pdf
│   ├── retrovirology_2014_summary.pdf
│   └── mbe_2009_summary.pdf
│
├── data/
│   ├── mutants_cyclinT1/
│   ├── yeast_screens/
│   ├── microscopy/
│   └── sequences_7SK_HEXIM/
│
├── scripts/
│   ├── analysis/
│   │   ├── rnabob_automaton.py
│   │   ├── sequence_alignment.sh
│   │   └── secondary_structure_viennaRNA.R
│   └── plotting/
│       └── visualize_mutations.ipynb
│
└── notebooks/
    ├── 01_mutagenesis_screen.ipynb
    ├── 02_coIP_and_kinase.ipynb
    ├── 03_microscopy_analysis.ipynb
    ├── 04_7SK_evolution.ipynb
    └── 05_summary_figures.ipynb
```







