# **P-TEFb / HEXIM / 7SK — Interaction Mapping & Evolutionary Analysis**

### **Doctoral Project – UPMC / ENS Paris**

---

## **1. Project Overview**

This repository summarizes a multi-year research project focused on the **structural and functional dissection of P-TEFb regulation**, combining:

* **High-throughput yeast genetics (reverse two-hybrid)**
* **Mutagenesis of Cyclin T1 and HEXIM1**
* **Cell-based assays in human and mouse cells**
* **Fluorescence microscopy (localization, FLIP, single-molecule tracking)**
* **Biochemical assays (co-IP, CTD kinase assays)**
* **RNA level analyses (7SK detection, phylogeny, structural modeling)**
* **Comparative genomics and evolutionary reconstruction of 7SK RNA and its protein partners**

The core objective was to **map the interaction surfaces between Cyclin T1 (CycT1), CDK9, HEXIM1, LARP7**, and the **7SK snRNA**, and to understand how this regulatory module evolved across Metazoa.

---

## **2. Biological Background**

### **2.1 P-TEFb (CDK9/Cyclin T)**

P-TEFb is the central transcription elongation factor that phosphorylates:

* **NELF, DSIF** (pausing release)
* **RNAPII CTD** (productive elongation)
* **RNA processing factors**

### **2.2 Inhibition by HEXIM1/7SK snRNA**

In cells, most P-TEFb is kept inactive by:

* **HEXIM1** bound to
* **7SK snRNA** assembled with **LARP7** and **MePCE**.

HEXIM1 and HIV-1 Tat compete for a **partially overlapping binding surface on Cyclin T1**.

### **2.3 Viral hijacking**

HIV-1 Tat binds the same groove of Cyclin T1 and recruits P-TEFb onto TAR RNA, strongly activating transcription.

---

## **3. Experimental Strategy**

### **3.1 Large-scale random mutagenesis of Cyclin T1**

* **Error-prone PCR** of defined CycT1 subdomains
* **Recombination in yeast** to rebuild full-length BD-CycT1
* **Reverse two-hybrid selection** on 5-FOA
* Identification of **single amino-acid substitutions disrupting HEXIM1 interaction**

### **3.2 Targeted mutagenesis of HEXIM1**

* Alanine scanning of conserved motifs (basic region, PYNT motif, dimerization helix)
* Identification of **key aromatic residues** essential for Cyclin-binding (F262, F267, H275)

### **3.3 Two-hybrid validation pipeline**

* HEXIM1 interaction
* CDK9 binding control
* Mutant selection based on **Hexim(–) / Cdk9(++)** profiles

### **3.4 Cell-based validation**

* **Co-immunoprecipitation** (CycT1 mutants + endogenous HEXIM1 or CDK9)
* **CTD kinase assays** with/without RNAse
* **Reporter assays** (Gal4-CycT1 to recruit P-TEFb)
* **HIV Tat rescue assays in mouse 3T3 cells**

### **3.5 Microscopy**

* Subcellular localization of Cyclin T1 mutants
* FLIP experiments to quantify mobility
* Generation of **stable Dendra-CycT1 cell lines** for single-molecule tracking

### **3.6 Comparative genomics**

* Exhaustive search of **7SK RNA**, **HEXIM**, **LARP7**, **MePCE**
* Discovery of **functional 7SK RNA in nematodes** (C. elegans)
* Reconstruction of **secondary structure variants (M1–M8 helices)**
* Evidence for **ancient metazoan origin** of the 7SK-HEXIM system

---

## **4. Key Results**

### **4.1 Mapping the P-TEFb–HEXIM1 binding interface**

* Mutations disrupting HEXIM1 binding cluster along a **groove between the two cyclin folds**.
* The **central hotspot is residue Y175** of Cyclin T1.
* Y175 establishes a **π-hydrogen bond network** bridging helices H1–H2/H1’–H2’.

### **4.2 Overlap with HIV Tat binding**

* Y175 mutants **cannot bind Tat**
* → cannot support Tat-mediated HIV LTR activation
* → cannot bind TAR RNA
* → do not recruit P-TEFb to HIV promoter

### **4.3 HEXIM1 structural determinants**

* Three conserved aromatic residues (F262, F267, H275) are **essential for Cyclin T1 docking**.
* They occur in an **intrinsically disordered region**, similar to regions in other CDK regulators (p21, p27).

### **4.4 Cellular phenotypes**

* CycT1-Y175 mutants show:

  * Loss of HEXIM1 binding
  * Normal CDK9 binding
  * Increased CTD-kinase activity (insensitive to 7SK/HEXIM repression)
  * Enhanced transcription activation when tethered to DNA

### **4.5 Single-molecule mobility**

* C-terminal truncations of CycT affect:

  * RNAPII association
  * Nuclear mobility (measured by FLIP and photoactivation)

### **4.6 Evolutionary findings (MBE 2009)**

* 7SK/HEXIM system present across most Metazoa.
* **Two alternative structural conformations** of 7SK M2 region → potential RNA switch.
* Discovery and experimental validation of **functional 7SK homolog in C. elegans**.
* HEXIM and LARP7 underwent **lineage-specific duplications and losses** (HEXIM1/2 in mammals).

---

## **5. ExperimentalMethods Summary**

* **Yeast genetics**: standard 2-hybrid / reverse 2-hybrid
* **Molecular cloning**: mutagenesis, recombination
* **Cell culture**: HEK293T, 3T3, stable Dendra lines
* **Protein assays**: IP, co-IP, western blot
* **Kinase assays**: CTD4 phosphorylation
* **Microscopy**: confocal, FRAP/FLIP, single-molecule tracking

---

## **6. Project Pipeline**

```mermaid
flowchart TD

A[Random mutagenesis of Cyclin T1] --> B[Reverse two-hybrid selection]
B --> C[Sequencing & identification of point mutations]
C --> D[Two-hybrid validation: HEXIM(-), Cdk9(+)]
D --> E[Cell-based validation]
E --> E1[Co-IP: Cdk9 / HEXIM binding]
E --> E2[CTD kinase assays]
E --> E3[Reporter assays]
E --> E4[Tat rescue assays]

F[HEXIM1 targeted mutagenesis] --> G[Two-hybrid mapping]
G --> H[Definition of HEXIM interaction domain]

I[Microscopy experiments] --> J[Localization, FLIP, mobility]

K[Comparative genomics pipeline] --> L[7SK/HEXIM/LARP7 discovery]
L --> M[Secondary structure modeling]
M --> N[Nematode 7SK validation]

E & H & J & N --> O[Integrated model of P-TEFb regulation]
```

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

---

## **8. Publications**

* **N. Verstraete et al. (2014)**
  *A Cyclin T1 point mutation that abolishes P-TEFb binding to HEXIM1 and HIV Tat.* Retrovirology.

* **M. Marz, A. Donath, N. Verstraete et al. (2009)**
  *Evolution of 7SK RNA and its protein partners in Metazoa.* Mol Biol Evol.

---

## **9. Contact**

For questions or collaboration inquiries, please contact: **xxxxxxxxx**

---

