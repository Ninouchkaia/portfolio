# Evolutionary Constraints on Amino Acid Usage in Proteomes

**Affiliation:** INQUIMAE – CONICET, University of Buenos Aires  
**Period:** 2013–2015  
**Publication:** [Molecular Biology and Evolution, 2014](https://doi.org/10.1093/molbev/msu228)  

---

## 🧭 Context
Protein composition reflects a balance between evolutionary pressures, biochemical stability, and energetic cost.  
This project aimed to understand how **metabolic constraints** influence amino acid usage across proteomes, and how this trade-off affects protein diversity across species.

---

## 🎯 Objectives
- Quantify amino acid usage across species in relation to biochemical properties and synthesis cost.  
- Evaluate evolutionary trade-offs between energetic efficiency, chemical stability, and proteome diversity.  
- Develop a quantitative framework describing these multi-objective constraints.  

---

## 🧪 Methods
- **Data acquisition:** Protein abundance data retrieved from **PaxDB** for 17 model organisms.  
- **Processing:** Cleaning, structuring, and cross-validation of abundance-weighted proteome datasets.  
- **Analysis:** Calculation of amino acid usage frequencies weighted by abundance.  
- **Comparative modeling:** Statistical evaluation of correlations between amino acid cost, stability, and frequency.  
- **Visualization:** Multidimensional scaling and regression analyses representing trade-offs between energy and diversity.  

---

## 💡 Contributions
- Retrieved and processed proteome-wide abundance data from PaxDB.  
- Quantified amino acid usage patterns and computed cost–diversity correlations.  
- Contributed to statistical modeling and graphical representation of trade-offs.  
- Participated in manuscript preparation and data validation.  

---

## 📘 Key Skills
Comparative genomics · Evolutionary modeling · Proteomics data analysis · Statistical inference · Python / R scripting  

---

## 🔗 Reference
*Krick T., Verstraete N.*, Alonso L.G., Shub D.A., Ferreiro D.U., Sanchez I.E.  
*Amino Acid Metabolism Conflicts with Protein Diversity.*  
*Molecular Biology and Evolution*, 2014. [DOI:10.1093/molbev/msu228](https://doi.org/10.1093/molbev/msu228)

## Figures

Amino acid relative abundances in proteomes

We estimate amino acid relative abundances in proteomes in two datasets. Dataset DS1 was derived from108 fully sequenced and annotated genomes from the three domains of life (Tekaia and Yeramian, 2006).We translated coding regions into protein sequences and counted the frequency of occurrence of each aminoacid, assuming that all proteins are equally abundant (Table E1). Dataset DS2 was derived from the PaxDBdatabase for protein abundances (Wang et al , 2012). We considered 17 organisms for which protein sequenceand relative abundance data are available for more than 50 per cent of the proteome. We used integrateddatasets for the whole organism whenever possible (Table E2).For both datasets, we tested several models for amino acid relative abundances. The results are shownin Table I, Figure 1 and Figure E1 below and explained in the next sections.

Correlation of amino acid relative abundances with the genetic code modelThe genetic code model relates amino acid relative abundance with the transcription and translation ofrandom DNA sequences of a given GC content (Dyer, 1971; Gupta, 2005). To evaluate this model with DS1and DS2 we retrieved the genomic GC content for each genome from (Kryukov et al , 2012) and used it tocalculate the expected relative abundances for all 61 amino acid coding triplets. We then translated thetriplets into amino acids and obtained the expected amino acid relative abundances in each proteome. Thismetabolism-agnostic model shows a good correlation between calculated and observed amino acid relativeabundances (Table I and Figure 1, panels C and F). The rvalues are 0.71 and 0.62 for DS1 and DS2.The correlation is also observed for individual organisms in the database regardless of genomic GC content(Figure 2, dashed lines in panels A and B). However, the rvalues are worse than for our metabolism-basedmodel when amino acid costs are measured in units of ATP/time (Table I). This holds regardless of genomicGC content (Figure 2). The rvalue is better for our model in 105 of the 108 organisms in DS1 (Figure 2,Panel A) and for the 17 organisms in DS2 (Figure 2, Panel B). This conclusion is also valid if the amino acidcysteine is excluded from the calculations (Table I). We interpret that amino acid relative abundances arebetter explained when we take into account the simultaneous maximization of entropy and minimization ofcost


Additional Expanded View Figure Legends

Table E1 Dataset DS1 was derived from 108 fully sequenced and annotated genomes from the threedomains of life (Tekaia and Yeramian, 2006). Coding regions were translated into protein sequences and the frequency of occurrence of each amino acid was calculated, assuming that all proteins are equally abundant.The table shows values of amino acid relative abundances, predicted abundances from the genetic codemodel, genomic GC content and correlation R-values.

Table E2 Dataset DS2 was derived from 17 organisms from the PaxDB database for protein abundances(Wang et al, 2012). We considered organisms for which protein sequence and relative abundance data areavailable for more than 50 per cent of the proteome and used integrated datasets for the whole organismwhenever possible. The table shows values of amino acid relative abundances, predicted abundances fromthe genetic code model, genomic GC content and correlation R-values
