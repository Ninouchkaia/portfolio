# **Drug Repurposing for COVID-19 through Network Medicine**

**Affiliations:** INSERM U1037 – Centre de Recherches en Cancérologie de Toulouse (CRCT)  
**Supervision:** Manlio De Domenico & Vera Pancaldi   
**Period:** 2020–2021  
**Publication:** [Network and Systems Medicine, 2020](https://www.liebertpub.com/doi/10.1089/nsm.2020.0011)  

<p align="center"><img width="1000" alt="multilayerbig" src="https://github.com/user-attachments/assets/a2ab551c-a40a-4419-af00-7c19d61c343a"></p>
<p align="left">
  <em>CovMulNet19 COVID-19 genotype–phenotype–drug interaction network. Result of the data integration and processing procedures. 
    (A) Nodes and schematic map of interdependencies among different layers encoding diseases, symptoms, drugs, GO terms, human proteins, and viral proteins. 
    (B) Map of the reconstructed structural interactions (e.g., protein–protein) and functional interdependencies (e.g., protein–disease, protein–GO term, or disease–symptom). Overall, the network consists of 1999 protein–protein, 19,755 protein–disease, 10,152 protein–symptom, 13,018 drug–target, 9210 protein–GO, and 3056 disease–symptom relationships.
  </em>
</p>

## Context
At the onset of the COVID-19 pandemic, identifying potential therapeutic candidates required integrative strategies beyond single-target screening. This project used **network medicine** approaches to explore interactions between SARS-CoV-2 proteins, host cellular pathways, and drug targets, with the goal of repositioning existing compounds.  


## Objectives
- Integrate multi-omics and molecular interaction data to construct a **virus–host–drug network**.  
- Identify biologically plausible drug candidates through **topological proximity** and **pathway enrichment**.  
- Test robustness of network-based predictions using simulated perturbations.  


## Methods
- **Data integration:** Host–virus interactome from public datasets (BioGRID, IntAct), drug–target relationships from DrugBank and ChEMBL.  
- **Network modeling:** Weighted graph representation of molecular associations.  
- **Simulation:** Random rewiring and node removal to assess prediction stability.  
- **Analysis:** Centrality and community detection to highlight key druggable modules.  
- **Validation:** Cross-checking candidate lists with published clinical data and ongoing trials.  

## Contributions
- Implemented random network simulations to evaluate robustness of predicted drug–disease associations.  
- Automated analysis of node connectivity and topological metrics for ranking candidate drugs.  
- Contributed to visualization and reporting of systemic network perturbations.  
- Participated in manuscript review and interpretation of results.  

## Reference
*Verstraete N.*, et al. *CovMulNet19, Integrating Proteins, Diseases, Drugs, and Symptoms: A Network Medicine Approach to COVID-19.*  
*Network and Systems Medicine*, 2020. [DOI:10.1089/nsm.2020.0011](https://www.liebertpub.com/doi/10.1089/nsm.2020.0011)

## Multilayer and Bootstrap Analysis Workflow
See [here](README_DETAILS.md) for the analysis workflow procedure.

