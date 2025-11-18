# Deconvolution Scripts (R & Python)

## Repository Structure

```

gemdecan-deconvolution/
├── pipeline/        # R scripts (2020–2021) + sample parser (Python)
├── revision/        # Python scripts (2021–2022) for ICB prediction
├── signatures/      # RNA / Methylation / Hybrid signatures
└── results/         # Deconvolution outputs + QC + figures

````

---

## Methods

### **R Pipeline**
- `count_to_tpm.R` — convert counts → TPM  
- `compute_geneLength.R` — extract gene lengths (tximport)  
- `deconvolution_epidish.R` — EpiDISH / RPC  
- `deconvolution_deconrnaseq.R` — DeconRNASeq  
- `deconvolution_quantiseq.R` — quanTIseq  
- `deconvolution_mcpcounter.R` — MCP-counter  
- `merge_deconv.R` / `merge_quantif.R` — unify all outputs  
- `RNAsign_functions.R` — utility functions for signature handling  
- `deconvolution_algorithms.R` — wrapper combining all methods

### **Python (Revision 2021–2022)**
- `predict_from_deconv_revision*.py` — ElasticNet models, cross-cohort testing  
- `heatmap_ROC.py` — ROC comparison across signatures  
- `sample_parser.py` — sample and matrix formatting  
- SHAP, coefficients and predictions via the revision scripts

---

## Pipeline Overview

```mermaid
flowchart TD

    %% RAW → TPM
    A[Counts / Salmon quant] --> B[count_to_tpm.R]
    B --> C[TPM matrix]

    %% TPM → SIGNATURES
    C --> D1[BPRNA / BPRNACan]
    C --> D2[BPmet / BPmetCan]
    C --> D3[Hybrid RNA+Meth+3D]

    %% SIGNATURES → DECONV
    D1 --> E1[EpiDISH\n deconvolution_epidish.R ]
    D1 --> E2[DeconRNASeq\n deconvolution_deconrnaseq.R ]
    D1 --> E3[quanTIseq\n deconvolution_quantiseq.R ]
    D1 --> E4[MCP-counter\n deconvolution_mcpcounter.R ]

    D2 --> E1
    D2 --> E2
    D2 --> E3
    D2 --> E4

    D3 --> E1
    D3 --> E2

    %% MERGING
    E1 --> F[merge_deconv.R]
    E2 --> F
    E3 --> F
    E4 --> F

    F --> G[Unified cell-type proportions]

    %% REVISION (PYTHON)
    G --> H1[predict_from_deconv_revision.py\nElasticNet]
    G --> H2[heatmap_ROC.py]
    G --> H3[SHAP\n+ cohort transfer]

    H1 --> I[ICB Response Models]
    H2 --> I
    H3 --> I
````

---

## Dependencies

**R**: limma, EpiDISH, DeconRNASeq, MCPcounter, immunedeconv, xCell, tximport, tidyverse
**Python**: numpy, pandas, scikit-learn, matplotlib, seaborn, shap

---