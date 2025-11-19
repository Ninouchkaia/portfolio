# Ankyrin Modularity Pipeline

**Domain enrichment • Conservation • ELM/SLiM analysis • Interacting partners**

This repository contains a modular and reproducible Python implementation of the analysis pipeline used to study:

> **Functional modules built by association of domains and linear motifs in Ankyrin proteins and their binding partners**

* PFAM parsing (query / homologs)
* Pfam → clan mapping
* Domain conservation across homologs
* Log(obs/exp) domain enrichment + conservation-colored scores
* Automated extraction of domain FASTAs and independent regions
* Full ELM/SLiM processing (counts, enrichment, domain–SLiM co-occurrence)
* Interaction-level domain/ELM association analysis

Each analysis stage is activated through **explicit flags in `pipeline.py`**.

---

# Project Structure

```
ankyrin_modularity/
│
├── config.py
├── pipeline.py
│
├── io.py
├── clans.py
├── conservation.py
├── enrichment.py
├── elm.py
│
└── data/
    ├── raw/
    │   ├── Ank_uniref50_string_mapped_522.fasta
    │   ├── binding_partners_2038.fasta
    │   ├── Pfam-A.clans.txt
    │   ├── FrequenciesPfam_mapped.txt
    │   ├── Pfam_domains_in_Ank1234_Zscores.txt
    │   └── Pfam_domains_in_BD_2038_Zscores.txt
    │
    ├── intermediate/
    │   ├── pfam_queries/
    │   ├── pfam_homologs/
    │   └── conservation / concatenated PFAM / intermediate tables
    │
    └── results/
        ├── sequences/
        │   ├── domains/
        │   └── independent/
        └── enrichment / ELM tables / co-occurrences / plots
```

---

# Running the Pipeline

All steps are controlled from the top of `pipeline.py`:

```python
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
```

Activate the steps you want by switching flags to **True**, then run:

```
python -m ankyrin_modularity.pipeline
```

Logs are automatically written to:

```
ankyrin_modularity/ankyrin_modularity.log
```

---

# Pipeline Diagram

### **Simple ASCII overview**

```
          +---------------+
          |  config.py    |
          +-------+-------+
                  |
     +------------+--------------+
     |                           |
 +---v---+                  +----v-----+
 | io.py |                  | clans.py |
 +---+---+                  +----+-----+
     |                           |
 +---v---------------+   +-------v-----------+
 | conservation.py   |   | enrichment.py     |
 +--------+----------+   +----------+--------+
          \                     /
           \                   /
            +--------v--------+
            |     elm.py     |
            +--------+--------+
                     |
                 +---v----+
                 |pipeline|
                 +--------+
```

---

### **Detailed pipeline (GitHub-friendly ASCII)**

```
 ┌──────────────────────── Ankyrin Modularity Pipeline ─────────────────────────┐

   Raw FASTA      PFAM scan outputs        Pfam-A.clans.txt
        │                │                         │
        ▼                ▼                         ▼
┌──────────────┐   ┌──────────────┐       ┌────────────────┐
│    io.py     │→→ │ conservation │       │    clans.py     │
│ FASTA index  │    domain %      │ → →   │ PFAM→clan map   │
│ PFAM parser  │    cons tables   │       │ clan freq table │
└──────┬───────┘   └──────┬───────┘       └────────┬────────┘
       │                  │                         │
       │                  ▼                         │
       │          ┌──────────────┐                  │
       │          │ enrichment.py│←─────────────────┘
       │          │ log(obs/exp) │
       │          │ color by cons│
       │          └───────┬──────┘
       │                  │
       ▼                  ▼
 ┌──────────────┐   ┌────────────────┐
 │ domain FASTA │   │  ELM tables     │
 │ extraction   │   │  Ank / Partners │
 └──────────────┘   └──────┬─────────┘
                            ▼
                      ┌──────────────┐
                      │    elm.py     │
                      │ count ELMs    │
                      │ ELM enrich.   │
                      │ ELM–domain    │
                      └──────┬────────┘
                             ▼
                        ┌────────────┐
                        │ pipeline.py│
                        │ orchestration
                        └────────────┘

 └──────────────────────────────────────────────────────────────────────────────┘
```

---

# Module Summary

### **io.py**

FASTA indexing, PFAM parsing, automatic PFAM file detection, domain FASTA extraction, independent region extraction, PFAM concatenation.

### **clans.py**

Pfam → clan mapping, aggregated clan frequencies.

### **conservation.py**

Per-domain conservation over homologs, conservation state tables, percentages.

### **enrichment.py**

Domain log(obs/exp) enrichment, combination with conservation + color assignment.

### **elm.py**

SLiM/ELM reading, counting, enrichment, domain–ELM co-occurrence (interacting pairs).

---

The pipeline avoids overwriting results unless explicitly intended and organizes them under:

```
data/intermediate/
data/results/
data/results/sequences/
```

---


