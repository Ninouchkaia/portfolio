# %% [markdown]
# # Barcode signatures and drug relationships
#
# This notebook visualizes the clonal barcoding drug screen:
#
# - Loads log2 fold-change (log2FC) per barcode per condition (DESeq2 output)
# - Explores distributions and basic QC
# - Recreates:
#   - Fig.3: drug clustering based on barcode signatures (log2FC heatmap)
#   - Fig.4: drug–drug correlation clustered heatmap
#   - Fig.5: drug network based on correlation
#
# ## File assumptions
#
# This notebook assumes the following project layout:
#
# ```text
# 2022_Barcodes/
# ├── data/
# │   ├── merged_pval_filtered_deseq2_fillna.csv
# │   ├── correl_matrix1_merged_pval_filtered_deseq2.csv      (optional, can be recomputed)
# │   ├── color_mapping.tsv                                   (optional)
# │   └── ...
# ├── results/
# │   └── paper_figures/
# └── notebooks/
#     └── 01_visualize_barcode_signatures.ipynb   (this notebook)
# ```
#
# If your paths differ, update the configuration cell below.

# %%
# Imports and configuration

from pathlib import Path
import os

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import networkx as nx

sns.set(style="white", context="notebook")

ROOT = Path("..").resolve()
DATA_DIR = ROOT / "data"
RESULTS_DIR = ROOT / "results"
FIG_DIR = RESULTS_DIR / "paper_figures"
FIG_DIR.mkdir(parents=True, exist_ok=True)

LOGFC_FILE = DATA_DIR / "merged_pval_filtered_deseq2_fillna.csv"
CORR_FILE = DATA_DIR / "correl_matrix1_merged_pval_filtered_deseq2.csv"
COLOR_MAP_FILE = DATA_DIR / "color_mapping.tsv"

print("Project root  :", ROOT)
print("Data dir      :", DATA_DIR)
print("Results dir   :", RESULTS_DIR)
print("Figures dir   :", FIG_DIR)

# %% [markdown]
# ## 1. Helper functions
#
# We define:
# - a parser to simplify condition names into `<drug>_<dose>`
# - a loader for the color mapping (drug → numeric color code)
# - basic loaders for log2FC matrix and correlation matrix

# %%
def condition_to_drug_label(cond_name: str) -> str:
    """
    Reduce a condition name like:
      'GefitinibA_006u_exp200921_run1_sample1'
    to something like:
      'Gefitinib_006u'
    
    Controls and time zeros get labels like:
      'CtrlMs_000u', 'Temps0_000u'
    """
    parts = cond_name.split("_")
    if len(parts) < 2:
        return cond_name

    head = parts[0]
    dose = parts[1]

    # Controls
    if head.startswith("Ctrl") or "Contro" in head:
        drug = "CtrlMs"
    # Time zeros
    elif head.startswith("Temps0") or "Temps0" in head:
        drug = "Temps0"
    else:
        # DrugA / DrugB → drug = all but last char, replicate = last char
        drug = head[:-1]

    return f"{drug}_{dose}"


def load_color_mapping(path: Path) -> dict:
    """
    color_mapping.tsv:
       <int_color_code>\t<drug_label>
    Example:
       1   Gefitinib_006u
    
    Returns:
        dict: drug_name (without replicate, dose included) → float in [0, 1]
    """
    cmap = {}

    if not path.exists():
        print(f"[info] Color mapping file not found: {path}")
        return cmap

    with path.open("r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) != 2:
                continue
            color_int, drug_label = parts
            color_int = int(color_int)

            # Harmonize special cases (as in networks.py)
            if drug_label == "Fluor_006u":
                drug_key = "5Fluor_006u"
            elif drug_label == "Azacyt_1,5u":
                drug_key = "Azacyt_1.5u"
            elif drug_label == "Bafilo_1,2n":
                drug_key = "Bafilo_1.2n"
            else:
                drug_key = drug_label

            # Scale to [0, 1] for colormap
            cmap[drug_key] = (color_int + 1) / 100.0

    print(f"[info] Loaded {len(cmap)} entries from color mapping.")
    return cmap


def load_logfc_matrix(path: Path) -> pd.DataFrame:
    """
    Load log2FC matrix (barcodes x conditions) from DESeq2 output.
    Assumes rows = barcodes, columns = conditions.
    """
    df = pd.read_csv(path, sep=";", header=0, index_col=0)
    print("[log2FC] Loaded matrix with shape:", df.shape)
    df = df.astype(float)
    df = df.dropna(axis=0, how="all").dropna(axis=1, how="all")

    # Safety: check for duplicate columns
    if df.columns.duplicated().any():
        dup = df.columns[df.columns.duplicated()].tolist()
        raise ValueError(f"Duplicate condition columns in logFC matrix: {dup}")

    print("[log2FC] After dropping all-NaN rows/cols:", df.shape)
    return df


def load_corr_matrix(logfc_df: pd.DataFrame, corr_file: Path) -> pd.DataFrame:
    """
    Load correlation matrix from file if present,
    otherwise compute Pearson correlation on the logFC matrix.
    """
    if corr_file.exists():
        corr = pd.read_csv(corr_file, sep=";", header=0, index_col=0)
        corr = corr.astype(float)
        print("[corr] Loaded precomputed correlation matrix:", corr.shape)
    else:
        print("[corr] Computing Pearson correlation from logFC matrix...")
        corr = logfc_df.corr(method="pearson", min_periods=1)
        print("[corr] Computed correlation matrix:", corr.shape)

    # Ensure perfect diagonal
    np.fill_diagonal(corr.values, 1.0)
    return corr


# %% [markdown]
# ## 2. Load data
#
# Here we load:
# - the log2 fold-change matrix (DESeq2 results merged across conditions)
# - the correlation matrix between conditions
# - an optional drug color mapping

# %%
logfc = load_logfc_matrix(LOGFC_FILE)
corr = load_corr_matrix(logfc, CORR_FILE)
color_mapping = load_color_mapping(COLOR_MAP_FILE)

logfc.head()

# %% [markdown]
# Quick sanity checks: basic statistics of log2FC values and fraction of missing values.

# %%
logfc_values = logfc.values.flatten()
nan_mask = np.isnan(logfc_values)
print("Total values      :", logfc_values.size)
print("NaN values        :", nan_mask.sum())
print("Non-NaN values    :", (~nan_mask).sum())

print("\nSummary statistics (non-NaN):")
print(pd.Series(logfc_values[~nan_mask]).describe())

# %% [markdown]
# ### Distribution of log2FC values
#
# We look at the global distribution (all barcodes, all conditions).

# %%
plt.figure(figsize=(6, 4))
sns.histplot(logfc_values[~nan_mask], bins=100, kde=True)
plt.xlabel("log2 fold-change")
plt.ylabel("Count")
plt.title("Global distribution of log2FC (all barcodes, all conditions)")
plt.tight_layout()
plt.show()

# %% [markdown]
# You can optionally clip extreme values to see the main bulk more clearly.

# %%
plt.figure(figsize=(6, 4))
clipped_vals = np.clip(logfc_values[~nan_mask], -5, 5)
sns.histplot(clipped_vals, bins=100, kde=True)
plt.xlabel("log2 fold-change (clipped to [-5, 5])")
plt.ylabel("Count")
plt.title("Distribution of log2FC (clipped)")
plt.tight_layout()
plt.show()

# %% [markdown]
# ## 3. Figure 3 – Drug clustering based on barcode signatures
#
# Each column is a condition; each row is a barcode.
# We cluster both barcodes and conditions based on log2FC profiles.
#
# To keep the figure readable, we:
# - clip log2FC values to [-4, 4]
# - optionally restrict to the most variable barcodes

# %%
# Compute per-barcode variance across conditions
barcode_var = logfc.var(axis=1, skipna=True)
barcode_var.describe()

# %%
# Select a subset of most variable barcodes for plotting
# (otherwise the heatmap can be huge)
N_BARCODES = 3000  # adjust as needed (e.g. 1000, 5000)
top_barcodes = barcode_var.sort_values(ascending=False).head(N_BARCODES).index

logfc_subset = logfc.loc[top_barcodes]
logfc_subset_clipped = logfc_subset.clip(lower=-4, upper=4)

logfc_subset_clipped.shape

# %%
# Clustered heatmap (Figure 3 style)
sns.set(font_scale=0.6)

g = sns.clustermap(
    logfc_subset_clipped,
    method="complete",
    metric="euclidean",
    cmap="RdBu_r",
    center=0.0,
    xticklabels=False,
    yticklabels=False,
    robust=False,
    figsize=(10, 12),
)

g.fig.suptitle(
    f"Drug clustering based on barcode signatures (top {N_BARCODES} variable barcodes)",
    y=1.02,
)

fig3_pdf = FIG_DIR / "Fig3_drug_signature_heatmap_topBarcodes.pdf"
fig3_png = FIG_DIR / "Fig3_drug_signature_heatmap_topBarcodes.png"

g.fig.savefig(fig3_pdf, bbox_inches="tight")
g.fig.savefig(fig3_png, dpi=300, bbox_inches="tight")

plt.show()

print("Saved:", fig3_pdf)
print("Saved:", fig3_png)

# %% [markdown]
# ### Sanity check: which drugs / conditions cluster together?
#
# We can inspect the ordered column labels from the clustermap.

# %%
ordered_conditions = g.dendrogram_col.reordered_ind
ordered_cols = [logfc_subset_clipped.columns[i] for i in ordered_conditions]
ordered_cols[:20]

# %% [markdown]
# You can search for specific drugs (e.g. EGFR inhibitors, chemotherapies) in the ordered condition list.

# %%
# Example: find positions of some manually chosen drugs (adjust names to your dataset)
keywords = ["Gefitinib", "Lazertinib", "Osimertinib", "Carboplatin", "Cisplatin"]

for kw in keywords:
    matches = [c for c in ordered_cols if kw in c]
    print(f"\nKeyword: {kw}")
    for m in matches:
        print("  ", m)

# %% [markdown]
# ## 4. Figure 4 – Drug–drug correlation clustered heatmap
#
# Each cell shows the Pearson correlation between two conditions
# (based on the log2FC profiles across barcodes).
#
# - Blue: correlation = 1
# - Red: correlation = -1
# - White: correlation ~ 0

# %%
sns.set(font_scale=0.5)

g_corr = sns.clustermap(
    corr,
    method="complete",
    metric="euclidean",
    cmap="RdBu_r",
    vmin=-1.0,
    vmax=1.0,
    center=0.0,
    xticklabels=False,
    yticklabels=False,
    figsize=(10, 10),
)

g_corr.fig.suptitle("Drug–drug correlation clustered heatmap (Pearson r)", y=1.02)

fig4_pdf = FIG_DIR / "Fig4_drug_correlation_heatmap.pdf"
fig4_png = FIG_DIR / "Fig4_drug_correlation_heatmap.png"

g_corr.fig.savefig(fig4_pdf, bbox_inches="tight")
g_corr.fig.savefig(fig4_png, dpi=300, bbox_inches="tight")

plt.show()

print("Saved:", fig4_pdf)
print("Saved:", fig4_png)

# %% [markdown]
# ### Inspect correlation of a given condition
#
# Helper: display the conditions most correlated with a selected condition.

# %%
def top_correlated_conditions(corr_matrix: pd.DataFrame, condition: str, n: int = 10):
    if condition not in corr_matrix.columns:
        raise ValueError(f"Condition '{condition}' not in correlation matrix.")
    series = corr_matrix[condition].dropna().sort_values(ascending=False)
    return series.head(n)


# Example: replace with a real condition name from your data
example_condition = corr.columns[0]
print("Example condition:", example_condition)
top_correlated_conditions(corr, example_condition, n=10)

# %% [markdown]
# ## 5. Figure 5 – Drug network based on correlation
#
# We build a graph where:
#
# - nodes = conditions
# - edges connect conditions with Pearson correlation ≥ threshold (default 0.8)
# - node colors reflect drug annotation (if `color_mapping.tsv` is available)
#
# Positioning is done with a Fruchterman–Reingold force-directed layout
# (`nx.spring_layout`).

# %%
CORR_THRESHOLD = 0.8

def build_correlation_graph(corr_matrix: pd.DataFrame, threshold: float = 0.8) -> nx.Graph:
    """
    Build graph from correlation matrix:
    - Add edge between (i, j) if corr(i, j) >= threshold
    - Ignore self-correlations
    """
    G = nx.Graph()
    conditions = corr_matrix.index.tolist()

    for i, s in enumerate(conditions):
        for j in range(i + 1, len(conditions)):
            t = conditions[j]
            r = corr_matrix.iloc[i, j]
            if np.isnan(r):
                continue
            if r >= threshold:
                G.add_edge(s, t, weight=float(r))

    return G


G = build_correlation_graph(corr, threshold=CORR_THRESHOLD)
print("Number of nodes:", G.number_of_nodes())
print("Number of edges:", G.number_of_edges())

# %%
# Build node colors using drug labels + color_mapping
node_colors = []
for n in G.nodes():
    drug_label = condition_to_drug_label(n)
    c = color_mapping.get(drug_label, 0.5)  # default mid grey if not found
    node_colors.append(c)

# Node labels: keep only the drug name without dose for readability
labels = {n: condition_to_drug_label(n).split("_")[0] for n in G.nodes()}

# Compute layout
pos = nx.spring_layout(G, k=0.15, iterations=200, seed=42)

plt.figure(figsize=(10, 10))
nx.draw_networkx_edges(
    G,
    pos,
    edgelist=G.edges(),
    width=0.5,
    alpha=0.6,
    edge_color="black",
)

nodes = nx.draw_networkx_nodes(
    G,
    pos,
    node_size=300,
    node_color=node_colors,
    cmap="viridis",
    alpha=0.9,
)

nx.draw_networkx_labels(G, pos, labels=labels, font_size=6)
plt.title(f"Drug network (Pearson r ≥ {CORR_THRESHOLD})")
plt.axis("off")

fig5_pdf = FIG_DIR / "Fig5_drug_network.pdf"
fig5_png = FIG_DIR / "Fig5_drug_network.png"

plt.savefig(fig5_pdf, bbox_inches="tight")
plt.savefig(fig5_png, dpi=300, bbox_inches="tight")
plt.show()

print("Saved:", fig5_pdf)
print("Saved:", fig5_png)

# %% [markdown]
# ### Optional: legend for drug classes
#
# If your `color_mapping.tsv` encodes drug *classes* rather than individual drugs,
# you can build a small legend by hand (e.g. manually mapping class → color).
# This is project-specific and can be customized later.
#
# For now, you can inspect the mapping between node labels and numerical colors:

# %%
drug_label_to_color = {}
for n in G.nodes():
    drug_label = condition_to_drug_label(n)
    drug_label_to_color[drug_label] = color_mapping.get(drug_label, 0.5)

list(drug_label_to_color.items())[:20]

# %% [markdown]
# ## 6. Summary
#
# In this notebook, we:
#
# - Loaded DESeq2-derived log2FC per barcode and per condition
# - Visualized barcode signatures across drugs (heatmap + clustering)
# - Computed and visualized drug–drug correlations
# - Built a correlation-based drug network with node colors based on annotation
#
# This notebook is a good "showcase" complement to the more scripted pipeline,
# and can be linked from your project README as an illustrative entry point.

