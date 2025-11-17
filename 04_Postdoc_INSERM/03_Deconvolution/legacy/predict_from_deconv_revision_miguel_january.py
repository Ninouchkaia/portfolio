# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.13.0
#   kernelspec:
#     display_name: ml_interactive
#     language: python
#     name: ml_interactive
# ---

# %% [markdown]
# # Predict from deconvolution - revision

# %% [markdown]
# Train ElasticNet penalized Logistic regression on deconvolution data
# only to assess how important / usefull new signatures are.

# %% [markdown]
# ## Imports

# %%
import pandas as pd
import numpy as np
from pathlib import Path

from sklearn.model_selection import train_test_split
import sklearn.metrics as metrics

from sklearn.pipeline import Pipeline
from sklearn.linear_model import LogisticRegression

# from sklearn.feature_selection import SelectKBest
from sklearn.linear_model import LogisticRegressionCV

from sklearn.preprocessing import StandardScaler

import pprint

from tqdm import tqdm


# %%
def write_log(file, message):
    print(message)
    file.write(message + "\n")


# %% [markdown]
# ## Data loading and preprocessing

# %%
ONLY_MELANOMA = True
DROP_SCORES = (
    True  # discard XCELL's 'ImmuneScore', 'StromaScore' and 'MicroenvironmentScore'
)
TEST_CLOUG = True  # models trained on all datasets are tested on the Cloughesy dataset
CV_ADAPT = False  # increase k, the number of CV folds, if training didn't converged with previous k
CV_TRAIN = 5  # default number of splits for CV
CV_MAX = 10  # max splits if previous training failed

path_interm = Path("../data/intermediate/Deconvolution_paper")
data_version = "all_deconvolutions_2021-08-10"
data_version = "all_deconvolutions_2021-08-30"
data_version = "all_deconvolutions_2021-09-16"
data_version = "all_deconvolutions_2021-10-15"
data_version = "all_deconvolutions_2021-10-19"
data_version = "all_deconvolutions_2021-11-02"
path_deconv = path_interm / "revision_deconv" / data_version
dir_save = Path("../data/processed/Deconvolution_paper_revision")
dir_save = dir_save / "data_version-{}".format(data_version)

if ONLY_MELANOMA:
    dir_save = dir_save / "only_melanoma"
else:
    dir_save = dir_save / "all_cancer_types"
if DROP_SCORES:
    dir_save = dir_save / "drop_scores"
else:
    dir_save = dir_save / "keep_scores"
if CV_ADAPT:
    dir_save = dir_save / "cv_adapt-{}-{}".format(CV_TRAIN, CV_MAX)
else:
    dir_save = dir_save / "cv-{}".format(CV_TRAIN)
dir_save_orig = dir_save
dir_save.mkdir(parents=True, exist_ok=True)

response = pd.read_csv(path_interm / "response_train.csv", index_col=0)

datasets = ["Gide", "Riaz", "Hugo"]  # , 'Cloughesy']
deconv = None
for dataset in datasets:
    path_data = path_deconv / ("all_deconvolutions_" + dataset + ".txt")
    if deconv is None:
        deconv = pd.read_csv(path_data, sep="\t", index_col=0)
    else:
        new_deconv = pd.read_csv(path_data, sep="\t", index_col=0)
        deconv = deconv.append(new_deconv)

common_ids = set(response.index).intersection(deconv.index)
df_all = deconv.loc[common_ids, :]
response = response.loc[df_all.index, :]

# drop variables with unique values
drop_col = [x for x in df_all.columns if df_all[x].unique().size == 1]
if len(drop_col) > 0:
    print("dropping these columns because they have a unique value:")
else:
    print("there is no column with unique value")
for i in drop_col:
    print("    ", i)
matches = []  # ['XCELL', 'TAMs', 'CD226']
if DROP_SCORES:
    matches = matches + ["Score"]
drop_col = drop_col + [x for x in df_all.columns if any(y in x for y in matches)]
if len(drop_col) > 0:
    df_all.drop(columns=drop_col, inplace=True)
# drop Maha samples
drop_row = [x for x in df_all.index if "Maha" in x]
if ONLY_MELANOMA:
    drop_row = drop_row + [x for x in df_all.index if "Snyder" in x]
    drop_row = drop_row + [x for x in df_all.index if "Cloughesy" in x]
df_all.drop(labels=drop_row, axis=0, inplace=True)
response.drop(labels=drop_row, inplace=True)
# set inf values to nan
to_nan = ~np.isfinite(df_all.values)
nb_nan = to_nan.sum()
if nb_nan != 0:
    print(f"There are {nb_nan} nan values")
    df_all[to_nan] = np.nan
    # impute non finite values (nan, +/-inf)
    # imputer = IterativeImputer(max_iter=100, random_state=0)
    imputer = KNNImputer(n_neighbors=5, weights="distance")
    df_all.loc[:, :] = imputer.fit_transform(df_all.values)
common_col = df_all.columns.values
nb_var = df_all.shape[1]

# rename variables
if data_version == "all_deconvolutions_2021-11-02":
    dic_rename = {
        "Epidish": "EpiDISH",
        "Quantiseq": "quanTIseq",
        "CBSX__": "CIBERSORTx_CBSX_",
        "CBSX_scRNA-Seq_melanoma_Tirosh_sigmatrix_SuppFig_3-b": "CBSX_melanoma",
        "CBSX_sigmatrix_HNSCC_Fig2cd": "CBSX_HNSCC",
        "CBSX_Fig2ab-NSCLC_PBMCs_scRNAseq_sigmatrix": "CBSX_NSCLC",
        "__": "_",
        # 'CBSX_LM22': 'CBSX_LM22', new signature, already well mentioned
    }
else:
    dic_rename = {
        "Epidish": "EpiDISH",
        "Quantiseq": "quanTIseq",
        "CBSX_Fig2ab-NSCLC_PBMCs_scRNAseq_refsample_single-cell": "CBSX_NSCLC",
        "CBSX_scRNA-Seq_reference_HNSCC_Puram_et_al_Fig2cd_single-cell": "CBSX_HNSCC",
        "CBSX_scRNA-Seq_reference_melanoma_Tirosh_SuppFig_3b-d_single-cell": "CBSX_melanoma",
    }
for key, val in dic_rename.items():
    new_cols = [x.replace(key, val) for x in df_all.columns]
    df_all.columns = new_cols
# add CIBERSORTx method name where it's missing
# not needed with dic_rename's key CBSX__
# new_cols = ['CIBERSORTx_' + x if ('CBSX' in x and not ('EpiDISH' in x or 'DeconRNASeq' in x)) else x for x in df_all.columns]
# df_all.columns = new_cols

if TEST_CLOUG:
    path_data = path_deconv / ("all_deconvolutions_" + "Cloughesy" + ".txt")
    deconv_cloug = pd.read_csv(path_data, sep="\t", index_col=0)
    deconv_cloug.index = ["Cloughesy_" + x for x in deconv_cloug.index]
    response_cloug = pd.read_csv(path_interm / "response_train.csv", index_col=0)
    common_ids_cloug = set(response_cloug.index).intersection(deconv_cloug.index)
    deconv_cloug = deconv_cloug.loc[common_ids_cloug, :]
    response_cloug = response_cloug.loc[deconv_cloug.index, :]
    if len(drop_col) > 0:
        deconv_cloug.drop(columns=drop_col, inplace=True)
    # set inf values to nan
    to_nan = ~np.isfinite(deconv_cloug.values)
    nb_nan = to_nan.sum()
    if nb_nan != 0:
        print(f"There are {nb_nan} nan values")
        deconv_cloug[to_nan] = np.nan
        # impute non finite values (nan, +/-inf)
        # imputer = IterativeImputer(max_iter=100, random_state=0)
        imputer = KNNImputer(n_neighbors=5, weights="distance")
        deconv_cloug.loc[:, :] = imputer.fit_transform(deconv_cloug.values)
    for key, val in dic_rename.items():
        new_cols = [x.replace(key, val) for x in deconv_cloug.columns]
        deconv_cloug.columns = new_cols
    # add CIBERSORTx method name where it's missing
    # new_cols = ['CIBERSORTx_' + x if ('CBSX' in x and not ('EpiDISH' in x or 'DeconRNASeq' in x)) else x for x in deconv_cloug.columns]
    # deconv_cloug.columns = new_cols

# %% [markdown]
# ## Setup training parameters

# %%
# we use trailing underscores to avoid including derived signatures
# like EpiDISH_BPRNACan3Dprom --> EpiDISH_BPRNACan3Dprom-enhan

if data_version == "all_deconvolutions_2021-11-02":
    signatures = [
        "quanTIseq",
        "MCP",
        "XCELL",
        "EpiDISH_BPRNACan_",
        "DeconRNASeq_BPRNACan_",
        "EpiDISH_BPRNACanProMet_",
        "DeconRNASeq_BPRNACanProMet_",
        "EpiDISH_BPRNACan3DProMet_",
        "DeconRNASeq_BPRNACan3DProMet_",
        "EpiDISH_CBSX_NSCLC_",
        "DeconRNASeq_CBSX_NSCLC_",
        "EpiDISH_CBSX_HNSCC_",
        "DeconRNASeq_CBSX_HNSCC_",
        "EpiDISH_CBSX_melanoma_",
        "DeconRNASeq_CBSX_melanoma_",
        "EpiDISH_CBSX_LM22_",
        "DeconRNASeq_CBSX_LM22_",
        "CIBERSORTx_CBSX_NSCLC_",
        "CIBERSORTx_CBSX_HNSCC_",
        "CIBERSORTx_CBSX_melanoma_",
        "CIBERSORTx_CBSX_LM22_",
    ]
else:
    signatures = [
        "quanTIseq",
        "MCP",
        "XCELL",
        "EpiDISH_BPRNACan_",
        "DeconRNASeq_BPRNACan_",
        "EpiDISH_BPRNACanProMet_",
        "DeconRNASeq_BPRNACanProMet_",
        "EpiDISH_BPRNACan3DProMet_",
        "DeconRNASeq_BPRNACan3DProMet_",
        "EpiDISH_CBSX_NSCLC_",
        "DeconRNASeq_CBSX_NSCLC_",
        "EpiDISH_CBSX_HNSCC_",
        "DeconRNASeq_CBSX_HNSCC_",
        "EpiDISH_CBSX_melanoma_",
        "DeconRNASeq_CBSX_melanoma_",
        "CIBERSORTx_CBSX_NSCLC_",
        "CIBERSORTx_CBSX_HNSCC_",
        "CIBERSORTx_CBSX_melanoma_",
    ]

if ONLY_MELANOMA or MERGE_OLD_SIG:
    datasets = ["Gide", "Hugo", "Riaz"]
else:
    datasets = ["Gide", "Hugo", "Riaz", "Snyder"]
# l1_ratio = 0 the penalty is an L2 penalty [Ridge]
# l1_ratio = 1 the penalty is an L1 penalty (Lasso)

# %% [markdown]
# # Exploratory data analysis

# %% [markdown]
# We look at how clinical variables and cell types proportions correlate with response.

# %% [markdown]
# ## Clinical data

# %%
path_clin = path_interm / "clinic_train.csv"
clin = pd.read_csv(path_clin, index_col=0)
drop_rows = [x for x in clin.index if x.startswith("Cloug")]
clin.drop(index=drop_rows, inplace=True)

# %%
# df_clinical_data = response.join(clin)
df_clinical_data = clin

# %%
# print(df)

target_clinical_columns = ["TMB", "SEX", "AAGE"]

np.random.seed(0)

score_labels = [
    "ROC AUC",  # Receiver Operating Characteristic Area Under the Curve
    "AP",  # Average Precision
    "MCC",  # Matthews Correlation Coefficient
]


def one_column_roc_auc_score(df_clinical_data, x_column, y_column="BOR - binary"):

    df = response.join(df_clinical_data)
    non_na_samples = ~df[x_column].isna()

    X = df.loc[:, x_column].values
    y = response["BOR - binary"].values

    if x_column == "SEX":
        X[X == "M"] = 0
        X[X == "F"] = 1

    X = X[non_na_samples].reshape(-1, 1)
    y = y[non_na_samples]

    X_train, X_test, y_train, y_test = train_test_split(
        X,
        y,
        test_size=0.25,
        random_state=0,
        shuffle=True,
    )
    # reverse labels to match "equaling progressive disease" objective
    y_train = -(y_train - 1)
    y_test = -(y_test - 1)

    pipe = Pipeline(
        steps=[
            # ('select', SelectKBest(k=1)),
            ("clf", LogisticRegression())
        ]
    )

    pipe.fit(X_train, y_train)

    y_pred = pipe.predict_proba(X_test)[:, 1]
    score = metrics.roc_auc_score(y_test, y_pred)
    return score


scores = list(
    map(
        lambda x_column: one_column_roc_auc_score(df_clinical_data, x_column),
        target_clinical_columns,
    )
)
score_pairs = dict(zip(target_clinical_columns, scores))

for target_clinical_column, score in score_pairs.items():
    print(
        f"ROC AUC model for 'BOR - binary' trained only on {target_clinical_column}: {score}"
    )


l1_ratios_list = [
    ["default", [0.5]],
    ["naive", np.linspace(0, 1, 21)],  # naive param grid
    [
        "advised",
        [0.1, 0.5, 0.7, 0.9, 0.95, 0.99, 1],
    ],  # advised in scikit-learn documentation
]


def deconvolution_roc_auc_score(
    df_deconvolutions,
    df_clinical_data,
    signature_beggining,
    clinical_column,
    join_clinical_data,
    l1_ratios,
    y_column="BOR - binary",
):
    all_signatures_index = [
        x for x in df_deconvolutions.columns if x.startswith(signature_beggining)
    ]

    good_index = list(
        set(
            df_clinical_data.index[~df_clinical_data.loc[:, clinical_column].isna()]
        ).intersection(set(df_all.index))
    )
    # merge clinical data to deconv variables
    selected_clinical_data = df_clinical_data.loc[good_index, clinical_column].copy()
    if clinical_column == "SEX":
        selected_clinical_data = selected_clinical_data.map({"M": 0, "F": 1})
    if join_clinical_data:
        X = (
            df_deconvolutions.loc[good_index, all_signatures_index]
            .join(selected_clinical_data)
            .values
        )
    else:
        X = df_deconvolutions.loc[good_index, all_signatures_index].values

    y = response["BOR - binary"].loc[good_index].values

    # reverse labels to match "equaling progressive disease" objective
    y = -(y - 1)
    # stratify train / test by response
    X_train, X_test, y_train, y_test = train_test_split(
        X,
        y,
        test_size=0.25,
        random_state=0,
        shuffle=True,
    )
    # print(X)

    cv_used = CV_TRAIN
    training_succeeded = False

    while not training_succeeded and cv_used <= CV_MAX:
        pipe = Pipeline(
            steps=[
                ("scaler", StandardScaler()),
                (
                    "clf",
                    LogisticRegressionCV(
                        cv=cv_used,
                        Cs=20,
                        penalty="elasticnet",
                        # scoring='neg_log_loss',
                        scoring="roc_auc",
                        solver="saga",
                        l1_ratios=l1_ratios,
                        max_iter=10000,
                        n_jobs=-1,  # or n_jobs-1 to leave one core available
                    ),
                ),
            ]
        )

        pipe.fit(X_train, y_train)
        training_succeeded = not np.all(pipe.named_steps["clf"].coef_ == 0)

        if not training_succeeded:
            if CV_ADAPT:
                cv_used += 1
            else:
                break

    if training_succeeded:
        y_pred_proba = pipe.predict_proba(X_test)[:, 1]
        y_pred = pipe.predict(X_test)
        score = {
            "ROC AUC": metrics.roc_auc_score(y_test, y_pred_proba),
            "AP": metrics.average_precision_score(y_test, y_pred_proba),
            "MCC": metrics.matthews_corrcoef(y_test, y_pred),
        }
    else:
        score = {
            "ROC AUC": np.nan,
            "AP": np.nan,
            "MCC": np.nan,
        }

    score_split = {}

    score_column = clinical_column
    if not join_clinical_data:
        score_column += " deconv only"

    for metric in score_labels:
        score_split[score_column + " - " + metric] = score[metric]

    return score_split


df_deconvolutions = df_all


all_scores = {}

# score_with_clinical_data = deconvolution_roc_auc_score(
#     df_deconvolutions, df_clinical_data, signatures[0], "TMB", True, l1_ratios_list[0][1]
# )
# score_only_deconvolution = deconvolution_roc_auc_score(
#     df_deconvolutions, df_clinical_data, signatures[0], "TMB", False, l1_ratios_list[0][1]
# )
# all_scores[signatures[0]] = {**score_with_clinical_data , **score_only_deconvolution}
# print(all_scores)

# Test either one of those combinations
l1_ratios_list = [
    ["default", [0.5]],
    ["naive", np.linspace(0, 1, 21)],  # naive param grid
    [
        "advised",
        [0.1, 0.5, 0.7, 0.9, 0.95, 0.99, 1],
    ],  # advised in scikit-learn documentation
]

# print("Doing l1_ratios models...")
# for l1_name, l1_ratios in tqdm(l1_ratios_list):
#     print(l1_name)
#     for sig in tqdm(signatures):
#         signature_columns = [x for x in df_all.columns if x.startswith(sig)]
#         str_sig = sig.strip("_")
#         # print(str_sig)
#         for clinical_column in target_clinical_columns:
#             score_with_clinical_data = deconvolution_roc_auc_score(
#                 df_deconvolutions,
#                 df_clinical_data,
#                 str_sig,
#                 clinical_column,
#                 True,
#                 l1_ratios,
#             )
#             score_only_deconvolution = deconvolution_roc_auc_score(
#                 df_deconvolutions,
#                 df_clinical_data,
#                 str_sig,
#                 clinical_column,
#                 False,
#                 l1_ratios,
#             )
#
#             if not str_sig in all_scores:
#                 all_scores[str_sig] = {}
#
#             if not l1_name in all_scores[str_sig]:
#                 all_scores[str_sig][l1_name] = {}
#
#             all_scores[str_sig][l1_name] = {**all_scores[str_sig][l1_name], **score_with_clinical_data , **score_only_deconvolution}
#
# pprint.pprint(all_scores)


from scipy import stats

y = response['BOR - binary'].values
# reverse labels to match "equaling progressive disease" objective
y = -(y-1)

rhos = []
pvals = []
for col in df_deconvolutions.columns:
    rho, pval = stats.spearmanr(df_deconvolutions[col], y)
    rhos.append(rho)
    pvals.append(pval)

# %%
from statsmodels.stats.multitest import fdrcorrection

df_corr = pd.DataFrame({'rho':rhos, 'pval': pvals}, index=df_deconvolutions.columns)
df_corr['abs rho'] = df_corr['rho'].abs()
for sig in signatures:
    var_idx = [x for x in df_corr.index if x.startswith(sig)]
    rejected, pval_corr = fdrcorrection(df_corr.loc[var_idx, 'pval'], method='indep')
    df_corr.loc[var_idx, 'pval_corr'] = pval_corr
df_corr.sort_values(by=['abs rho'], ascending=False, inplace=True)

df_corr.loc[df_corr['pval'] <= 0.05]

# %%
df_corr.loc[df_corr['pval_corr'] <= 0.05]

print(df_corr.loc[df_corr['pval_corr'] <= 0.05])

df_corr.to_csv('df_corr_miguel_january.csv', sep=';', index=True)


# %%
counts = []
for sig in signatures:
    var_idx = [x for x in df_all.columns if x.startswith(sig)]
    counts.append(len(var_idx))
sig_counts = pd.DataFrame({'# cell types': counts}, index=signatures)
sig_counts.index.name = 'signatures'
sig_counts
# print(sig_counts)


def deconvolution_only_roc_auc_score(
    df_deconvolutions,
    signature_beggining,
    dataset_column,
    l1_ratios,
    y_column="BOR - binary",
):
    test_ids = [x for x in df_deconvolutions.index if dataset_column in x]
    train_ids = [x for x in df_deconvolutions.index if dataset_column not in x]
    X_train = df_deconvolutions.loc[train_ids, var_idx].values
    X_test= df_deconvolutions.loc[test_ids, var_idx].values
    y_train = response[y_column].loc[train_ids].values
    y_test = response[y_column].loc[test_ids].values
    # reverse labels to match "equaling progressive disease" objective
    y_train = -(y_train-1)
    y_test = -(y_test-1)

    training_succeeded = False
    cv_used = CV_TRAIN

    while not training_succeeded and cv_used <= CV_MAX:
        pipe = Pipeline(
            steps=[
                ("scaler", StandardScaler()),
                (
                    "clf",
                    LogisticRegressionCV(
                        cv=cv_used,
                        Cs=20,
                        penalty="elasticnet",
                        # scoring='neg_log_loss',
                        scoring="roc_auc",
                        solver="saga",
                        l1_ratios=l1_ratios,
                        max_iter=10000,
                        n_jobs=-1,  # or n_jobs-1 to leave one core available
                    ),
                ),
            ]
        )

        pipe.fit(X_train, y_train)
        training_succeeded = not np.all(pipe.named_steps["clf"].coef_ == 0)

        if not training_succeeded:
            if CV_ADAPT:
                cv_used += 1
            else:
                break

    if training_succeeded:
        y_pred_proba = pipe.predict_proba(X_test)[:, 1]
        y_pred = pipe.predict(X_test)
        score = {
            "ROC AUC": metrics.roc_auc_score(y_test, y_pred_proba),
            "AP": metrics.average_precision_score(y_test, y_pred_proba),
            "MCC": metrics.matthews_corrcoef(y_test, y_pred),
        }
    else:
        score = {
            "ROC AUC": np.nan,
            "AP": np.nan,
            "MCC": np.nan,
        }

    score_split = {}
    for metric in score_labels:
        score_split[dataset_column + " - " + metric] = score[metric]

    return score_split


print("Doing l1_ratios LODO models...")
all_deconvolution_only_scores = {}
for l1_name, l1_ratios in tqdm(l1_ratios_list):
    print(l1_name)
    for sig in tqdm(signatures):
        signature_columns = [x for x in df_deconvolutions.columns if x.startswith(sig)]
        str_sig = sig.strip("_")
        for dataset_column in tqdm(datasets):
            score_deconvolution_only_clinical_data = deconvolution_only_roc_auc_score(
                df_deconvolutions,
                str_sig,
                dataset_column,
                l1_ratios,
            )

            if not str_sig in all_deconvolution_only_scores:
                all_deconvolution_only_scores[str_sig] = {}

            if not l1_name in all_deconvolution_only_scores[str_sig]:
                all_deconvolution_only_scores[str_sig][l1_name] = {}

            all_deconvolution_only_scores[str_sig][l1_name] = {**all_deconvolution_only_scores[str_sig][l1_name],**score_deconvolution_only_clinical_data}

pprint.pprint(all_deconvolution_only_scores)
