import torch
import anndata as ad
import pandas as pd
import numpy as np

from scvi.model import SCANVI


# ============================================================
# PATHS
# ============================================================

MODEL_FILE = (
    "/storage/lemus_g/roldan/ARMD/"
    "HRCA_reference/model.pt"
)

QUERY_FILE = (
    "/storage/lemus_g/roldan/ARMD/results/"
    "SRP413248/HRCA/"
    "SRP413248_HRCA_query_symbols_10000.h5ad"
)

OUTDIR = (
    "/storage/lemus_g/roldan/ARMD/results/"
    "SRP413248/HRCA"
)


# ============================================================
# 1. LOAD CHECKPOINT
# ============================================================

print("=" * 70)
print("1. Loading HRCA checkpoint")
print("=" * 70)

checkpoint = torch.load(
    MODEL_FILE,
    map_location="cpu",
    weights_only=False,
)

state_dict = checkpoint["model_state_dict"]
attr_dict = checkpoint["attr_dict"]
var_names = list(checkpoint["var_names"])

registry = attr_dict["registry_"]

setup_args = registry["setup_args"]


# ============================================================
# 2. REFERENCE CATEGORIES
# ============================================================

batch_registry = (
    registry["field_registries"]["batch"]
)

reference_batches = list(
    batch_registry[
        "state_registry"
    ]["categorical_mapping"]
)

label_registry = (
    registry["field_registries"]["labels"]
)

reference_labels = list(
    label_registry[
        "state_registry"
    ]["categorical_mapping"]
)


# ============================================================
# 3. LOAD QUERY
# ============================================================

print("\n" + "=" * 70)
print("2. Loading ARMD query")
print("=" * 70)

query = ad.read_h5ad(
    QUERY_FILE
)

print(
    "Cells:",
    query.n_obs
)

print(
    "Genes:",
    query.n_vars
)


# ============================================================
# 4. VERIFY GENES
# ============================================================

assert list(query.var_names) == var_names, (
    "Query genes do not exactly match "
    "the HRCA reference."
)

print(
    "Gene order matches HRCA reference."
)


# ============================================================
# 5. QUERY METADATA
# ============================================================

query.obs["celltype"] = "Unknown"

if "sample" in query.obs.columns:

    query.obs["sampleid"] = (
        query.obs["sample"]
        .astype(str)
    )

elif "GSM" in query.obs.columns:

    query.obs["sampleid"] = (
        query.obs["GSM"]
        .astype(str)
    )

else:

    raise ValueError(
        "No sample or GSM column found."
    )


print(
    "\nQuery batches:",
    query.obs[
        "sampleid"
    ].nunique()
)


# ============================================================
# 6. CREATE DUMMY DATA
# ============================================================

print("\n" + "=" * 70)
print("3. Reconstructing HRCA model")
print("=" * 70)

n_dummy = len(
    reference_batches
)


dummy_batches = pd.Categorical(
    reference_batches,
    categories=reference_batches,
)

dummy_labels = pd.Categorical(
    ["Unknown"] * n_dummy,
    categories=reference_labels,
)


dummy = ad.AnnData(

    X=query.X[
        :n_dummy,
        :
    ].copy(),

    obs=pd.DataFrame(

        {
            "sampleid":
                dummy_batches,

            "celltype":
                dummy_labels,
        },

        index=[
            f"HRCA_dummy_{i}"
            for i in range(n_dummy)
        ],
    ),

    var=query.var.copy(),
)


if "counts" in query.layers:

    dummy.layers["counts"] = (
        query.layers[
            "counts"
        ][
            :n_dummy,
            :
        ].copy()
    )


# ============================================================
# 7. SCANVI SETUP
# ============================================================

SCANVI.setup_anndata(

    dummy,

    labels_key=(
        setup_args[
            "labels_key"
        ]
    ),

    unlabeled_category=(
        setup_args[
            "unlabeled_category"
        ]
    ),

    batch_key=(
        setup_args[
            "batch_key"
        ]
    ),
)


# ============================================================
# 8. RECREATE EXACT ARCHITECTURE
# ============================================================

init_params = (
    attr_dict["init_params_"]
)

non_kwargs = dict(
    init_params[
        "non_kwargs"
    ]
)

model_kwargs = dict(
    init_params[
        "kwargs"
    ][
        "model_kwargs"
    ]
)


model = SCANVI(

    dummy,

    **non_kwargs,

    **model_kwargs,
)


# ============================================================
# 9. LOAD HRCA WEIGHTS
# ============================================================

model.module.load_state_dict(
    state_dict,
    strict=True,
)

print(
    "HRCA reference weights loaded."
)


# ============================================================
# 10. PREPARE QUERY
# ============================================================

print("\n" + "=" * 70)
print("4. Preparing ARMD query")
print("=" * 70)

query_prepared = (
    SCANVI.prepare_query_anndata(

        query,

        reference_model=model,

        inplace=False,
    )
)


print(
    "Prepared:",
    query_prepared.shape
)


# ============================================================
# 11. LOAD QUERY
# ============================================================

print("\n" + "=" * 70)
print("5. Loading query into HRCA")
print("=" * 70)

query_model = (
    SCANVI.load_query_data(

        query_prepared,

        model,

        freeze_dropout=True,
    )
)


print(
    "Query loaded successfully."
)


# ============================================================
# 12. PREDICT CELL TYPES
# ============================================================

print("\n" + "=" * 70)
print("6. Predicting HRCA cell types")
print("=" * 70)

predictions = (
    query_model.predict()
)

probabilities = (
    query_model.predict(
        soft=True
    )
)


# ============================================================
# 13. MAXIMUM PROBABILITY
# ============================================================

max_probability = (
    probabilities.max(
        axis=1
    )
)


# ============================================================
# 14. ADD ANNOTATIONS
# ============================================================

query.obs[
    "HRCA_celltype"
] = predictions.astype(str)

query.obs[
    "HRCA_prediction_probability"
] = max_probability.values


# ============================================================
# 15. SAVE PREDICTIONS
# ============================================================

prediction_file = (
    f"{OUTDIR}/"
    "SRP413248_HRCA_predictions.csv"
)

query.obs[
    [
        "sample",
        "GSM",
        "Disease",
        "Age",
        "HRCA_celltype",
        "HRCA_prediction_probability",
    ]
].to_csv(
    prediction_file
)


# ============================================================
# 16. SAVE FULL ANNOTATED H5AD
# ============================================================

annotated_file = (
    f"{OUTDIR}/"
    "SRP413248_HRCA_annotated.h5ad"
)

query.write(
    annotated_file
)


# ============================================================
# 17. PRINT RESULTS
# ============================================================

print("\n" + "=" * 70)
print("HRCA ANNOTATION RESULTS")
print("=" * 70)

print(
    "\nTotal cells:",
    query.n_obs
)

print(
    "\nCell-type predictions:"
)

print(
    query.obs[
        "HRCA_celltype"
    ].value_counts()
)


print(
    "\nPrediction probability:"
)

print(
    query.obs[
        "HRCA_prediction_probability"
    ].describe()
)


# ============================================================
# 18. SAVE CELL-TYPE COUNTS
# ============================================================

counts_file = (
    f"{OUTDIR}/"
    "SRP413248_HRCA_celltype_counts.csv"
)

celltype_counts = (
    query.obs[
        "HRCA_celltype"
    ]
    .value_counts()
    .rename_axis(
        "HRCA_celltype"
    )
    .reset_index(
        name="n_cells"
    )
)

celltype_counts[
    "percentage"
] = (
    celltype_counts["n_cells"]
    / query.n_obs
    * 100
)

celltype_counts.to_csv(
    counts_file,
    index=False
)


# ============================================================
# 19. SAVE PROBABILITIES
# ============================================================

probability_file = (
    f"{OUTDIR}/"
    "SRP413248_HRCA_prediction_probabilities.csv"
)

probabilities.to_csv(
    probability_file
)


# ============================================================
# DONE
# ============================================================

print("\n" + "=" * 70)
print("DONE")
print("=" * 70)

print(
    "\nPredictions:"
)

print(
    prediction_file
)

print(
    "\nAnnotated H5AD:"
)

print(
    annotated_file
)

print(
    "\nCell-type counts:"
)

print(
    counts_file
)

print(
    "\nProbabilities:"
)

print(
    probability_file
)

print(
    "\nHRCA annotation completed."
)

print("=" * 70)
