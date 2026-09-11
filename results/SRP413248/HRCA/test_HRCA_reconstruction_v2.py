import torch
import anndata as ad
import pandas as pd
import pprint

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


# ============================================================
# 1. LOAD HRCA CHECKPOINT
# ============================================================

print("\n" + "=" * 70)
print("[1] Loading HRCA checkpoint")
print("=" * 70)

print("Model:")
print(MODEL_FILE)

checkpoint = torch.load(
    MODEL_FILE,
    map_location="cpu",
    weights_only=False,
)

print("\nCheckpoint loaded successfully.")

print("Checkpoint type:")
print(type(checkpoint))

print("\nCheckpoint keys:")
for key in checkpoint.keys():
    print("  ", key)


state_dict = checkpoint["model_state_dict"]

attr_dict = checkpoint["attr_dict"]

var_names = list(
    checkpoint["var_names"]
)

print("\nReference genes:")
print(len(var_names))

print("\nState_dict entries:")
print(len(state_dict))


# ============================================================
# 2. READ REFERENCE REGISTRY
# ============================================================

print("\n" + "=" * 70)
print("[2] Reading reference registry")
print("=" * 70)

registry = attr_dict["registry_"]

print("\nReference scvi version:")
print(registry.get("scvi_version"))

print("\nReference model:")
print(registry.get("model_name"))

setup_args = registry["setup_args"]

print("\nReference setup args:")
pprint.pp(setup_args)


# ============================================================
# 3. REFERENCE BATCHES
# ============================================================

print("\n" + "=" * 70)
print("[3] Extracting reference batch categories")
print("=" * 70)

batch_registry = (
    registry["field_registries"]["batch"]
)

reference_batches = list(
    batch_registry[
        "state_registry"
    ]["categorical_mapping"]
)

reference_n_batch = (
    batch_registry[
        "summary_stats"
    ]["n_batch"]
)

print("\nReference batch categories:")
print(
    "  Number:",
    len(reference_batches)
)

print(
    "  Registry n_batch:",
    reference_n_batch
)

assert len(reference_batches) == 229
assert reference_n_batch == 229


# ============================================================
# 4. REFERENCE CELL TYPES
# ============================================================

print("\n" + "=" * 70)
print("[4] Extracting reference cell-type categories")
print("=" * 70)

label_registry = (
    registry["field_registries"]["labels"]
)

reference_labels = list(
    label_registry[
        "state_registry"
    ]["categorical_mapping"]
)

reference_n_labels = (
    label_registry[
        "summary_stats"
    ]["n_labels"]
)

print("\nReference cell-type categories:")
print(
    "  Number:",
    len(reference_labels)
)

print(
    "  Registry n_labels:",
    reference_n_labels
)

print("\nReference cell types:")

for i, label in enumerate(reference_labels):
    print(
        f"  {i:3d}  {label}"
    )

assert len(reference_labels) == 124
assert reference_n_labels == 124


# ============================================================
# 5. LOAD ARMD QUERY
# ============================================================

print("\n" + "=" * 70)
print("[5] Loading ARMD query")
print("=" * 70)

print("Query:")
print(QUERY_FILE)

query = ad.read_h5ad(
    QUERY_FILE
)

print("\nQuery loaded.")

print("\nQuery shape:")
print(
    "  Cells:",
    query.n_obs
)

print(
    "  Genes:",
    query.n_vars
)

print("\nQuery obs columns:")

print(
    list(query.obs.columns)
)

assert query.n_vars == 10000


# ============================================================
# 6. CHECK GENES
# ============================================================

print("\n" + "=" * 70)
print("[6] Checking query/reference genes")
print("=" * 70)

query_genes = list(
    query.var_names
)

reference_gene_set = set(
    var_names
)

query_gene_set = set(
    query_genes
)

missing = (
    reference_gene_set
    - query_gene_set
)

extra = (
    query_gene_set
    - reference_gene_set
)

print("\nReference genes:")
print(
    len(reference_gene_set)
)

print("Query genes:")
print(
    len(query_gene_set)
)

print("Missing:")
print(
    len(missing)
)

print("Extra:")
print(
    len(extra)
)

if len(missing) > 0:

    print("\nMissing genes:")

    for gene in sorted(missing):
        print(
            "  ",
            gene
        )

if len(extra) > 0:

    print("\nExtra genes:")

    for gene in sorted(extra):
        print(
            "  ",
            gene
        )

if (
    len(missing) > 0
    or len(extra) > 0
):

    raise ValueError(
        "\nQuery genes do not exactly "
        "match HRCA reference genes."
    )

print(
    "\nGenes match exactly."
)


# ------------------------------------------------------------
# Reorder genes exactly like HRCA
# ------------------------------------------------------------

query = query[
    :,
    var_names
].copy()

print(
    "Query reordered:",
    query.shape
)


# ============================================================
# 7. QUERY METADATA
# ============================================================

print("\n" + "=" * 70)
print("[7] Preparing query metadata")
print("=" * 70)


# ------------------------------------------------------------
# Cell type
# ------------------------------------------------------------

query.obs["celltype"] = "Unknown"


# ------------------------------------------------------------
# Batch/sample
# ------------------------------------------------------------

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


print("\nQuery batches:")

print(
    "Number:",
    query.obs[
        "sampleid"
    ].nunique()
)

print("\nCells per batch:")

print(
    query.obs[
        "sampleid"
    ]
    .value_counts()
    .sort_index()
)


# ============================================================
# 8. CREATE DUMMY REFERENCE DATA
# ============================================================

print("\n" + "=" * 70)
print("[8] Preparing reference-compatible dummy AnnData")
print("=" * 70)

print(
    "\nThe dummy AnnData is only used "
    "to reconstruct the original HRCA "
    "model architecture."
)

print(
    "It is NOT used for annotation."
)


n_reference_batches = (
    len(reference_batches)
)


# ------------------------------------------------------------
# Batch categories
# ------------------------------------------------------------

dummy_batches = pd.Categorical(
    reference_batches,
    categories=reference_batches,
)


# ------------------------------------------------------------
# Label categories
#
# IMPORTANT:
# Keep all 124 HRCA categories.
# Actual dummy cells are Unknown.
# ------------------------------------------------------------

dummy_labels = pd.Categorical(
    ["Unknown"] *
    n_reference_batches,

    categories=reference_labels,
)


# ------------------------------------------------------------
# Create dummy object
# ------------------------------------------------------------

dummy = ad.AnnData(

    X=query.X[
        :n_reference_batches,
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
            for i in range(
                n_reference_batches
            )
        ],
    ),

    var=query.var.copy(),
)


# ------------------------------------------------------------
# Copy counts layer
# ------------------------------------------------------------

if "counts" in query.layers:

    dummy.layers["counts"] = (
        query.layers[
            "counts"
        ][
            :n_reference_batches,
            :
        ].copy()
    )


print("\nDummy shape:")
print(
    dummy.shape
)

print(
    "\nDummy batch categories:",
    dummy.obs[
        "sampleid"
    ].cat.categories.size
)

print(
    "Dummy label categories:",
    dummy.obs[
        "celltype"
    ].cat.categories.size
)


# ============================================================
# 9. SCANVI SETUP
# ============================================================

print("\n" + "=" * 70)
print("[9] Running SCANVI.setup_anndata()")
print("=" * 70)

print(
    "\nlabels_key:",
    setup_args["labels_key"]
)

print(
    "unlabeled_category:",
    setup_args[
        "unlabeled_category"
    ]
)

print(
    "batch_key:",
    setup_args["batch_key"]
)


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


print(
    "\nsetup_anndata completed."
)


# ============================================================
# 10. READ INITIALIZATION PARAMETERS
# ============================================================

print("\n" + "=" * 70)
print("[10] Reading SCANVI initialization parameters")
print("=" * 70)

init_params = (
    attr_dict["init_params_"]
)

print("\nRaw init_params:")

pprint.pp(
    init_params
)


# ------------------------------------------------------------
# Non-kwargs
# ------------------------------------------------------------

non_kwargs = dict(
    init_params[
        "non_kwargs"
    ]
)


# ------------------------------------------------------------
# IMPORTANT:
#
# The checkpoint stores:
#
# kwargs:
#   model_kwargs:
#       latent_distribution
#       use_batch_norm
#       use_layer_norm
#       deeply_inject_covariates
#       encode_covariates
#
# We need the CONTENTS of model_kwargs.
#
# We must NOT pass:
#
# model_kwargs={...}
#
# directly to SCANVI.
# ------------------------------------------------------------

model_kwargs = dict(
    init_params[
        "kwargs"
    ][
        "model_kwargs"
    ]
)


print("\nFinal non_kwargs:")

pprint.pp(
    non_kwargs
)


print("\nFinal model kwargs:")

pprint.pp(
    model_kwargs
)


# ============================================================
# 11. INITIALIZE REFERENCE ARCHITECTURE
# ============================================================

print("\n" + "=" * 70)
print("[11] Initializing SCANVI architecture")
print("=" * 70)


model = SCANVI(

    dummy,

    **non_kwargs,

    **model_kwargs,
)


print(
    "\nSCANVI architecture "
    "initialized successfully."
)


# ============================================================
# 12. CHECK SUMMARY STATISTICS
# ============================================================

print("\n" + "=" * 70)
print("[12] Checking reconstructed model")
print("=" * 70)


model_n_batch = (
    model.summary_stats[
        "n_batch"
    ]
)

model_n_labels = (
    model.summary_stats[
        "n_labels"
    ]
)


print(
    "\nModel n_batch:",
    model_n_batch
)

print(
    "Expected n_batch:",
    reference_n_batch
)

print(
    "\nModel n_labels:",
    model_n_labels
)

print(
    "Expected n_labels:",
    reference_n_labels
)


if (
    model_n_batch
    != reference_n_batch
):

    raise RuntimeError(

        "\nBATCH ARCHITECTURE MISMATCH\n"

        f"Model: {model_n_batch}\n"

        f"Reference: {reference_n_batch}"
    )


if (
    model_n_labels
    != reference_n_labels
):

    raise RuntimeError(

        "\nLABEL ARCHITECTURE MISMATCH\n"

        f"Model: {model_n_labels}\n"

        f"Reference: {reference_n_labels}"
    )


print(
    "\nBatch and label architecture "
    "match the HRCA registry."
)


# ============================================================
# 13. COMPARE STATE DICT
# ============================================================

print("\n" + "=" * 70)
print("[13] Comparing checkpoint/model architecture")
print("=" * 70)


model_state = (
    model.module.state_dict()
)


mismatches = []

missing_in_model = []

missing_in_checkpoint = []

shape_mismatches = []


# ------------------------------------------------------------
# Check checkpoint keys
# ------------------------------------------------------------

for key, checkpoint_tensor in (
    state_dict.items()
):

    if key not in model_state:

        missing_in_model.append(
            key
        )

        mismatches.append(

            (
                key,

                "MISSING_IN_MODEL",

                tuple(
                    checkpoint_tensor.shape
                ),

                None,
            )
        )

        continue


    model_tensor = (
        model_state[key]
    )


    checkpoint_shape = (
        tuple(
            checkpoint_tensor.shape
        )
    )

    model_shape = (
        tuple(
            model_tensor.shape
        )
    )


    if (
        checkpoint_shape
        != model_shape
    ):

        shape_mismatches.append(

            (
                key,

                checkpoint_shape,

                model_shape,
            )
        )

        mismatches.append(

            (
                key,

                "SHAPE_MISMATCH",

                checkpoint_shape,

                model_shape,
            )
        )


# ------------------------------------------------------------
# Check model keys
# ------------------------------------------------------------

for key in model_state:

    if key not in state_dict:

        missing_in_checkpoint.append(
            key
        )

        mismatches.append(

            (
                key,

                "MISSING_IN_CHECKPOINT",

                None,

                tuple(
                    model_state[
                        key
                    ].shape
                ),
            )
        )


# ============================================================
# 14. PRINT COMPARISON
# ============================================================

print("\nCheckpoint parameters:")
print(
    len(state_dict)
)

print("\nModel parameters:")
print(
    len(model_state)
)

print("\nMissing in model:")
print(
    len(missing_in_model)
)

print("\nMissing in checkpoint:")
print(
    len(missing_in_checkpoint)
)

print("\nShape mismatches:")
print(
    len(shape_mismatches)
)


# ------------------------------------------------------------
# Shape mismatches
# ------------------------------------------------------------

if len(shape_mismatches) > 0:

    print(
        "\n" + "-" * 70
    )

    print(
        "SHAPE MISMATCHES"
    )

    print(
        "-" * 70
    )

    for (
        key,
        checkpoint_shape,
        model_shape,
    ) in shape_mismatches:

        print(
            "\nKEY:",
            key
        )

        print(
            "  checkpoint:",
            checkpoint_shape
        )

        print(
            "  model:",
            model_shape
        )


# ------------------------------------------------------------
# Missing in model
# ------------------------------------------------------------

if len(missing_in_model) > 0:

    print(
        "\n" + "-" * 70
    )

    print(
        "MISSING IN MODEL"
    )

    print(
        "-" * 70
    )

    for key in missing_in_model:

        print(
            "  ",
            key
        )


# ------------------------------------------------------------
# Missing in checkpoint
# ------------------------------------------------------------

if len(missing_in_checkpoint) > 0:

    print(
        "\n" + "-" * 70
    )

    print(
        "MISSING IN CHECKPOINT"
    )

    print(
        "-" * 70
    )

    for key in missing_in_checkpoint:

        print(
            "  ",
            key
        )


# ============================================================
# 15. STOP IF ARCHITECTURE DOES NOT MATCH
# ============================================================

if len(mismatches) > 0:

    raise RuntimeError(

        "\n\nARCHITECTURE MISMATCH\n"

        f"Total mismatches: "
        f"{len(mismatches)}\n\n"

        "The reconstructed model does not "
        "yet exactly match the HRCA checkpoint."
    )


print(
    "\nALL CHECKPOINT/MODEL PARAMETERS MATCH."
)


# ============================================================
# 16. LOAD REFERENCE WEIGHTS
# ============================================================

print("\n" + "=" * 70)
print("[16] Loading HRCA reference state_dict")
print("=" * 70)


result = (
    model.module.load_state_dict(
        state_dict,
        strict=True,
    )
)


print(
    "\nReference state_dict "
    "loaded successfully!"
)

print(
    result
)


# ============================================================
# 17. PREPARE ARMD QUERY
# ============================================================

print("\n" + "=" * 70)
print("[17] Preparing ARMD query")
print("=" * 70)


print(
    "\nOriginal query:"
)

print(
    "  Cells:",
    query.n_obs
)

print(
    "  Genes:",
    query.n_vars
)


query_prepared = (
    SCANVI.prepare_query_anndata(

        query,

        reference_model=model,

        inplace=False,
    )
)


print(
    "\nPrepared query:"
)

print(
    "  Cells:",
    query_prepared.n_obs
)

print(
    "  Genes:",
    query_prepared.n_vars
)


# ============================================================
# 18. LOAD QUERY INTO REFERENCE
# ============================================================

print("\n" + "=" * 70)
print("[18] Loading ARMD query into HRCA reference")
print("=" * 70)


print(
    "\nAdapting ARMD query to HRCA..."
)


query_model = (
    SCANVI.load_query_data(

        query_prepared,

        model,

        freeze_dropout=True,
    )
)


print(
    "\nARMD query successfully "
    "loaded into HRCA reference."
)


# ============================================================
# 19. SUCCESS
# ============================================================

print("\n")

print("=" * 70)
print("SUCCESS")
print("=" * 70)


print(
    "\nQuery cells:",
    query_model.adata.n_obs
)

print(
    "Query genes:",
    query_model.adata.n_vars
)

print(
    "Model:",
    type(query_model)
)


print(
    "\nHRCA reference reconstruction "
    "completed successfully."
)

print(
    "\nARMD query successfully adapted "
    "to the HRCA reference."
)


print(
    "\nNext step:"
)

print(
    "    query_model.predict()"
)

print(
    "\nThis will generate HRCA "
    "cell-type predictions."
)

print(
    "\n" + "=" * 70
)
