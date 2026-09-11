import copy
import torch
import anndata as ad
import pandas as pd

from scvi.model import SCANVI

# ============================================================
# PATHS
# ============================================================

MODEL_FILE = "/storage/lemus_g/roldan/ARMD/HRCA_reference/model.pt"
QUERY_FILE = "SRP413248_HRCA_query_symbols_10000.h5ad"

# ============================================================
# SETTINGS
# ============================================================

N_CELLS = 500
SEED = 12345

print("=" * 70)
print("HRCA SCANVI RECONSTRUCTION TEST")
print("=" * 70)

# ============================================================
# LOAD CHECKPOINT
# ============================================================

print("\n[1] Loading HRCA checkpoint...")

checkpoint = torch.load(
    MODEL_FILE,
    map_location="cpu",
    weights_only=False
)

print("Checkpoint keys:")
print(list(checkpoint.keys()))

attr_dict = copy.deepcopy(checkpoint["attr_dict"])
var_names = pd.Index(checkpoint["var_names"])
state_dict = checkpoint["model_state_dict"]

print("\nReference genes:", len(var_names))
print("State dict entries:", len(state_dict))

# ============================================================
# LOAD QUERY
# ============================================================

print("\n[2] Loading query...")

adata_full = ad.read_h5ad(QUERY_FILE)

print("Full query:", adata_full.shape)

# ------------------------------------------------------------
# Random 500 cells
# ------------------------------------------------------------

adata = adata_full[
    adata_full.obs.sample(
        n=N_CELLS,
        random_state=SEED
    ).index
].copy()

print("Test query:", adata.shape)

# ============================================================
# PREPARE OBS
# ============================================================

print("\n[3] Preparing query metadata...")

# SCANVI expects these names from the reference registry
if "celltype" in adata.obs.columns:
    adata.obs["celltype_original"] = adata.obs["celltype"]

adata.obs["celltype"] = "Unknown"

if "sampleid" not in adata.obs.columns:
    if "sample" in adata.obs.columns:
        adata.obs["sampleid"] = adata.obs["sample"].astype(str)
    elif "GSM" in adata.obs.columns:
        adata.obs["sampleid"] = adata.obs["GSM"].astype(str)
    else:
        adata.obs["sampleid"] = "ARMD_query"

print("celltype:")
print(adata.obs["celltype"].value_counts())

print("\nsampleid:")
print(adata.obs["sampleid"].value_counts().head())

# ============================================================
# REFERENCE REGISTRY
# ============================================================

print("\n[4] Reading reference registry...")

registry = copy.deepcopy(attr_dict["registry_"])

print("Model name:", registry.get("model_name"))
print("scvi version:", registry.get("scvi_version"))
print("Setup method:", registry.get("setup_method_name"))

print("\nSetup args:")
print(registry["setup_args"])

# ============================================================
# SETUP ANNDATA USING REFERENCE REGISTRY
# ============================================================

print("\n[5] Running SCANVI.setup_anndata()...")

setup_method_name = registry["setup_method_name"]
setup_method = getattr(SCANVI, setup_method_name)

setup_args = copy.deepcopy(registry["setup_args"])

setup_method(
    adata,
    source_registry=registry,
    extend_categories=True,
    allow_missing_labels=True,
    **setup_args,
)

print("setup_anndata: SUCCESS")

# ============================================================
# REMOVE REGISTRY FROM INIT ATTRIBUTES
# ============================================================

init_attr_dict = copy.deepcopy(attr_dict)

# _initialize_model normally removes registry_ before constructing model
init_attr_dict.pop("registry_", None)

# ============================================================
# INITIALIZE SCANVI
# ============================================================

print("\n[6] Initializing SCANVI architecture...")

init_params = init_attr_dict["init_params_"]

print("init_params:")
print(init_params)

# Extract kwargs used during original initialization
kwargs = copy.deepcopy(init_params["kwargs"])
non_kwargs = copy.deepcopy(init_params["non_kwargs"])

# model_kwargs are stored inside kwargs
model_kwargs = kwargs.get("model_kwargs", {})

model = SCANVI(
    adata,
    **non_kwargs,
    **model_kwargs
)

print("\nSCANVI architecture initialized.")
print("Model:", type(model))

# ============================================================
# LOAD REFERENCE WEIGHTS
# ============================================================

print("\n[7] Loading reference state_dict...")

model.module.load_state_dict(
    state_dict,
    strict=True
)

print("STATE DICT: SUCCESS")

model.module.eval()
model.is_trained_ = True

# ============================================================
# REFERENCE MODEL SUMMARY
# ============================================================

print("\n[8] Reference model summary")

print("n_vars:", model.summary_stats.n_vars)
print("n_batch:", model.summary_stats.n_batch)
print("n_labels:", model.summary_stats.n_labels)

# ============================================================
# PREPARE QUERY AGAINST REFERENCE
# ============================================================

print("\n[9] Testing prepare_query_anndata()...")

query_prepared = SCANVI.prepare_query_anndata(
    adata=adata,
    reference_model=model,
    inplace=False,
)

print("prepare_query_anndata: SUCCESS")
print("Prepared query shape:", query_prepared.shape)

# ============================================================
# LOAD QUERY DATA / SCARCHES SURGERY
# ============================================================

print("\n[10] Testing load_query_data()...")

query_model = SCANVI.load_query_data(
    query_prepared,
    reference_model=model,
    accelerator="cpu",
    device=1,
    freeze_dropout=True,
)

print("\nload_query_data: SUCCESS")
print("Query model:", type(query_model))

print("\n==============================================")
print("ALL RECONSTRUCTION TESTS PASSED")
print("==============================================")
