# ------------------------------------------------------------------------------
# hdWGCNA on cell type subsets of the data
# High dimensional weighted gene co-expression analysis (hdWGCNA) of RNA and 
#   ATAC data 
# TS O'Leary
# ------------------------------------------------------------------------------

# Load libraries
library(tidyverse)
library(Seurat)
library(Signac)
library(hdWGCNA)
library(cowplot)
library(patchwork)

# Define base output directory
out_dir <- "output/wgcna/cell_type"

# Load data
dat_all <- readRDS(here::here("data/processed/seurat_object/10_dat_linked.rds"))

# Define cell types of interest
cell_types <- c(
  "ectoderm prim.", "endoderm prim.",
  "mesoderm prim.", "ventral nerve cord prim.",
  "peripheral nervous system prim.", "tracheal prim.", 
  "foregut & hindgut prim.", "amnioserosa", "germ cell"
)


# Part 1: Pre-network analysis up to soft power testing for all cell types -----

# Create lists to store pre-network processed objects and soft power results
dat_list <- list()
softpower_results_list <- list()

# Loop through each cell type
for (ct in cell_types) {
  message("Processing cell type: ", ct)
  
  # Subset the data for the current cell type
  dat_sub <- subset(dat_all, subset = cell_type == ct)
  
  # Create an output directory for the current cell type
  ct_out_dir <- file.path(out_dir, gsub(" ", "_", ct))
  if (!dir.exists(ct_out_dir)) {
    dir.create(ct_out_dir, recursive = TRUE)
  }
  
  # Set the default assay and prepare data for hdWGCNA
  DefaultAssay(dat_sub) <- "SCT"
  dat <- SetupForWGCNA(
    dat_sub,
    gene_select = "fraction",
    fraction = 0.05,
    wgcna_name = "wgcna"
  )
  
  # Construct metacells using the "acc_temp" grouping and "umap" reduction
  dat <- MetacellsByGroups(
    seurat_obj = dat,
    group.by = c("acc_temp", "cell_type"),
    reduction = "umap",
    k = 25,
    max_shared = 10,
    ident.group = "cell_type"
  )
  
  # Normalize the metacell expression matrix
  dat <- NormalizeMetacells(dat)
  
  # Create a new assay using just the counts
  new_SCT <- CreateAssayObject(counts = GetAssayData(dat, assay = "SCT", slot = "counts"))
  
  # Then set the 'data' slot manually
  new_SCT@data <- GetAssayData(dat, assay = "SCT", slot = "data")
  
  # Add the new assay to your object and set it as default
  dat[["SCT"]] <- new_SCT
  DefaultAssay(dat) <- "SCT"
  
  # Set expression data for a chosen group (e.g., "18°C")
  dat <- SetDatExpr(
    dat,
    group_name = "18°C",
    assay = "SCT",
    group.by = "acc_temp",
    slot = "data"
  )
  
  # Test different soft powers; results are stored in the object
  dat <- TestSoftPowers(dat, networkType = "signed")
  
  # Store the processed object and its soft power results
  dat_list[[ct]] <- dat
  softpower_results_list[[ct]] <- dat@wgcna_params$softpower_results
}

# Optionally, inspect the soft power results for each cell type:
print(softpower_results_list)

# Part 2: Continue analysis with soft power for each cell type -----------------
selected_softpowers <- list(
  "ectoderm prim." = 8,
  "endoderm prim." = 7,
  "mesoderm prim." = 9,
  "ventral nerve cord prim." = 8,
  "peripheral nervous system prim." = 7,
  "tracheal prim." = 8,
  "foregut & hindgut prim." = 7,
  "amnioserosa" = 9,
  "germ cell" = 8
)

# Loop through each cell type to complete the analysis with its specified soft power
for (ct in cell_types) {
  message("Completing analysis for cell type: ", ct)
  
  # Retrieve the pre-processed data for this cell type
  dat <- dat_list[[ct]]
  
  # Create (or verify) output directory for this cell type
  ct_out_dir <- file.path(out_dir, gsub(" ", "_", ct))
  if (!dir.exists(ct_out_dir)) {
    dir.create(ct_out_dir, recursive = TRUE)
  }
  
  # Retrieve the chosen soft power for this cell type
  sp <- selected_softpowers[[ct]]
  message("Using soft power ", sp, " for ", ct)
  
  # Construct the co-expression network using the chosen soft power
  dat <- ConstructNetwork(
    dat, 
    soft_power = sp,
    setDatExpr = FALSE,
    tom_outdir = ct_out_dir,
    tom_name = "wgcna",
    overwrite_tom = TRUE
  )
  
  # Scale data using variable features
  dat <- ScaleData(dat, features = VariableFeatures(dat))
  
  # Compute module eigengenes
  dat <- ModuleEigengenes(dat, group.by.vars = "acc_temp")
  
  # Identify differential module eigengenes between two conditions (adjust conditions as needed)
  DMEs <- FindDMEs(
    dat,
    barcodes1 = dat@meta.data |> filter(acc_temp == "18°C") |> rownames(),
    barcodes2 = dat@meta.data |> filter(acc_temp == "25°C") |> rownames(),
    test.use = "wilcox",
    wgcna_name = "wgcna"
  )
  
  # Compute module connectivity (kME) for a chosen condition (e.g., "25°C")
  dat <- ModuleConnectivity(
    dat,
    group.by = "acc_temp", 
    group_name = "25°C"
  )
  
  # Get the gene module assignments
  modules <- GetModules(dat)
  
  # Save outputs for this cell type
  saveRDS(DMEs, file.path(ct_out_dir, "DMEs.rds"))
  saveRDS(modules, file.path(ct_out_dir, "modules.rds"))
  saveRDS(dat, file.path("data/processed/seurat_object", paste0("dat_wgcna_", gsub(" ", "_", ct), ".rds")))
}
