# ------------------------------------------------------------------------------
# Annotating Seurat Clusters using Fisher's enrichment of BDGP in situ data
# TS O'Leary
# ------------------------------------------------------------------------------

# Load libraries
library(hdWGCNA)
library(tidyverse)
library(Seurat)
library(Signac)
library(cowplot)
library(patchwork)

# Output dir
out_dir <- "output/wgcna"

# Load data
dat <- readRDS(here::here("data/processed/seurat_object/10_dat_linked.rds"))
DefaultAssay(dat) <- "SCT"

# Set up the data -- exclude genes expressed in less than 5% of cells
dat <- SetupForWGCNA(
  dat,
  gene_select = "fraction",
  fraction = 0.05,
  wgcna_name = "wgcna"
)

# Construct metacells  in each group
dat <- MetacellsByGroups(
  seurat_obj = dat,
  group.by = c("cell_type", "acc_temp"), 
  reduction = "umap",
  k = 25,
  max_shared = 10,
  ident.group = "cell_type"
)

# Normalize meta cell expression matrix
dat <- NormalizeMetacells(
  dat
)

# Set dat expr
dat <- SetDatExpr(
  dat,
  group_name = "18°C",
  assay = "SCT",
  group.by = "acc_temp",
  slot = "data"
)

# Test different soft powers
dat <- TestSoftPowers(
  dat,
  networkType = "signed"
)

# Construct co-expression network:
dat <- ConstructNetwork(
  dat, 
  soft_power = 8,
  setDatExpr = FALSE,
  tom_outdir = out_dir,
  tom_name = "wgcna",
  overwrite_tom = TRUE
)

# Scale data first
dat <- ScaleData(
  dat,
  features = VariableFeatures(dat)
)

# Compute all MEs in the full single-cell dataset
dat <- ModuleEigengenes(
  dat,
  group.by.vars = "acc_temp"
)

# Find differential module eigengenes
DMEs <- FindDMEs(
  dat,
  barcodes1 = dat@meta.data |> filter(acc_temp == "18°C") |> rownames(),
  barcodes2 = dat@meta.data |> filter(acc_temp == "25°C") |> rownames(),
  test.use = "wilcox",
  wgcna_name = "wgcna"
)

# Save DMEs
saveRDS(DMEs, here::here(out_dir, "DMEs.rds"))

# Gene module assignment
modules <- GetModules(dat)

# Get genes from significant modules
brown <- modules |> 
  filter(module == "brown") |> 
  pull(gene_name)

# Get genes from significant modules
red <- modules |> 
  filter(module == "red") |> 
  pull(gene_name)

# Save modules
saveRDS(modules, here::here(out_dir, "modules.rds"))
write_csv(brown |> as_tibble(), here::here(out_dir, "brown.csv"))
write_csv(red |> as_tibble(), here::here(out_dir, "red.csv"))

# Save dat with hdWGCNA in it for plotting
saveRDS(dat, here::here("data/processed/seurat_object/13_dat_wgcna.rds"))

# ATAC -------------------------------------------------------------------------

# Output dir
out_dir <- "output/wgcna/atac"

# Load data
dat <- readRDS(here::here("data/processed/seurat_object/10_dat_linked.rds"))
DefaultAssay(dat) <- "peaks"

# Set up the data -- exclude genes expressed in less than 5% of cells
dat <- SetupForWGCNA(
  dat,
  gene_select = "fraction",
  fraction = 0.05,
  wgcna_name = "wgcna"
)

# Construct metacells  in each group
dat <- MetacellsByGroups(
  seurat_obj = dat,
  group.by = c("cell_type", "acc_temp"), 
  reduction = "umap",
  k = 25,
  max_shared = 10,
  ident.group = "cell_type"
)

# Normalize meta cell expression matrix
dat <- NormalizeMetacells(
  dat
)

# Set dat expr
dat <- SetDatExpr(
  dat,
  group_name = "18°C",
  assay = "peaks",
  group.by = "acc_temp",
  slot = "data"
)

# Test different soft powers
dat <- TestSoftPowers(
  dat,
  networkType = "signed"
)

# Construct co-expression network
dat <- ConstructNetwork(
  dat, 
  soft_power = 5,
  setDatExpr = FALSE,
  tom_outdir = out_dir,
  tom_name = "wgcna_atac",
  overwrite_tom = TRUE
)

# Scale data first
dat <- ScaleData(
  dat,
  features = VariableFeatures(dat)
)

# Compute all MEs in the full single-cell dataset
dat <- ModuleEigengenes(
  dat,
  assay = "peaks",
  group.by.vars = "acc_temp"
)

# Find differential module eigengenes
DMEs <- FindDMEs(
  dat,
  barcodes1 = dat@meta.data |> filter(acc_temp == "18°C") |> rownames(),
  barcodes2 = dat@meta.data |> filter(acc_temp == "25°C") |> rownames(),
  test.use = "wilcox",
  wgcna_name = "wgcna"
)

# Save DMEs
saveRDS(DMEs, here::here(out_dir, "DMEs.rds"))

# Gene module assignment
modules <- GetModules(dat)

# Get regions from significant modules
brown <- modules |> 
  filter(module == "brown") |> 
  pull(gene_name)

# Get regions from significant modules
turquoise <- modules |> 
  filter(module == "turquoise") |> 
  pull(gene_name)

# Save modules
saveRDS(modules, here::here(out_dir, "modules.rds"))
write_csv(brown |> as_tibble(), here::here(out_dir, "brown.csv"))
write_csv(turquoise |> as_tibble(), here::here(out_dir, "turquoise.csv"))

# Save dat with hdWGCNA in it for plotting
saveRDS(dat, here::here("data/processed/seurat_object/13_dat_wgcna_atac.rds"))

