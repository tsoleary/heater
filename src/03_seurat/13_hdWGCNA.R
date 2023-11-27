# ------------------------------------------------------------------------------
# High dimensional weighted gene co-expression analysis (hdWGCNA) of RNA and 
#   ATAC data 
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

# Module connectivity
# compute eigengene-based connectivity (kME):
dat <- ModuleConnectivity(
  dat,
  group.by = "acc_temp", 
  group_name = "25°C"
)

PlotKMEs(
  dat,
  ncol = 3,
  plot_widths = c(3, 1),
  text_size = 3,
  n_hubs = 8
)

# Save
ggsave(here::here(fig_dir, "kMEs.png"),
       height = 20,
       width = 25,
       units = "cm")

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

# Get genes from green
green <- modules |> 
  filter(module == "green") |> 
  pull(gene_name)

# Save modules
saveRDS(modules, here::here(out_dir, "modules.rds"))
write_csv(brown |> as_tibble(), here::here(out_dir, "brown.csv"))
write_csv(red |> as_tibble(), here::here(out_dir, "red.csv"))
write_csv(green |> as_tibble(), here::here(out_dir, "green.csv"))

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

# Module connectivity
# compute eigengene-based connectivity (kME):
dat <- ModuleConnectivity(
  dat,
  group.by = "acc_temp", 
  group_name = "25°C"
)

PlotKMEs(
  dat,
  ncol = 3,
  plot_widths = c(3, 1),
  text_size = 3,
  n_hubs = 8
)

# Gene module assignment
modules <- GetModules(dat)

# Get regions from significant modules
brown <- modules |> 
  filter(module == "brown") |> 
  pull(gene_name)

# Get regions from significant modules
blue <- modules |> 
  filter(module == "blue") |> 
  pull(gene_name)

# Save modules
saveRDS(modules, here::here(out_dir, "modules.rds"))
write_csv(brown |> as_tibble(), here::here(out_dir, "brown.csv"))
write_csv(blue |> as_tibble(), here::here(out_dir, "blue.csv"))

# Get genes linked to the blue peaks
blue_genes <- GetLinkedGenes(
  dat, 
  features = blue,
  min.abs.score = 0.1
)

# Get genes linked to the brown peaks
brown_genes <- GetLinkedGenes(
  dat, 
  features = brown,
  min.abs.score = 0.1
)

# Save the list of linked peaks
write_csv(brown_genes |> as_tibble(), here::here(out_dir, "brown_genes_0.1.csv"))
write_csv(blue_genes |> as_tibble(), here::here(out_dir, "blue_genes_0.1.csv"))

# Save dat with hdWGCNA in it for plotting
saveRDS(dat, here::here("data/processed/seurat_object/13_dat_wgcna_atac.rds"))

