# ------------------------------------------------------------------------------
# Dimension reduction and clustering
# TS O'Leary
# ------------------------------------------------------------------------------

# Load libraries
library(tidyverse)
library(Seurat)
library(Signac)
library(BSgenome.Dmelanogaster.UCSC.dm6)
library(SeuratDisk)

# Load data --------------------------------------------------------------------
dat <- readRDS(
  here::here("data/processed/seurat_object/06_dat_qc.rds")
)
# Save S phase genes in a vector
s_genes <- read_csv(
  here::here("data/raw/annot/dmel_cell-cycle_genes.csv")) |> 
  filter(phase == "S") |> 
  select(gene) |> 
  deframe()
# Save G2 or M phase genes in a vector
g2m_genes <- read_csv(
  here::here("data/raw/annot/dmel_cell-cycle_genes.csv")) |> 
  filter(phase == "G2/M") |> 
  select(gene) |> 
  deframe()

# Output dir
out_dir <- "data/processed/seurat_object"

# Joint clustering based on Signac tutorial ------------------------------------

# RNA data processing -----

# RNA data pre-processing before cell-cycle scoring
DefaultAssay(dat) <- "RNA"
dat <- SCTransform(dat)

# Score cell cycle
dat <- CellCycleScoring(
  dat,
  s.features = s_genes,
  g2m.features = g2m_genes,
  set.ident = TRUE
)

# Set assay back to RNA and re-run with S.Score and G2M.Score regressed out 
DefaultAssay(dat) <- "RNA"
dat <- SCTransform(
  dat, 
  vars.to.regress = c("S.Score", "G2M.Score"),
  variable.features.n = 3000
)
dat <- RunPCA(dat)

# ATAC data processing -----
DefaultAssay(dat) <- "peaks"
dat <- FindTopFeatures(dat, min.cutoff = 5)
dat <- RunTFIDF(dat)
dat <- RunSVD(dat)

# Build a joint neighbor graph using both assays -----
dat <- FindMultiModalNeighbors(
  object = dat,
  k.nn = 100,
  reduction.list = list("pca", "lsi"), 
  dims.list = list(1:50, 2:40),
  verbose = TRUE
)

# Build a joint UMAP
dat <- RunUMAP(
  object = dat,
  n.neighbors = 100L,
  nn.name = "weighted.nn",
  assay = "RNA",
  verbose = TRUE
)

# Build a joint t-SNE
dat <- RunTSNE(
  object = dat,
  nn.name = "weighted.nn",
  assay = "RNA",
  verbose = TRUE
)

# Find clusters on the wknn at resolutions ranging from 0.1 to 1
dat <- FindClusters(
  dat, 
  resolution = seq(from = 0.1, to = 1, by = 0.1),
  graph.name = "wknn", 
  algorithm = 3
)

# Save data
saveRDS(dat, here::here(out_dir, "07_dat_cluster.rds"))


################################################################################
# Creating RNA-only and ATAC-only dimension reductions -------------------------
# To see how and if the dimensional reduction look different 

# RNA-only dim-reduction -------------------------------------------------------
# Load data
dat <- readRDS(
  here::here("data/processed/seurat_object/07_dat_cluster.rds")
)

# Write over UMAP reduction with RNA-only UMAP
dat <- RunUMAP(
  dat, 
  assay = "SCT",
  dims = c(1:50)
)
# Write over t-SNE reduction with RNA-only t-SNE
dat <- RunTSNE(
  dat, 
  assay = "SCT",
  dims = c(1:50)
)

# Save dat with RNA only UMAP
saveRDS(dat,
        here::here(out_dir, "07_dat_rna_only_dim.rds")
)

# ATAC-only dim-reduction ------------------------------------------------------
# Load data
dat <- readRDS(
  here::here("data/processed/seurat_object/07_dat_cluster.rds")
)

DefaultAssay(dat) <- "peaks"
dat <- ScaleData(dat)
dat <- RunPCA(dat)
dat <- RunUMAP(dat, dims = c(1:50))
dat <- RunTSNE(dat, dims = c(1:50))

# Save dat with ATAC-only UMAP
saveRDS(dat,
        here::here(out_dir, "07_dat_atac_only_dim.rds")
)

# Cluster markers --------------------------------------------------------------

# Load libraries
library(tidyverse)
library(Seurat)
library(Signac)

# Load data
dat <- readRDS(here::here("data/processed/seurat_object/07_dat_cluster.rds"))

# Output dir
out_dir <- "output/markers"

# Set resolution for the seurat_clusters
dat@meta.data$seurat_clusters <- dat@meta.data$wknn_res.0.7
Idents(dat) <- "seurat_clusters"

# Find markers of each cluster that are conserved between the two conditions ---

# Create empty list to populate in the for loop
markers <- list()

# Loop through each cluster individually
for (i in levels(dat)) {
  markers[[i]] <- FindConservedMarkers(dat, 
                                       ident.1 = i,
                                       grouping.var = "acc_temp",
                                       only.pos = TRUE) |> 
    rownames_to_column("gene")
}

# Combine markers for all clusters
markers <- bind_rows(markers, 
                     .id = "cluster")

# Number of markers per cluster
markers |> 
  mutate(cluster = as.numeric(cluster)) |> 
  filter(max_pval < 0.05) |> 
  group_by(cluster) |> 
  tally()

# Save markers object and csv for additional manual annotation
saveRDS(markers, here::here(out_dir, "cluster_markers.rds"))
write_csv(markers, here::here(out_dir, "cluster_markers_all.csv"))
write_csv(
  markers |> 
    mutate(cluster = as.numeric(cluster)) |> 
    group_by(cluster) |>
    slice_min(max_pval, n = 10),
  here::here(out_dir, "cluster_markers_top_10.csv")
)

# Link peaks to genes ----------------------------------------------------------


# # Load data and set seq style
# dat <- readRDS(here::here("data/processed/seurat_object/09_dat_annot.rds"))
# seqlevelsStyle(BSgenome.Dmelanogaster.UCSC.dm6) <- "NCBI"
# seqlevelsStyle(dat@assays$peaks@annotation) <- "NCBI"

# Set assay
DefaultAssay(dat) <- "peaks"

# Compute the GC content for each peak
dat <- RegionStats(
  dat, 
  genome = BSgenome.Dmelanogaster.UCSC.dm6
)

# Link peaks to genes
dat <- LinkPeaks(
  object = dat,
  peak.assay = "peaks",
  expression.assay = "SCT"
)

# Save data
saveRDS(
  dat,
  here::here("data/processed/seurat_object/10_dat_linked.rds")
)
