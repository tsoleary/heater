# ------------------------------------------------------------------------------
# Initial data processing
# Basic quality control metrics and filtering out non-cell barcodes from the 10X 
# Cell Ranger ARC called cells.
# Dimension reduction and clustering to be used for DoubletFinder
# AMULET to detect doublets in scATAC-Seq > 2 overlapping reads
# DoubletFinder for detecting multiplets in the scRNA libraries
# Remove multiplets


# TS O'Leary
# ------------------------------------------------------------------------------

# Load libraries
library(tidyverse)
library(Seurat)
library(Signac)
library(clustree)
library(scDblFinder)
library(DoubletFinder)

# Load the RNA-Seq and ATAC-Seq data -------------------------------------------
data_dir <- here::here("data/processed/seq/all/outs/")
raw_feature_dir <- paste0(data_dir, "raw_feature_bc_matrix")
frag_path <- paste0(data_dir, "atac_fragments.tsv.gz")
counts <- Read10X(raw_feature_dir)

# Output dir
out_dir <- "data/processed/seurat_object"

# Get gene annotations for dm6 -------------------------------------------------
BDGP6.32 <- AnnotationHub::query(AnnotationHub::AnnotationHub(), 
                                 c("EnsDb", "Drosophila melanogaster", "109"))[[1]]
annotation <- GetGRangesFromEnsDb(BDGP6.32)

# Meta data --------------------------------------------------------------------
# Create meta data for all cells based on the GEM well suffix. From 10X Genomics
# "This number, which indicates which GEM well the barcode sequence came from, 
# is called the GEM well suffix. The numbering of the GEM wells will reflect the 
# order that the GEM wells were provided in the Aggregation CSV."
meta_data <- tibble(cellnames = colnames(counts[[1]])) |> 
  mutate(orig.ident = str_extract(cellnames, "[1-4]")) |> 
  full_join(tibble::tribble(
    ~orig.ident, ~sample_name, ~acc_temp,
    "1", "18C_Rep1", "18°C",
    "2", "18C_Rep2", "18°C",
    "3", "25C_Rep1", "25°C",
    "4", "25C_Rep2", "25°C"
  ), 
  by = "orig.ident") |> 
  column_to_rownames("cellnames")

# Create a Seurat object containing the RNA data
dat <- CreateSeuratObject(
  project = "heater",
  counts = counts$`Gene Expression`,
  assay = "RNA",
  names.delim = "-",
  meta.data = meta_data
)

# Create ATAC assay and add it to the object
dat[["ATAC"]] <- CreateChromatinAssay(
  counts = counts$Peaks,
  sep = c(":", "-"),
  fragments = frag_path,
  annotation = annotation,
  validate.fragments = FALSE
)

# Save raw Seurat object with all barcodes included for plotting
saveRDS(dat, here::here(out_dir, "00_dat_raw.rds"))

# Create Seurat object from the filtered matrix from 10X Cell Ranger ARC -------

# Load filtered data
filtered_feature_dir <- paste0(data_dir, "filtered_feature_bc_matrix")
counts <- Read10X(filtered_feature_dir)

# Meta data for counts from filtered counts
meta_data <- tibble(cellnames = colnames(counts[[1]])) |>
  mutate(orig.ident = str_extract(cellnames, "[1-4]")) |>
  full_join(tibble::tribble(
    ~orig.ident, ~sample_name, ~acc_temp,
    "1", "18C_Rep1", "18°C",
    "2", "18C_Rep2", "18°C",
    "3", "25C_Rep1", "25°C",
    "4", "25C_Rep2", "25°C"
  ), 
  by = "orig.ident") |>
  column_to_rownames("cellnames")

# Create a Seurat object containing the RNA data
dat <- CreateSeuratObject(
  project = "heater",
  counts = counts$`Gene Expression`,
  assay = "RNA",
  names.delim = "-",
  meta.data = meta_data
)

# Create ATAC assay and add it to the object
dat[["ATAC"]] <- CreateChromatinAssay(
  counts = counts$Peaks,
  sep = c(":", "-"),
  fragments = frag_path,
  annotation = annotation
)

# Save 10X Cell Ranger ARC filtered Seurat object
saveRDS(dat, here::here(out_dir, "00_dat_10x_cells.rds"))

# Calculate Nucleosome Signal and TSS Enrichment for ATAC QC metrics
DefaultAssay(dat) <- "ATAC"
dat <- NucleosomeSignal(dat)
dat <- TSSEnrichment(dat)

# Calculate the percentage of RNA reads that map to mitochondrial genes 
DefaultAssay(dat) <- "RNA"
dat[["percent.mt"]] <- PercentageFeatureSet(
  dat,
  pattern = "^mt:"
)
# Calculate the percentage of RNA reads that map to ribosomal genes
dat[["percent.ribo"]] <- PercentageFeatureSet(
  dat, 
  pattern = "^Rp[S|L]"
)

# Save object with Nucleosome, TSS, and mitochondria meta.data added 
saveRDS(
  dat, 
  here::here(out_dir, "01_dat_10x_cells.rds")
)

# Filter out barcodes with low counts or high mitochondrial content ------------

# Define filtering thresholds
low_ATAC <- 800
low_RNA <- 200
max_mt <- 5
max_rb <- 25

# Subset 10x data based on thresholds
dat <- dat |>
  subset(
    nCount_RNA >= low_RNA &
      nCount_ATAC >= low_ATAC & 
      percent.mt < max_mt &
      percent.ribo < max_rb
  )

# Save object
saveRDS(
  dat, 
  here::here(out_dir, "00_dat_10x_filtered.rds")
)

# Create a filtered fragment file to be used in the Amulet doublet finder ------

# Write the file with fragments from only the filtered barcodes
FilterCells(
  frag_path,
  Cells(dat),
  outfile = here::here(data_dir, "atac_fragments_cells.tsv.gz")
)

# ------------------------------------------------------------------------------
# Dimension reduction and clustering to be used for DoubletFinder
# TS O'Leary
# ------------------------------------------------------------------------------

# Set RNA to the default assay set for the following set of operations 
DefaultAssay(dat) <- "RNA"

# Log-Normalize data
dat <- NormalizeData(dat)

# Find variable features
dat <- FindVariableFeatures(dat)

# Scale and center data
dat <- ScaleData(dat)

# Run principle component analysis 
dat <- RunPCA(dat)

# Run UMAP
dat <- RunUMAP(dat, dims = 1:30)

# Construct nearest-neighbor graph
dat <- FindNeighbors(dat, dims = 1:30)

# Finding the correct resolution of cluster -----
# Run on a bunch of different resolutions
dat <- FindClusters(dat, resolution = seq(0.01, 0.1, 0.01))
# Look at the cluster tree
clustree::clustree(dat)

# Remove the RNA_snn columns added
dat@meta.data <- dat@meta.data |>
  dplyr::select(!contains("RNA_snn"))

# Set the final clusters
dat <- FindClusters(dat, resolution = 0.03)

# Save data
saveRDS(
  dat, 
  here::here(out_dir, "02_dat_initial_cluster.rds")
)

# ------------------------------------------------------------------------------
# AMULET to detect doublets in scATAC-Seq > 2 overlapping reads
# TS O'Leary
# ------------------------------------------------------------------------------

# Regions to exclude: sex chromosomes, mitochondrial genome, non-chromosomal 
#   fragments.
exclude_regions <- GenomicRanges::GRanges(
  c("X", 
    "Y", 
    "mitochondrion_genome", 
    "Unmapped_Scaffold_8_D1580_D1567",
    "211000022280328", 
    "211000022278760"),
  IRanges::IRanges(1L, width = 10^9)
)

# Run amulet on fragments excluding regions
doublet_amulet <- amulet(frag_path, regionsToExclude = exclude_regions)

# Save the list of doublets from amulet
saveRDS(doublet_amulet, here::here("data/processed/qc", "doublet_amulet.rds"))

# ------------------------------------------------------------------------------
# DoubletFinder for detecting multiplets in the scRNA libraries -----------
# TS O'Leary
# ------------------------------------------------------------------------------

# DoubletFinder for RNA doublets -----------------------------------------------

# Define expected doublet rate from 10X estimates
doublet_rate <- 0.008/1000
target_total_nuclei <- 4000
exp_double_rate <- doublet_rate*target_total_nuclei

# Each sample must be done separately per the instructions of DoubletFinder

# 18C_Rep1 ---------------------------------------------------------------------

# Load data
dat <- readRDS(
  here::here("data/processed/seurat_object/01_dat_10x_filtered.rds")) |>
  base::subset(sample_name == "18C_Rep1")

# Pre-process Seurat object ----------------------------------------------------
DefaultAssay(dat) <- "RNA"
dat <- NormalizeData(dat)
dat <- FindVariableFeatures(dat, selection.method = "vst", nfeatures = 2000)
dat <- ScaleData(dat)
dat <- RunPCA(dat)
dat <- RunUMAP(dat, dims = 1:10)
dat <- FindNeighbors(dat, dims = 1:10)
dat <- FindClusters(dat, resolution = 0.1)

# pK Identification ------------------------------------------------------------
sweep.res.list_dat <- paramSweep_v3(dat, PCs = 1:10, sct = FALSE)
sweep.stats_dat <- summarizeSweep(sweep.res.list_dat, GT = FALSE)
bcmvn_dat <- find.pK(sweep.stats_dat)

# Homotypic Doublet Proportion Estimate ----------------------------------------
homotypic.prop <- modelHomotypic(dat@meta.data$seurat_clusters)
nExp_poi <- round(exp_double_rate * nrow(dat@meta.data))
nExp_poi.adj <- round(nExp_poi*(1 - homotypic.prop))

# Run DoubletFinder with varying classification stringencies -------------------
dat <- doubletFinder_v3(dat, 
                        PCs = 1:10,
                        pK = 0.03,
                        nExp = nExp_poi.adj,
                        reuse.pANN = FALSE, 
                        sct = FALSE)

# Save the classified doublets
doublets_18C_Rep1 <- dat@meta.data |>
  base::subset(DF.classifications_0.25_0.03_38 == "Doublet") |>
  base::rownames()

# 18C_Rep2 ---------------------------------------------------------------------

# Load data
dat <- readRDS(
  here::here("data/processed/seurat_object/01_dat_10x_filtered.rds")) |>
  base::subset(sample_name == "18C_Rep2")

# Pre-process Seurat object ----------------------------------------------------
DefaultAssay(dat) <- "RNA"
dat <- NormalizeData(dat)
dat <- FindVariableFeatures(dat, selection.method = "vst", nfeatures = 2000)
dat <- ScaleData(dat)
dat <- RunPCA(dat)
dat <- RunUMAP(dat, dims = 1:10)
dat <- FindNeighbors(dat, dims = 1:10)
dat <- FindClusters(dat, resolution = 0.1)

# pK Identification ------------------------------------------------------------
sweep.res.list_dat <- paramSweep_v3(dat, PCs = 1:10, sct = FALSE)
sweep.stats_dat <- summarizeSweep(sweep.res.list_dat, GT = FALSE)
bcmvn_dat <- find.pK(sweep.stats_dat)

# Homotypic Doublet Proportion Estimate ----------------------------------------
homotypic.prop <- modelHomotypic(dat@meta.data$seurat_clusters)
nExp_poi <- round(exp_double_rate * nrow(dat@meta.data))
nExp_poi.adj <- round(nExp_poi*(1 - homotypic.prop))

# Run DoubletFinder with varying classification stringencies -------------------
dat <- doubletFinder_v3(dat, 
                        PCs = 1:10,
                        pK = 0.005,
                        nExp = nExp_poi.adj,
                        reuse.pANN = FALSE, 
                        sct = FALSE)

# Save the classified doublets
doublets_18C_Rep2 <- dat@meta.data |>
  base::subset(DF.classifications_0.25_0.005_107 == "Doublet") |>
  base::rownames()


# 25C_Rep1 ---------------------------------------------------------------------

# Load data
dat <- readRDS(
  here::here("data/processed/seurat_object/01_dat_10x_filtered.rds")) |>
  base::subset(sample_name == "25C_Rep1")

# Pre-process Seurat object ----------------------------------------------------
DefaultAssay(dat) <- "RNA"
dat <- NormalizeData(dat)
dat <- FindVariableFeatures(dat, selection.method = "vst", nfeatures = 2000)
dat <- ScaleData(dat)
dat <- RunPCA(dat)
dat <- RunUMAP(dat, dims = 1:10)
dat <- FindNeighbors(dat, dims = 1:10)
dat <- FindClusters(dat, resolution = 0.1)

# pK Identification ------------------------------------------------------------
sweep.res.list_dat <- paramSweep_v3(dat, PCs = 1:10, sct = FALSE)
sweep.stats_dat <- summarizeSweep(sweep.res.list_dat, GT = FALSE)
bcmvn_dat <- find.pK(sweep.stats_dat)

# Homotypic Doublet Proportion Estimate ----------------------------------------
homotypic.prop <- modelHomotypic(dat@meta.data$seurat_clusters)
nExp_poi <- round(exp_double_rate * nrow(dat@meta.data))
nExp_poi.adj <- round(nExp_poi*(1 - homotypic.prop))

# Run DoubletFinder with varying classification stringencies -------------------
dat <- doubletFinder_v3(dat, 
                        PCs = 1:10,
                        pK = 0.2,
                        nExp = nExp_poi.adj,
                        reuse.pANN = FALSE, 
                        sct = FALSE)

# Save the classified doublets
doublets_25C_Rep1 <- dat@meta.data |>
  base::subset(DF.classifications_0.25_0.2_59 == "Doublet") |>
  base::rownames()

# 25C_Rep2 ---------------------------------------------------------------------

# Load data
dat <- readRDS(
  here::here("data/processed/seurat_object/01_dat_10x_filtered.rds")) |>
  base::subset(sample_name == "25C_Rep2")

# Pre-process Seurat object ----------------------------------------------------
DefaultAssay(dat) <- "RNA"
dat <- NormalizeData(dat)
dat <- FindVariableFeatures(dat, selection.method = "vst", nfeatures = 2000)
dat <- ScaleData(dat)
dat <- RunPCA(dat)
dat <- RunUMAP(dat, dims = 1:10)
dat <- FindNeighbors(dat, dims = 1:10)
dat <- FindClusters(dat, resolution = 0.1)

# pK Identification ------------------------------------------------------------
sweep.res.list_dat <- paramSweep_v3(dat, PCs = 1:10, sct = FALSE)
sweep.stats_dat <- summarizeSweep(sweep.res.list_dat, GT = FALSE)
bcmvn_dat <- find.pK(sweep.stats_dat)

# Homotypic Doublet Proportion Estimate ----------------------------------------
homotypic.prop <- modelHomotypic(dat@meta.data$seurat_clusters)
nExp_poi <- round(exp_double_rate * nrow(dat@meta.data))
nExp_poi.adj <- round(nExp_poi*(1 - homotypic.prop))

# Run DoubletFinder with varying classification stringencies -------------------
dat <- doubletFinder_v3(dat, 
                        PCs = 1:10,
                        pK = 0.005,
                        nExp = nExp_poi.adj,
                        reuse.pANN = FALSE, 
                        sct = FALSE)

# Save the classified doublets
doublets_25C_Rep2 <- dat@meta.data |>
  base::subset(DF.classifications_0.25_0.005_67 == "Doublet") |>
  base::rownames()

doublets <- c(doublets_18C_Rep1,
              doublets_18C_Rep2,
              doublets_25C_Rep1,
              doublets_25C_Rep2)

# Save doublets
saveRDS(doublets, here::here("data/processed/qc/doublet_finder.rds"))


# ------------------------------------------------------------------------------
# Remove multiplets
# TS O'Leary
# ------------------------------------------------------------------------------



# Load multiplet barcodes
doublet_amulet <- readRDS(
  here::here("data/processed/qc/doublet_amulet.rds")
) |>
  dplyr::filter(q.value < 0.05) |>
  base::rownames()
doublet_finder <- readRDS(here::here("data/processed/qc/doublet_finder.rds"))

# Load seurat object
dat <- readRDS(
  here::here("data/processed/seurat_object/02_dat_initial_cluster.rds")
) 

# Output dir
out_dir <- "data/processed/seurat_object"

# Mark amulet described multiplet barcodes
dat$doublet_amulet <- ifelse(colnames(dat) %in% doublet_amulet, 
                             "Multiplet",
                             "Singlet")

# Mark amulet described multiplet barcodes
dat$doublet_finder <- ifelse(colnames(dat) %in% doublet_finder,
                             "Multiplet",
                             "Singlet")

dat@meta.data <- dat@meta.data |> 
  mutate(multiplet_call = case_when(doublet_amulet == "Singlet" & 
                                      doublet_finder == "Singlet" ~ "Singlet",
                                    .default = "Multiplet")) 



# Remove putative multiplets from the data set -----

# Make clusters at a more fine resolution
dat <- FindClusters(dat, resolution = 1.0)

# Plot of subclusters
DimPlot(dat, label = TRUE) + 
  theme_void() +
  theme(title = element_blank(),
        legend.position = "none")

# Plot of individual multiplets marked
DimPlot(dat, group = "multiplet_call") +
  scale_color_manual(name = element_blank(),
                     values = c("firebrick", "grey80")) +
  theme_void() +
  theme(title = element_blank(),
        legend.position = "bottom")

# Identify proportion of multiplets within subclusters
multiplet_prop <- dat@meta.data |> 
  group_by(seurat_clusters, multiplet_call) |> 
  tally() |> 
  pivot_wider(names_from = multiplet_call,
              values_from = n) |> 
  mutate(percent_multiplet = Multiplet/(Singlet + Multiplet)) |> 
  arrange(desc(percent_multiplet))

# Identify all clusters with greater than 15% multiplets
multiplet_clusters <- multiplet_prop |> 
  filter(percent_multiplet > 0.15) |> 
  pull(seurat_clusters)

# Identify all cells marked individually as multiplets or cells that are in 
#   clusters with at least 15% multiplets
dat@meta.data <- dat@meta.data |> 
  mutate(multiplet_rm = 
           case_when(multiplet_call == "Multiplet" |
                       seurat_clusters %in% multiplet_clusters ~ "Multiplet",
                     .default = "Singlet"))

# Save dat with multiplets marked
saveRDS(dat, here::here(out_dir, "05_dat_multiplets.rds"))

# Plot all multiplets that will be removed
DimPlot(dat, group = "multiplet_rm") +
  scale_color_manual(name = element_blank(),
                     values = c("firebrick", "grey80")) +
  theme_void() +
  theme(title = element_blank(),
        legend.position = "bottom")

# Remove all those multiplets
dat <- dat |>
  base::subset(multiplet_rm != "Multiplet")

# Remove cluster info from that initial clustering
dat@meta.data <- dat@meta.data |>
  dplyr::select(-seurat_clusters) |>
  dplyr::select(!contains("snn")) |>
  dplyr::select(!contains("doublet")) |> 
  dplyr::select(!contains("multiplet"))

