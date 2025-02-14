# ------------------------------------------------------------------------------
# Load Data: Load in raw RNA-Seq count data and the ATAC-Seq peak fragment data 
# and save as an rds file to be used for quality control and filtering
# TS O'Leary
# ------------------------------------------------------------------------------

# Load libraries
library(tidyverse)
library(Seurat)
library(Signac)

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

# Save raw Seurat object
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
