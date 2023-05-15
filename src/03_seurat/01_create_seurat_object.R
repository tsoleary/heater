# ------------------------------------------------------------------------------
# Load data and set up Seurat object
# May 11, 2023
# TS O'Leary
# ------------------------------------------------------------------------------

# Description 
# Load in the RNA-Seq count data and the ATAC-Seq peak fragment data and save as
# an rds file to be used for quality control and filtering.

# Load libraries
library(tidyverse)
library(Seurat)
library(Signac)
library(AnnotationHub)

# Load the RNA-Seq and ATAC-Seq data -------------------------------------------
data_dir <- here::here("data/processed/seq/all/outs/")
raw_feature_dir <- paste0(data_dir, "raw_feature_bc_matrix")
frag_path <- paste0(data_dir, "atac_fragments.tsv.gz")
counts <- Read10X(raw_feature_dir)

# Get gene annotations for dm6 -------------------------------------------------
BDGP6.32 <- query(
  AnnotationHub(), 
  c("EnsDb", "Drosophila melanogaster", "109")
)
annotation <- GetGRangesFromEnsDb(BDGP6.32[[1]])
seqlevelsStyle(annotation) <- "UCSC"

# Meta data --------------------------------------------------------------------
# Create meta data for all cells based on the GEM well suffix. From 10X Genomics
# "This number, which indicates which GEM well the barcode sequence came from, 
# is called the GEM well suffix. The numbering of the GEM wells will reflect the 
# order that the GEM wells were provided in the Aggregation CSV."
meta_data <- tibble(cellnames = colnames(counts[[1]])) %>%
  mutate(orig.ident = str_extract(cellnames, "[1-4]")) %>%
  full_join(tibble::tribble(
    ~orig.ident, ~sample_name, ~acc_temp,
    "1", "18C_Rep1", "18C",
    "2", "18C_Rep2", "18C",
    "3", "25C_Rep1", "25C",
    "4", "25C_Rep2", "25C"
  ), 
  by = "orig.ident") %>%
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
  # Validate fragments arg makes sure all cells can be found in the frags file
  # However, there may be cells/barcodes that are present in the assay, but are 
  # not represented in the fragment file. For example, RNA-seq reads that are
  # from GEM barcodes that have no fragments associated with them. Check with 
  # Seth about this.
  validate.fragments = FALSE
)

# Save raw Seurat object
saveRDS(dat, here::here("data/processed/seurat_object/01_dat_raw.rds"))
