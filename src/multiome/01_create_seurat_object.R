# ------------------------------------------------------------------------------
# Load data and set up Seurat object
# March 27, 2023
# TS O'Leary
# ------------------------------------------------------------------------------

# Description 
# Load in the RNA-Seq count data and the ATAC-Seq peak fragment data

# Load libraries
library(Seurat)

# Load the RNA-Seq and ATAC-Seq data
counts <- Read10X_h5("../vignette_data/multiomic/pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix.h5")
fragpath <- "../vignette_data/multiomic/pbmc_granulocyte_sorted_10k_atac_fragments.tsv.gz"

# Get gene annotations for hg38
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotation) <- "UCSC"

# Create a Seurat object containing the RNA adata
dat <- CreateSeuratObject(
  counts = counts$`Gene Expression`,
  assay = "RNA"
)

# Create ATAC assay and add it to the object
dat[["ATAC"]] <- CreateChromatinAssay(
  counts = counts$Peaks,
  sep = c(":", "-"),
  fragments = fragpath,
  annotation = annotation
)

# Save data
saveRDS(dat, "data/processed/seurat_object_01_initial.rds")