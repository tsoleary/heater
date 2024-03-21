# ------------------------------------------------------------------------------
# Find enriched TF motifs from lists of DARs
# TS O'Leary
# ------------------------------------------------------------------------------

# Load libraries
library(tidyverse)
library(Signac)
library(Seurat)
library(JASPAR2022)
library(TFBSTools)
library(BSgenome.Dmelanogaster.UCSC.dm6)

# Load data
dat <- readRDS(here::here("data/processed/seurat_object/10_dat_linked.rds"))
dars <- readRDS(here::here("output/dars/dars.rds")) |>
  dplyr::filter(p_val_adj < 0.05) |>
  rownames()
dars_celltype <- readRDS(here::here("output/dars/dars_cell-type.rds")) |> 
  dplyr::filter(p_val_adj < 0.05) |>
  pull(region)
# meso <- readRDS(here::here("output/dars/dars_cell-type.rds")) |> 
#   dplyr::filter(
#     p_val_adj < 0.05,
#     cell_type == "mesoderm prim.") |> 
#   pull(region)
# vnc <- readRDS(here::here("output/dars/dars_cell-type.rds")) |> 
#   dplyr::filter(
#     p_val_adj < 0.05,
#     cell_type == "ventral nerve cord prim.") |> 
#   pull(region)
# blue <- read_csv(here::here("output/wgcna/atac/blue.csv")) |> pull()
# brown <- read_csv(here::here("output/wgcna/atac/brown.csv")) |> pull()

# Output dir
out_dir <- "output/tf_motif"

# # Linked peaks to genes
# GetLinkedGenes(dat, dars, min.abs.score = 0.1)

# Set up data
DefaultAssay(dat) <- "peaks"
seqlevelsStyle(BSgenome.Dmelanogaster.UCSC.dm6) <- "NCBI"

# Get a list of motif position frequency matrices from the JASPAR database
pfm <- getMatrixSet(
  x = JASPAR2022,
  opts = list(collection = "CORE", 
              species = "Drosophila melanogaster",
              tax_group = "insects", 
              all_versions = FALSE)
)

# Add motif information
dat <- AddMotifs(
  object = dat,
  genome = BSgenome.Dmelanogaster.UCSC.dm6,
  pfm = pfm
)

# Save dat with motifs
saveRDS(dat, here::here("data/processed/seurat_object/14_dat_tf_motif.rds"))

# Test motif enrichment -----

# DARs from pseudobulk
enriched_motifs <- FindMotifs(
  object = dat,
  background = 10800,
  assay = "peaks",
  features = dars) |> 
  dplyr::filter(p.adjust < 0.05)

# DARs per cell-type
enriched_motifs_celltype <- FindMotifs(
  object = dat,
  background = 10000,
  assay = "peaks",
  features = dars_celltype) |> 
  dplyr::filter(p.adjust < 0.05)

# DARs mesoderm
enriched_motifs_meso <- FindMotifs(
  object = dat,
  background = 10000,
  assay = "peaks",
  features = meso) |>
  dplyr::filter(p.adjust < 0.05)

# DARs mesoderm
enriched_motifs_vnc <- FindMotifs(
  object = dat,
  background = 10000,
  assay = "peaks",
  features = vnc) |>
  dplyr::filter(p.adjust < 0.05)

# DARs mesoderm
enriched_motifs_meso <- FindMotifs(
  object = dat,
  background = 10000,
  assay = "peaks",
  features = meso) |>
  dplyr::filter(p.adjust < 0.05)

# Turquoise
enriched_motifs_blue <- FindMotifs(
  object = dat,
  background = 10000,
  assay = "peaks",
  features = blue) |> 
  dplyr::filter(p.adjust < 0.05)

# Brown
enriched_motifs_brown <- FindMotifs(
  object = dat,
  background = 10000,
  assay = "peaks",
  features = brown) |> 
  dplyr::filter(p.adjust < 0.05)

# Motif plots -----

# Pseudobulk
MotifPlot(
  object = dat,
  motifs = rownames(enriched_motifs)) +
  theme(axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_blank())

# Cell-type
MotifPlot(
  object = dat,
  motifs = rownames(enriched_motifs_celltype)) +
  theme(axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_blank())

# Motif plots
MotifPlot(
  object = dat,
  motifs = enriched_motifs_turquoise$motif) +
  theme(axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_blank())

# Map peaks to enriched motifs
all_enriched_motifs <- bind_rows(list(bulk = enriched_motifs, 
                    cell_type = enriched_motifs_celltype, 
                    vnc = enriched_motifs_celltype,
                    brown = enriched_motifs_brown, 
                    turquoise = enriched_motifs_turquoise),
                    .id = "id") |> 
  tibble()

# Set of GRanges object of the pseudobulk dars
dars_granges <- dars |> 
  tibble() |> 
  separate(dars, into = c("seqnames", "start", "end"), sep = "-") |> 
  mutate(strand = "*") |> 
  makeGRangesFromDataFrame()

# Initialize dars column
all_enriched_motifs$dars <- NULL

# Loop through all motifs
for (i in 1:nrow(all_enriched_motifs)){
  
  # Get all genomic ranges for a motif
  motif_ranges <- 
    dat@assays$peaks@motifs@positions[all_enriched_motifs$motif[i]]

  # Find the dars that include the motifs
  all_enriched_motifs$dars[i] <- subsetByOverlaps(dars_granges, motif_ranges) |> 
    as_tibble() |> 
    mutate(dars = paste(seqnames, start, end, sep = "-")) |> 
    pull(dars) |> 
    str_flatten(",")
  
}

# Pivot DARs column to longer
all_enriched_motifs <- all_enriched_motifs |> 
  separate_longer_delim(dars, delim = ",")

# Get DARs data with fold-change
dars_df <- readRDS(here::here("output/dars/dars.rds")) |> 
  filter(p_val_adj < 0.05) |> 
  rownames_to_column("region") 

# Join motifs with their DAR statistics
dars_df <- full_join(
  all_enriched_motifs,
  dars_df, 
  by = c("dars" = "region")
)

# Save these dars with motif information
saveRDS(dars_df, here::here(out_dir, "dars_df_motif.rds"))

################################################################################
# # CHROM VAR NOT WORKING -----
# # Set up
# #BiocParallel::register(BiocParallel::MulticoreParam(8)) 
# # Per-cell motif activity score 
# DefaultAssay(dat) <- "peaks"
# dat <- RunChromVAR(
#   object = dat,
#   genome = BSgenome.Dmelanogaster.UCSC.dm6,
# )
# 
# # Set assay for this 
# DefaultAssay(dat) <- "chromvar"
# 
# # look at the activity of Mef2c
# p2 <- FeaturePlot(
#   object = mouse_brain,
#   features = "MA0497.1",
#   min.cutoff = 'q10',
#   max.cutoff = 'q90',
#   pt.size = 0.1
# )
