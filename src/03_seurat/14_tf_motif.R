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
  dplyr::filter(padj < 0.05) |> 
  pull(region)
meso <- readRDS(here::here("output/dars/dars_cell-type.rds")) |> 
  dplyr::filter(
    padj < 0.05,
    cell_type == "mesoderm prim.") |> 
  pull(region)
vnc <- readRDS(here::here("output/dars/dars_cell-type.rds")) |> 
  dplyr::filter(
    padj < 0.05,
    cell_type == "ventral nerve cord prim.") |> 
  pull(region)
turquoise <- read_csv(here::here("output/wgcna/atac/turquoise.csv")) |> pull()
brown <- read_csv(here::here("output/wgcna/atac/brown.csv")) |> pull()

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

# Save dat with hdWGCNA in it for plotting
saveRDS(dat, here::here("data/processed/seurat_object/14_dat_tf_motif.rds"))

# Test motif enrichment -----

# DARs from pseudobulk
enriched_motifs <- FindMotifs(
  object = dat,
  background = 10000,
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
enriched_motifs_turquoise <- FindMotifs(
  object = dat,
  background = 10000,
  assay = "peaks",
  features = turquoise) |> 
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
# Plotting ---------------------------------------------------------------------
# NEED TO MAKE THIS A SEPERATE SCRIPT!!
################################################################################

# Load libraries
require(cowplot)

# Load plot themes 
source(here::here("src/04_plots/00_plot_themes.R"))

# Output figure directory
fig_dir <- "output/figs/tf_motif"

# Bulk enriched motif
dars_df |>
  filter(id == "bulk") |> 
  ggplot() +
  geom_histogram(aes(x = avg_log2FC),
                 color = "grey20",
                 fill = "grey80") +
  geom_vline(xintercept = 0,
             color = "grey80") +
  facet_wrap(~motif.name) + 
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = "Average log2(fold-change)\n18°C vs 25°C",
       y = "Number of DARs") +
  theme_minimal_hgrid() +
  theme(strip.background = element_rect(fill = "grey90"))

# Save plot
ggsave(here::here(fig_dir, "pseudobulk_motif_dars_vfl_hist.png"),
       height = 10,
       width = 15,
       units = "cm")

# Brown histogram of enriched motifs DARs Log2FC
dars_df |>
  filter(id == "brown") |> 
  ggplot() +
  geom_histogram(aes(x = avg_log2FC),
                 color = "grey20",
                 fill = "grey80") +
  geom_vline(xintercept = 0,
             color = "grey80") +
  facet_wrap(~motif.name) + 
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = "Average log2(fold-change)\n18°C vs 25°C",
       y = "Number of DARs") +
  theme_minimal_hgrid() +
  theme(strip.background = element_rect(fill = "grey90"))

# Save plot
ggsave(here::here(fig_dir, "brown_motif_dars_hist.png"),
       height = 30,
       width = 30,
       units = "cm")

# Turquoise histogram of enriched motifs DARs Log2FC
dars_df |>
  filter(id == "turquoise") |> 
  ggplot() +
  geom_histogram(aes(x = avg_log2FC),
                 color = "grey20",
                 fill = "grey80") +
  geom_vline(xintercept = 0,
             color = "grey80") +
  facet_wrap(~motif.name) + 
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = "Average log2(fold-change)\n18°C vs 25°C",
       y = "Number of DARs") +
  theme_minimal_hgrid() +
  theme(strip.background = element_rect(fill = "grey90"))

# Save plot
ggsave(here::here(fig_dir, "turquoise_motif_dars_hist.png"),
       height = 35,
       width = 35,
       units = "cm")

# VNC histogram of enriched motifs DARs Log2FC

full_join(
  dars_df |>
    filter(id == "vnc"),
  readRDS(here::here("output/dars/dars_cell-type.rds")) |> 
    dplyr::filter(
      padj < 0.05,
      cell_type == "ventral nerve cord prim."),
  by = c("dars" =  "region")) |> 
  ggplot() +
  geom_histogram(aes(x = avg_log2FC_18_25),
                 color = "grey20",
                 fill = "grey80") +
  geom_vline(xintercept = 0,
             color = "grey80") +
  facet_wrap(~motif.name) + 
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = "Average log2(fold-change)\n18°C vs 25°C",
       y = "Number of DARs") +
  theme_minimal_hgrid() +
  theme(strip.background = element_rect(fill = "grey90"))

# Save plot
ggsave(here::here(fig_dir, "vnc_motif_dars_hist.png"),
       height = 35,
       width = 35,
       units = "cm")

# Cell-type
dars_df |>
  filter(id == "cell_type") |> 
  ggplot() +
  geom_histogram(aes(x = avg_log2FC),
                 color = "grey20",
                 fill = "grey80") +
  geom_vline(xintercept = 0,
             color = "grey80") +
  facet_wrap(~motif.name) + 
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = "Average log2(fold-change)\n18°C vs 25°C",
       y = "Number of DARs") +
  theme_minimal_hgrid() +
  theme(strip.background = element_rect(fill = "grey90"))

# Save plot
ggsave(here::here(fig_dir, "cell-type_motif_dars_hist.png"),
       height = 35,
       width = 35,
       units = "cm")


# Follow up on Zelda TF binding site differential accessibility ----------------

# Loop through pseudobulf vfl DARs to see increase in 18°C accessibility
zelda <- dars_df |> 
  filter(id == "bulk") |> 
  filter(motif.name == "vfl")

Idents(dat) <- "acc_temp"

for(i in 1:nrow(zelda)) {
  print(i)
  
  cov_plot <- CoveragePlot(
    object = dat,
    region = zelda$dars[i],
    assay = "peaks",
    annotation = FALSE,
    links = FALSE,
    peaks = FALSE,
    extend.upstream = 5000,
    extend.downstream = 5000) &
    scale_fill_manual(values = acc_colors) & 
    scale_y_continuous(expand = c(0, 0.05),
                       name = "Normalized signal") &
    theme_minimal_hgrid() &
    theme(legend.position = "none",
          strip.background = element_blank(),
          strip.text = element_blank())
  
  peak_rect <- annotate("rect",
                        xmin = as.numeric(str_remove_all(
                          str_extract(zelda$dars[i], "-[:digit:]+-"), "-")),
                        xmax = as.numeric(str_remove_all(
                          str_extract(zelda$dars[i], "-[:digit:]+$"), "-")),
                        ymin = -Inf,
                        ymax = Inf,
                        fill = "grey75",
                        alpha = 0.5)
  
  cov_plot$layers <- c(peak_rect, cov_plot$layers)
  
  gene_plot <- AnnotationPlot(
    object = dat,
    region = zelda$dars[i],
    extend.upstream = 5000,
    extend.downstream = 5000) & 
    theme_cowplot() &
    theme(axis.line.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.y = element_blank(),
          axis.title.y = element_blank())
  
  peak_plot <- PeakPlot(
    object = dat,
    region = zelda$dars[i],
    assay = "peaks",
    extend.upstream = 5000,
    extend.downstream = 5000) &
    theme_cowplot() &
    theme(axis.line.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.y = element_blank(),
          axis.title.y = element_blank(),
          legend.position = "none")
  
  CombineTracks(
    plotlist = list(
      cov_plot,
      peak_plot,
      gene_plot
    ),
    heights = c(8, 0.1, 2)
  )
  
  # Save plot
  ggsave(here::here(fig_dir, "zelda", paste0(zelda$dars[i], ".png")),
         height = 35,
         width = 35,
         units = "cm")
}




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
