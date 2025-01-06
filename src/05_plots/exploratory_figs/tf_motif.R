# ------------------------------------------------------------------------------
# Transcription factor motif related plots
# TS O'Leary
# ------------------------------------------------------------------------------

# Load libraries
library(Signac)
library(Seurat)
library(tidyverse)
library(cowplot)
library(IRanges)

# Load plot themes 
source(here::here("src/04_plots/00_plot_themes.R"))

# Load data
dat <- readRDS(here::here("data/processed/seurat_object/14_dat_tf_motif.rds"))
dars_df <- readRDS(here::here("output/tf_motif/dars_df_motif.rds"))
dars <- readRDS(here::here("output/dars/dars.rds")) |>
  dplyr::filter(p_val_adj < 0.05) |>
  rownames()

# Output figure directory
fig_dir <- "output/figs/tf_motif"

# Bulk enriched motif
dars_df |>
  filter(id == "bulk") |> 
  ggplot() +
  geom_histogram(aes(x = avg_log2FC),
                 color = "grey20",
                 fill = "grey80",
                 bins = 25) +
  geom_vline(xintercept = 0,
             color = "grey80") +
  #facet_wrap(~motif.name) + 
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(limits = c(-0.7, 0.7),
                     breaks = c(-0.5, -0.25, 0, 0.25, 0.5)) +
  labs(x = "log2(fold-change)\n18°C vs 25°C",
       y = "Number of DARs") +
  theme_minimal_hgrid() 

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


# Brown histogram of enriched motifs DARs Log2FC
dars_df |>
  filter(id == "turquoise") |> 
  filter(motif.name == "pan") |> 
  ggplot() +
  geom_histogram(aes(x = avg_log2FC),
                 color = "grey20",
                 fill = "grey80") +
  geom_vline(xintercept = 0,
             color = "grey80") +
  facet_wrap(~motif) + 
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = "Average log2(fold-change)\n18°C vs 25°C",
       y = "Number of DARs") +
  theme_minimal_hgrid() +
  theme(strip.background = element_rect(fill = "grey90"))

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


# Heat shock factor motifs

hsf_motifs <- dat@assays$peaks@motifs@positions["MA1458.1"]


hsf_dars <- IRanges::subsetByOverlaps(dars_granges, hsf_motifs) |> 
  as_tibble() |> 
  mutate(dars = paste(seqnames, start, end, sep = "-")) |> 
  pull(dars)

readRDS(here::here("output/dars/dars.rds")) |> 
  filter(p_val_adj < 0.05) |> 
  rownames_to_column("region") |> 
  filter(region %in% hsf_dars) |> 
  ggplot() +
  geom_histogram(aes(x = avg_log2FC),
                 color = "grey20",
                 fill = "grey80",
                 bins = 15) +
  geom_vline(xintercept = 0,
             color = "grey80") +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(limits = c(-0.7, 0.7),
                     breaks = c(-0.5, -0.25, 0, 0.25, 0.5)) +
  labs(x = "log2(fold-change)\n18°C vs 25°C",
       y = "Number of DARs") +
  theme_minimal_hgrid() 


# Save plot
ggsave(here::here(fig_dir, "pseudobulk_motif_dars_hsf_hist.png"),
       height = 10,
       width = 15,
       units = "cm")

# Pseudobulk
MotifPlot(
  object = dat,
  motifs = "MA1458.1") +
  theme(axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_blank())


Idents(dat) <- "acc_temp"


for(i in 1:length(hsf_dars)) {
  print(i)
  
  cov_plot <- CoveragePlot(
    object = dat,
    region = hsf_dars[i],
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
                          str_extract(hsf_dars[i], "-[:digit:]+-"), "-")),
                        xmax = as.numeric(str_remove_all(
                          str_extract(hsf_dars[i], "-[:digit:]+$"), "-")),
                        ymin = -Inf,
                        ymax = Inf,
                        fill = "grey75",
                        alpha = 0.5)
  
  cov_plot$layers <- c(peak_rect, cov_plot$layers)
  
  gene_plot <- AnnotationPlot(
    object = dat,
    region = hsf_dars[i],
    extend.upstream = 5000,
    extend.downstream = 5000) & 
    theme_cowplot() &
    theme(axis.line.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.y = element_blank(),
          axis.title.y = element_blank())
  
  peak_plot <- PeakPlot(
    object = dat,
    region = hsf_dars[i],
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
  ggsave(here::here(fig_dir, "hsf", paste0(hsf_dars[i], ".png")),
         height = 35,
         width = 35,
         units = "cm")
}


hsf_up_18 <- readRDS(here::here("output/dars/dars.rds")) |> 
  filter(p_val_adj < 0.05) |> 
  rownames_to_column("region") |> 
  filter(region %in% hsf_dars) |> 
  filter(avg_log2FC > 0) |> 
  pull(region)

hsf_up_25 <- readRDS(here::here("output/dars/dars.rds")) |> 
  filter(p_val_adj < 0.05) |> 
  rownames_to_column("region") |> 
  filter(region %in% hsf_dars) |> 
  filter(avg_log2FC < 0) |> 
  pull(region)

GetLinkedGenes(
  dat,
  features = hsf_up_18,
  assay = "peaks",
  min.abs.score = 0.1
) |> 
  str_flatten_comma()

GetLinkedGenes(
  dat,
  features = hsf_up_25,
  assay = "peaks",
  min.abs.score = 0.1
) |> 
  str_flatten_comma()

