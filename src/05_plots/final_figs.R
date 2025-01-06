# ------------------------------------------------------------------------------
# Final figures
# TS O'Leary
# ------------------------------------------------------------------------------

# Load libraries
library(tidyverse)
library(Seurat)
library(Signac)
library(cowplot)
library(hdWGCNA)

# Output fig dir
fig_dir <- "output/figs/final"

# Load plot themes & modified functions
source(here::here("src/05_plots/00_plot_themes.R"))
source(here::here("src/05_plots/plot_functions.R"))

# Phenotype data
pheno <- read_csv(
  here::here("data/raw/pheno/acc_hs_survival.csv")) |> 
  mutate(acc_temp_factor = factor(paste0(acc_temp, "°C")))

# Load raw data
dat_raw <- readRDS(
  here::here("data/processed/seurat_object/00_dat_raw.rds")
)

# Load final seurat object
dat <- readRDS(
  here::here("data/processed/seurat_object/14_dat_tf_motif.rds")
)

# Recode sample name for prettier plotting
dat@meta.data <- dat@meta.data |> 
  mutate(sample_name_fig = case_when(sample_name == "18C_Rep1" ~ "18°C Rep 1",
                                     sample_name == "18C_Rep2" ~ "18°C Rep 2",
                                     sample_name == "25C_Rep1" ~ "25°C Rep 1",
                                     sample_name == "25C_Rep2" ~ "25°C Rep 2")) |> 
  mutate(sample_name_fig = factor(sample_name_fig, levels = c(c("25°C Rep 2",
                                                        "25°C Rep 1",
                                                        "18°C Rep 2",
                                                        "18°C Rep 1"))))
# Create dataframe with cell names as column for plotting
meta <- dat@meta.data |> 
  rownames_to_column("cell")

# Load all cell-type specific acclimation comparisons
degs <- readRDS(here::here("output/degs/degs_cell-type_MAST.rds")) |> 
  mutate(sig = ifelse(abs(avg_log2FC_18_25) >= 0.25 & p_val_adj < 0.05 &
                        (pct.18 >= 0.10 | pct.25 >= 0.10), 
         TRUE, FALSE)) 
dars <- readRDS(here::here("output/dars/dars_cell-type_MAST.rds")) |> 
  mutate(sig = ifelse(abs(avg_log2FC_18_25) >= 0.25 & p_val_adj < 0.05 &
                        (pct.18 >= 0.10 | pct.25 >= 0.10), 
                      TRUE, FALSE)) 

dars_bulk <- readRDS(here::here("output/dars/dars_bulk.rds")) |> 
  as.data.frame() |> 
  mutate(sig = ifelse(abs(avg_log2FC) >= 0.25 & p_val_adj < 0.05 &
                        (pct.1 >= 0.10 | pct.2 >= 0.10), 
                      TRUE, FALSE)) 
degs_bulk <- readRDS(here::here("output/degs/pseudobulk_DESeq_res.rds")) 

# eRegulon data
eReg <- read_csv(here::here("data/processed/scenic/eRegulon_data.csv"))

# ------------------------------------------------------------------------------
# Manuscript figures -----------------------------------------------------------
# ------------------------------------------------------------------------------


# Figure 1 ---------------------------------------------------------------------

# Acclimation egg hatching phenotype
f1 <- pheno |>
  ggplot(aes(x = acc_temp,
             y = n_hatched/n_eggs*100,
             group = acc_temp_factor,
             fill = acc_temp_factor)) +
  geom_boxplot(width = 2,
               color = "grey20",
               outlier.shape = NA) +
  ggbeeswarm::geom_beeswarm(size = 2,
                            cex = 3.75,
                            shape = 21,
                            color = "grey20",
                            alpha = 0.9) +
  xlab("Acclimation temperature") +
  scale_y_continuous(name = "Hatching success\nafter 45 min @ 39°C",
                     limits = c(0, 100),
                     expand = expansion(mult = c(0, 0.05)),
                     labels = function(x) paste0(x, "%")) + 
  scale_x_continuous(breaks = c(18, 25, 30),
                     labels = c("18°C", "25°C","30°C" )) +
  scale_fill_manual(values = c(acc_colors, "#dc8383")) +
  cowplot::theme_minimal_hgrid() +
  theme(legend.position = "none",
        text = element_text(size = 8),
        axis.title = element_text(size = 8),        
        axis.title.x = element_text(size = 8),       
        axis.title.y = element_text(size = 8),     
        axis.text.x = element_text(size = 8),        
        axis.text.y = element_text(size = 8))

# Save
ggsave(here::here(fig_dir, "fig_1.png"),
       plot = f1,
       height = 6,
       width = 8.7,
       units = "cm")
ggsave(here::here(fig_dir, "fig_1.pdf"),
       plot = f1,
       height = 6,
       width = 8.7,
       units = "cm")

# Figure 2 ---------------------------------------------------------------------

# UMAP acclimation treatment overlay
p_overlay <- DimPlot(
  dat,
  group.by = "acc_temp",
  reduction = "umap") +
  scale_color_manual(name = element_blank(),
                     values = acc_colors) +
  labs(title = element_blank()) +
  theme_void() +
  theme(legend.position = c(0.8, 0.2)) +
  guides(color = guide_legend(byrow = TRUE,
                              nrow = 2,
                              override.aes = list(size = 3)))

# UMAP clusters 
p_clust <- DimPlot(
  dat,
  group.by = "seurat_clusters",
  reduction = "umap") +
  labs(title = element_blank()) +
  scale_color_manual(
    values = scales::hue_pal()(12) |> 
      colorspace::desaturate(amount = 0.5)
  ) +
  theme_void() +
  theme(legend.position = "none")

# Label seurat clusters
p_clust_labeled <- LabelClusters(
  p_clust, 
  id = "seurat_clusters",
  fontface = "bold",
  color = "grey20",
  alpha = 0.8,
  box = TRUE,
  repel = TRUE,
  family = "Arial"
)


# UMAP cell_type
# Set cell-type as the primary identity
Idents(dat) <- "cell_type"

# UMAP colored by cell-types
p_cell_type <- DimPlot(dat, 
                reduction = "umap",
                repel = TRUE,
                label.color = "grey20",
                pt.size = 0.2) +
  theme_void() +
  scale_color_manual(
    values = c("grey90",
               "#ADD9F4",
               "#57A4B2",
               "#D39C9E",
               "#FEF29A",
               "#F9DCEE",
               "#819FC5",
               "#A7BF9B",
               "#bfa3a4"))

# Get legend for saving later
legend_cell_types <- 
  get_plot_component(p_cell_type, "guide-box", return_all = TRUE)[[1]]

# Combining the top row
p_umaps <- cowplot::plot_grid(
  p_overlay,
  p_clust_labeled,
  p_cell_type + theme(legend.position = "none"), # Combine without legend
  rel_widths = c(1, 1, 1),
  nrow = 1,
  labels = c("A", "B", "C"),
  label_colour = "grey20"
)

# Heatmaps
Idents(dat) <- "cell_type"

# Aggregate marker genes across all cell_types
dat_agg <- AggregateExpression(
  dat,
  assay = c("SCT"),
  return.seurat = TRUE,
  features = markers$gene)

# Scale data for heatmap
dat_agg <- ScaleData(dat_agg)

# Run heatmap on scaled.data
hmap_exp <- pheatmap::pheatmap(
  t(dat_agg@assays$SCT$scale.data),
  colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name = "RdBu")))(100),
  treeheight_col = 20,
  cluster_rows = FALSE,
  treeheight_row = 10,
  border_color = NA,
  show_colnames = FALSE,
  show_rownames = FALSE) 

gene_order <- colnames(t(dat_agg@assays$SCT$scale.data))[hmap_exp$tree_col$order]

hmap_exp <- hmap_exp |> 
  ggplotify::as.ggplot() +
  theme(legend.key.size = unit(0.25, "cm"))

# For ATAC ----
DefaultAssay(dat) <- "peaks"

# Save all peaks that are linked to marker genes for plotting purposes
linked_marker_peaks <- Links(dat) |>
  as_tibble() |>
  filter(gene %in% markers$gene) |>
  group_by(gene) |> 
  slice_min(pvalue)

linked_peaks_order <- linked_marker_peaks |> 
  mutate(gene = factor(gene, levels = gene_order)) |> 
  arrange(gene) |> 
  pull(peak)

# Aggregate marker peaks across all cell_types
dat_agg_atac <- AggregateExpression(
  dat,
  assay = "peaks",
  return.seurat = TRUE,
  features = linked_peaks_order)

# Scale data
dat_agg_atac <- ScaleData(dat_agg_atac)

# Do min-max normalization
atac_min_max_norm <- caret::preProcess(dat_agg_atac@assays$peaks$scale.data, 
                                       method = c("range")) |> 
  predict(dat_agg_atac@assays$peaks$scale.data)

# Run heatmap on scaled.data
hmap_atac <- pheatmap::pheatmap(
  t(atac_min_max_norm),
  pals::ocean.thermal(n = 100),
  treeheight_col = 20,
  treeheight_row = 10,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  border_color = NA,
  show_rownames = FALSE,
  show_colnames = FALSE) |> 
  ggplotify::as.ggplot() &
  theme(legend.key.size = unit(0.25, "cm"),
        legend.position = "bottom")

# Heatmaps
p_heatmaps <- cowplot::plot_grid(
  hmap_exp,
  hmap_atac, 
  labels = c("D", "E"),
  hjust = 1,
  nrow = 2)

# Save individual parts of the figures ---
ggsave(here::here(fig_dir, "fig_2/umaps.svg"),
       plot = p_umaps,
       height = 10,
       width = 30,
       units = "cm")
ggsave(here::here(fig_dir, "fig_2/umap_cell_type_legend.svg"),
       plot = ggdraw(legend_cell_types),
       height = 10,
       width = 10,
       units = "cm")
ggsave(here::here(fig_dir, "fig_2/heatmaps.svg"),
       plot = p_heatmaps,
       height = 10,
       width = 30,
       units = "cm")



# Figure 3 ---------------------------------------------------------------------

# Load all DARs
dars_motif_df <- readRDS(here::here("output/tf_motif/dars_df_motif.rds"))

# Zelda and HSF motifs
zelda_motif <- "MA1462.1"
hsf_motif <- "MA1458.1"

# Save the ranges where these motifs appear
zelda_motif_ranges <- dat@assays$peaks@motifs@positions[zelda_motif]
hsf_motif_ranges <- dat@assays$peaks@motifs@positions[hsf_motif]

# Set of GRanges object of the pseudobulk dars
dars_granges <- dars_bulk |>
  rownames_to_column("region") |> 
  tibble() |> 
  distinct(region) |> 
  filter(str_detect(region, "mitochondrion", negate = TRUE)) |> 
  separate(region, into = c("seqnames", "start", "end"), sep = "-") |> 
  mutate(strand = "*") |> 
  GenomicRanges::makeGRangesFromDataFrame()

# Find the dars that include the motifs
zelda_target <- subsetByOverlaps(dars_granges, zelda_motif_ranges) |> 
  as_tibble() |> 
  mutate(peaks = paste(seqnames, start, end, sep = "-")) 
hsf_target <- subsetByOverlaps(dars_granges, hsf_motif_ranges) |> 
  as_tibble() |> 
  mutate(peaks = paste(seqnames, start, end, sep = "-")) 

# Motif histograms 
# Get all zelda peaks
zelda_dars <- dars_bulk |> 
  rownames_to_column("region") |> 
  mutate(zelda = ifelse(region %in% zelda_target$peaks, TRUE, FALSE)) |> 
  filter(zelda) |> 
  filter(sig) 

zelda_peaks <- dars_bulk |> 
  rownames_to_column("region") |> 
  mutate(zelda = ifelse(region %in% zelda_target$peaks, TRUE, FALSE)) |> 
  filter(zelda)

# Zelda dars
p_zelda <- dars_bulk |> 
  rownames_to_column("region") |> 
  mutate(zelda = ifelse(region %in% zelda_target$peaks, TRUE, FALSE)) |> 
  filter(zelda) |> 
  filter(sig) |> 
  ggplot() +
  geom_histogram(aes(x = -avg_log2FC),
                 color = "grey20",
                 fill = "grey80",
                 bins = 20) +
  geom_vline(xintercept = 0,
             color = "grey70") +
  scale_y_continuous(expand = c(0, 0),
                     breaks = c(0, 2,4,6,8, 10)) +
  scale_x_continuous(limits = c(-0.7, 0.7),
                     breaks = c(-0.6, -0.3, 0, 0.3, 0.6)) +
  labs(x = "log2(fold-change)",
       y = "Number of DARs") +
  theme_minimal_hgrid()

# HSF DARs
p_hsf <- dars_bulk |>
  rownames_to_column("region") |>
  mutate(hsf = ifelse(region %in% hsf_target$peaks, TRUE, FALSE)) |>
  filter(hsf) |>
  filter(sig) |>
  ggplot() +
  geom_histogram(aes(x = -avg_log2FC),
                 color = "grey20",
                 fill = "grey80",
                 bins = 20) +
  geom_vline(xintercept = 0,
             color = "grey70") +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(limits = c(-0.7, 0.7),
                     breaks = c(-0.6, -0.3, 0, 0.3, 0.6)) +
  labs(x = "log2(fold-difference)",
       y = "Number of DARs") +
  theme_minimal_hgrid()


p_zelda_hsf  <- cowplot::plot_grid(
  NULL,
  p_zelda,
  p_hsf,
  NULL,
  rel_widths = c(1, 1),
  nrow = 1
)



# Need to pick two of these!
zelda_rep_dars <- c(
  "3L-16773087-16773529",
  "2L-1517733-1518344",
  "3L-14944198-14944514",
  "3R-4662778-4663171",
  "3R-28631629-28632049",
  "3L-4408296-4410088"
)

hsf_rep_dars <- c("2R-10486272-10486985",
                 "2L-14409388-14410630")

cov_tracks_to_plot <- c(zelda_rep_dars[1:4], hsf_rep_dars[1:2])
Idents(dat) <- "acc_temp"

p_covs <- list()

# Loop through all six plots

for (i in 1:length(cov_tracks_to_plot)) {
  
  cov_plot <- CoveragePlot(
    object = dat,
    region = cov_tracks_to_plot[i],
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
                          str_extract(cov_tracks_to_plot[i], "-[:digit:]+-"), "-")),
                        xmax = as.numeric(str_remove_all(
                          str_extract(cov_tracks_to_plot[i], "-[:digit:]+$"), "-")),
                        ymin = -Inf,
                        ymax = Inf,
                        fill = "grey75",
                        alpha = 0.5)
  
  cov_plot$layers <- c(peak_rect, cov_plot$layers)
  
  gene_plot <- AnnotationPlot(
    object = dat,
    region = cov_tracks_to_plot[i],
    extend.upstream = 5000,
    extend.downstream = 5000) & 
    theme_cowplot() &
    theme(axis.line.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.y = element_blank(),
          axis.title.y = element_blank())
  
  peak_plot <- PeakPlot(
    object = dat,
    region = cov_tracks_to_plot[i],
    assay = "peaks",
    extend.upstream = 5000,
    extend.downstream = 5000) &
    theme_cowplot() &
    theme(axis.line.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.y = element_blank(),
          axis.title.y = element_blank(),
          legend.position = "none")
  
  p_covs[[i]] <- CombineTracks(
    plotlist = list(
      cov_plot,
      peak_plot,
      gene_plot
    ),
    heights = c(8, 0.1, 2)
  )
}




# Plot all coverage plots
p_coverage_plots <- cowplot::plot_grid(
  p_covs[[1]],
  p_covs[[2]],
  p_covs[[5]],
  p_covs[[3]],
  p_covs[[4]],
  p_covs[[6]],
  nrow = 2
)

f3 <- cowplot::plot_grid(
  p_zelda_hsf,
  p_coverage_plots,
  rel_heights = c(0.25, 1),
  nrow = 2
)

ggsave(here::here(fig_dir, "fig_3.png"),
       plot = f3,
       height = 30,
       width = 40,
       bg = "white",
       units = "cm")
ggsave(here::here(fig_dir, "fig_3.svg"),
       plot = f3,
       height = 30,
       width = 40,
       bg = "white",
       units = "cm")
ggsave(here::here(fig_dir, "fig_3.pdf"),
       plot = f3,
       height = 30,
       width = 40,
       bg = "white",
       units = "cm")

# Figure 4 ---------------------------------------------------------------------

# Load data
dotplot_dat <- read_csv(
  here::here("data/processed/scenic/dotplot.csv")) |> 
  mutate(eRegulon_type = factor(
    str_extract(eRegulon_name, "(?<=_)\\+|-(?=_)"), levels = c("+", "-"))) |> 
  mutate(Region_signature_name = eRegulon_name) |> 
  mutate(eRegulon_name = str_replace_all(eRegulon_name, "_[\\+|-]_", " ")) |> 
  mutate(eRegulon_type = ifelse(eRegulon_type == "+", "activator", "repressor"))

# Create eRegulon order based on TF expression per cell-type
eRegulon_order <- dotplot_dat |> 
  group_by(eRegulon_name, eRegulon_type) |> 
  distinct() |>
  slice_max(color_val, n = 1) |> 
  arrange(eRegulon_type, index) |> 
  pull(eRegulon_name)

# Make the plot for heatmap
p_ereg <- dotplot_dat |> 
  distinct() |> 
  separate(index, into = c("cell_type", "acclimation"), 
           sep = " (?=[[:digit:]])") |>
  mutate(cell_type = ordered(cell_type)) |> 
  mutate(eRegulon_name = factor(eRegulon_name, 
                                levels = eRegulon_order)) |> 
  ggplot(aes(y = interaction(acclimation, cell_type, sep = " "), 
             x = eRegulon_name,
             fill = color_val)) +
  geom_point(aes(size = size_val),
             shape = 21,
             stroke = 0.1) +
  scale_fill_gradient(low = "grey95", 
                      # mid = "grey85", 
                      high = "forestgreen",
                      #midpoint = 0.5,
                      name = "TF exp") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 0.98),
        axis.title = element_blank(),
        strip.background = element_rect(fill = "grey90", color = "grey90")) +
  scale_y_discrete(labels = rep(c("18°C", "25°C"), 9)) +
  facet_grid(rows = vars(cell_type),
             cols = vars(eRegulon_type),
             switch = "y",
             scales = "free",
             labeller = label_wrap_gen(multi_line = TRUE)) +
  ggh4x::force_panelsizes(cols = c(1, 0.4)) +
  scale_size_continuous(range = c(0.25, 6),
                        name = "Target region\nenrichment") +
  theme(strip.placement = "outside",
        strip.text.y.left = element_text(angle = 0, 
                                         hjust = 1, 
                                         margin = margin(r = 5)),
        strip.background.y = element_rect(fill ="grey98", color = "grey98"))

# Target site accessibility 
dars_tf <- eReg |> 
  distinct(TF, Region_signature_name, Region, Gene) |> 
  mutate(Region = str_remove_all(Region, "chr")) |> 
  mutate(Region = str_replace(Region, ":", "-")) |> 
  right_join(dars, by = c("Region" = "region"),
             relationship = "many-to-many")

# GRNs
GRN <- full_join(
  dars_tf |> dplyr::rename("gene" = "Gene"), 
  degs,
  by = c("gene", "cell_type"),
  suffix = c(".region", ".gene")) |> 
  right_join(degs |> 
               filter(sig) |> 
               filter(gene %in% unique(eReg$TF)) |> 
               rename_with(.cols = p_val:sig, function(x){paste0(x, ".TF")}),
             by = c("TF" = "gene", "cell_type")) |> 
  pivot_longer(cols = contains("avg_log2FC"),
               names_to = "group", 
               values_to = "L2FC") |> 
  mutate(group = str_remove_all(group, "avg_log2FC_18_25.")) |> 
  mutate(group = factor(group, levels = c("TF", "region", "gene")))


p_net <- GRN |>
  mutate(cell_type_eGRN = paste(Region_signature_name, cell_type, sep = "\n")) |> 
  mutate(L2FC = -L2FC) |> 
  filter(cell_type_eGRN %in% diffGRN) |> 
  filter(pct.18.gene >= 0.1 | pct.25.gene >= 0.1,
         pct.18.region >= 0.1 | pct.25.region >= 0.1) |>
  ggplot(
    aes(x = group,
        y = L2FC,
        color = L2FC,
        fill = L2FC)) +
  geom_hline(yintercept = 0,
             color = "grey20") +
  geom_line(aes(group = paste(TF, Region, gene)),
            size = 0.2) +
  geom_point(aes(size = group),
             color = "grey80",
             shape = 21) +
  scale_size_manual(values = c(4, 2, 2)) +
  scale_y_continuous(limits = c(-3, 3),
                     name = "log2(fold-difference)") +
  scale_x_discrete(name = element_blank(),
                   labels = c("Target gene exp", 
                              "Motif peak accessibility",
                              "TF exp"),
                   limits = rev) +
  cowplot::theme_minimal_grid() +
  theme(axis.line.x = element_blank(),
        legend.position = "bottom",
        strip.background = element_rect(fill = "grey95"),
        strip.text = element_text(size = 10)) +
  guides(size = FALSE) +
  scale_color_gradient2(low = "blue",
                        mid = "grey90",
                        high = "red",
                        name = element_blank()) +
  scale_fill_gradient2(low = "blue",
                        mid = "grey90",
                        high = "red",
                       name = element_blank()) +
  facet_wrap(~cell_type_eGRN,
             nrow = 2,
             labeller = label_wrap_gen(width = 20)) +
  coord_flip()


f4 <- cowplot::plot_grid(
  p_ereg,
  p_net,
  nrow = 2,
  labels = c("A", "B")
)

ggsave(here::here(fig_dir, "fig_4.png"),
       plot = f4,
       height = 35,
       width = 35,
       bg = "white",
       units = "cm"
)

ggsave(here::here(fig_dir, "fig_4.svg"),
       plot = f4,
       height = 35,
       width = 35,
       bg = "white",
       units = "cm"
)

# Figure 5 ---------------------------------------------------------------------

# Load data
dat <- readRDS(here::here("data/processed/seurat_object/13_dat_wgcna.rds"))
DMEs <- readRDS(here::here("output/wgcna/DMEs.rds"))

DMEs <- DMEs |>
  rownames_to_column("mod") |> 
  mutate(mod = fct_relevel(mod, DMEs |> arrange(avg_log2FC) |> pull(module)))

mod_colors <- c(
  "red", 
  "brown", 
  "blue", 
  "turquoise", 
  "black", 
  "yellow", 
  "green"
) |> colorspace::desaturate(amount = 0.5)

mod_tally <- GetModules(dat) |> 
  group_by(module) |> 
  tally()

p_wgcna <- left_join(DMEs, mod_tally) |> 
  ggplot(aes(y = mod,
             x = -avg_log2FC)) +
  geom_segment(aes(xend = 0, 
                   yend = mod,
                   color = mod),
               linewidth = 1) +
  # annotate(geom = "text", 
  #          y = "red", 
  #          x = -0.18, 
  #          color = "grey20",
  #          label = "*", 
  #          vjust = 0.8, 
  #          size = 10) +
  # annotate(geom = "text", 
  #          y = "brown", 
  #          x = -0.18, 
  #          color = "grey20",
  #          label = "*", 
  #          vjust = 0.8, 
  #          size = 10) +
  geom_point(aes(fill = mod, 
                 size = n), 
             color = "grey80",
             shape = 21) +
  geom_text(aes(label = n),
            color = c("grey90", "grey90", "grey90", "grey30",
                      "grey30", "grey95", "grey20")) +
  geom_vline(xintercept = 0, color = "grey50") +
  scale_x_continuous(name = "log2(fold-change)",
                     limits = c(-0.25, 0.25)) +
  scale_fill_manual(values = rev(mod_colors)) +
  scale_color_manual(values = rev(mod_colors)) +
  scale_size_continuous(range = c(7, 12),
                        name = "Number of\ngenes",
                        breaks = c(75, 100, 250, 500)) +
  labs(y = element_blank()) +
  cowplot::theme_minimal_vgrid() +
  theme(panel.grid.major.y = element_line(color = "grey95"),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        legend.position = "none") 


MEs_umap <- GetMEs(dat) |> 
  rownames_to_column("cell") |> 
  full_join(dat@reductions[["umap"]]@cell.embeddings |> 
              as_tibble(rownames = "cell")) |> 
  as_tibble() |> 
  pivot_longer(cols = green:black,
               values_to = "value",
               names_to = "module")

# Save module colors and desaturate for prettier plotting
mod_colors <- c(
  "turquoise", 
  "brown",
  "yellow", 
  "blue", 
  "green",
  "black", 
  "red"
) |> colorspace::desaturate(amount = 0.5)

# Reset the module colors for plotting
dat <- ResetModuleColors(dat, mod_colors)

# Create all module plots
wgcna_plots <- ModuleFeaturePlot(
  dat,
  features = "MEs",
  order = TRUE)

# Temporarily save all plots as a ggplot object
p <- wgcna_plots |> map(ggplotify::as.ggplot)

# Re-order UMAPs of modules
p_wgcna_umap <- cowplot::plot_grid(
  p$red, 
  p$brown,
  p$blue,
  p$turquoise,
  p$black,
  p$yellow,
  p$green,
  ncol = 2
)

# GO ORA dot plot 
red_go <- read_csv(here::here("output/wgcna/red_go.csv"))
brown_go <- read_csv(here::here("output/wgcna/brown_go.csv"))
green_go <- read_csv(here::here("output/wgcna/green_go.csv"))

go <- bind_rows(list("red" = red_go, 
                     "brown" = brown_go,
                     "green" = green_go), .id = "group") |> 
  mutate(group = factor(group, levels = c("red", "brown", "green")))
  
# GO dot plot
p_go <- go |> 
  group_by(group) |> 
  slice_min(FDR, n = 10) |> 
  mutate(Description = fct_reorder(Description, Ratio)) |> 
  ggplot() +
  geom_point(aes(y = Description,
                 x = Ratio,
                 size = Size,
                 fill = FDR),
             shape = 21,
             color = "grey30") +
  scale_size(range = c(1, 5),
             breaks = c(30, 50, 100, 150),
             name = "Number of\ngenes") +
  scale_fill_gradient(low = "red", high = "grey80") +
  scale_y_discrete(name = element_blank()) +
  scale_x_continuous(name = "Gene Ratio") +
  theme_minimal() +
  facet_wrap(~ group, nrow = 3, scales = "free") +
  ggh4x::force_panelsizes(rows = c(0.3, 1.2, 1.2)) 

p_2_col <- cowplot::plot_grid(
  p_wgcna, 
  p_go, 
  labels = c("B", "C"),
  ncol = 1
)


f5 <- cowplot::plot_grid(
  p_wgcna_umap,
  p_2_col,
  labels = c("A", ""),
  ncol = 2,
  rel_widths = c(1, 0.75)
)

# Save
ggsave(here::here(fig_dir, "fig_5.pdf"),
       plot = f5,
       height= 25,
       width = 40,
       bg = "white",
       units = "cm")
ggsave(here::here(fig_dir, "fig_5.png"),
       plot = f5,
       height= 25,
       width = 40,
       bg = "white",
       units = "cm")
ggsave(here::here(fig_dir, "fig_5.svg"),
       plot = f5,
       height= 25,
       width = 40,
       bg = "white",
       units = "cm")

# Figure 6 ---------------------------------------------------------------------

# Data frame for DARs and DEGS with tally for plotting 
degs_df <- degs |>  
  filter(sig == TRUE) |> 
  group_by(cell_type) |> 
  add_tally() |> 
  full_join(color_cell_type) 

dars_df <- dars |>  
  filter(sig == TRUE) |> 
  group_by(cell_type) |> 
  add_tally() |> 
  full_join(color_cell_type) 

# Define order for plotting
cell_type_order <- degs_df |>
  arrange(n) |> 
  distinct(cell_type) |> 
  pull()
  
# DEGs plot
p_degs <- degs_df |> 
  ggplot(aes(x = -avg_log2FC_18_25,
             y = factor(cell_type, levels = cell_type_order),
             fill = colors,
             color = colors)) +
  ggbeeswarm::geom_beeswarm(
    shape = 21,
    color = "grey60",
    size = 0.7,
    cex = 0.6,
    side = 1,
    alpha = 0.8) +
  geom_vline(xintercept = 0,
             color = "grey30",
             linetype = 5) +
  geom_label(data = degs_df |> distinct(cell_type, n, colors),
             aes(x = -2.75,
                 y = cell_type,
                 label = n),
             color = "grey20",
             vjust = -0.5,
             alpha = 0.4) +
  scale_y_discrete(name = element_blank(), 
                   expand = c(0, 0.05)) +
  scale_x_continuous(name = "Log2 fold-change",
                     limits = c(-3, 3.2),
                     breaks = -3:3) +
  scale_fill_identity() +
  scale_color_identity() +
  ggridges::geom_density_ridges(alpha = 0.5) +
  theme_minimal_hgrid() +
  theme(axis.text.y = element_blank())

# DARs plot
p_dars <- dars_df |> 
  ggplot(aes(x = -avg_log2FC_18_25,
             y = factor(cell_type, levels = cell_type_order),
             fill = colors,
             color = colors)) +
  ggbeeswarm::geom_beeswarm(
    shape = 21,
    color = "grey60",
    size = 1.5,
    cex = 1.05,
    side = 1,
    alpha = 0.8) +
  geom_vline(xintercept = 0,
             color = "grey30",
             linetype = 5) +
  geom_label(data = dars_df |> distinct(cell_type, n, colors) |> 
               mutate(n = ifelse(is.na(n), 0, n)),
             aes(x = -3,
                 y = cell_type,
                 label = n),
             color = "grey20",
             vjust = -0.5,
             alpha = 0.4) +
  scale_y_discrete(name = element_blank(),
                   expand = c(0, 0.05)) +
  scale_x_continuous(name = "Log2 fold-change",
                     limits = c(-3.25, 3.25),
                     breaks = c(-3:3)) +
  scale_fill_identity() +
  scale_color_identity() +
  ggridges::geom_density_ridges(alpha = 0.5) +
  theme_minimal_hgrid()

# Save
ggsave(here::here(fig_dir, "fig_6.svg"),
       plot = p_ridge,
       height = 20,
       width = 25,
       bg = "white",
       units = "cm")
ggsave(here::here(fig_dir, "fig_6.pdf"),
       plot = p_ridge,
       height = 20,
       width = 25,
       bg = "white",
       units = "cm")

# ------------------------------------------------------------------------------
# Supplementary figures --------------------------------------------------------
# ------------------------------------------------------------------------------

# Figure S1 --------------------------------------------------------------------

# Define count cutoffs
low_ATAC <- 800
low_RNA <- 200

# Knee-plot for ATAC 
p_knee_atac <- tibble(nCount_ATAC = sort(dat_raw$nCount_ATAC,
                                         decreasing = TRUE)) |>
  rownames_to_column("rank") |>
  mutate(rank = as.numeric(rank)) |>
  head(100000) |>
  ggplot() +
  geom_point(aes(x = rank,
                 y = nCount_ATAC),
             color = "grey50",
             shape = 21) +
  geom_hline(aes(yintercept = low_ATAC),
             linetype = 2,
             color = "grey20",
             alpha = 0.8,
             linewidth = 1.1) +
  scale_y_continuous(trans = "log10", 
                     breaks = c(100, 1000, 10000, 100000),
                     labels = c("100", "1k", "10k", "100k"),
                     name = "ATAC counts") +
  scale_x_continuous(trans = "log10",
                     breaks = c(1, 10, 100, 1000, 10000, 100000),
                     labels = c("1", "10", "100", "1k", "10k", "100k"),
                     name = "Barcodes in rank order") +
  theme_minimal_grid()

# Knee-plot for RNA 
p_knee_rna <- tibble(nCount_RNA = sort(dat_raw$nCount_RNA,
                                       decreasing = TRUE)) |>
  rownames_to_column("rank") |>
  mutate(rank = as.numeric(rank)) |>
  head(100000) |>
  ggplot() +
  geom_point(aes(x = rank,
                 y = nCount_RNA),
             color = "grey50",
             shape = 21) +
  geom_hline(aes(yintercept = low_RNA),
             linetype = 2,
             color = "grey20",
             alpha = 0.8,
             linewidth = 1.1) +
  scale_y_continuous(trans = "log10", 
                     breaks = c(100, 1000, 10000, 100000),
                     labels = c("100", "1k", "10k", "100k"),
                     name = "RNA counts") +
  scale_x_continuous(trans = "log10",
                     breaks = c(1, 10, 100, 1000, 10000, 100000),
                     labels = c("1", "10", "100", "1k", "10k", "100k"),
                     name = "Barcodes in rank order") +
  theme_minimal_grid()


# Knee plots together
p_knees <- cowplot::plot_grid(
  p_knee_atac, 
  p_knee_rna,
  labels = c("A", "B"),
  nrow = 1
)

# Scatter plot
p_scatter <- dat_raw@meta.data |>
  filter(nCount_RNA >= 1 &
           nCount_ATAC >= 1) |>
  rownames_to_column("cells") |>
  mutate(cells_10x = ifelse(cells %in% Cells(dat) &
                              nCount_RNA >= low_RNA &
                              nCount_ATAC >= low_ATAC, 
                            "Nuclei", 
                            "Non-nuclei barcodes")) |>
  mutate(cells_10x = factor(cells_10x,
                            levels = c("Non-nuclei barcodes",
                                       "Nuclei"))) |>
  arrange(cells_10x) |> 
  ggplot(aes(x = nCount_RNA,
             y = nCount_ATAC)) +
  geom_point(aes(color = cells_10x),
             shape = 21) +
  geom_segment(aes(x = low_RNA, xend = low_RNA,
                   y = low_ATAC, yend = 150000),
               linetype = 2,
               color = "grey50",
               alpha = 0.8,
               linewidth = 1.1) +
  geom_segment(aes(y = low_ATAC, yend = low_ATAC,
                   x = low_RNA, xend = 50000),
               linetype = 2,
               color = "grey50",
               alpha = 0.8,
               linewidth = 1.1) +
  scale_x_continuous(trans = "log10", 
                     breaks = c(1, 10, 100, 1000, 10000, 100000),
                     labels = c("1", "10", "100", "1k", "10k", "100k"),
                     name = "RNA counts") +
  scale_y_continuous(trans = "log10", 
                     breaks = c(1, 10, 100, 1000, 10000, 100000),
                     labels = c("1", "10", "100", "1k", "10k", "100k"),
                     name = "ATAC counts") +
  scale_color_manual(values = c("grey80", "#91AB73"),
                     name = element_blank()) +
  theme_minimal_grid() +
  theme(legend.position = "bottom")

# Number of cells per sample
p_cells_sample <- dat@meta.data |>
  group_by(sample_name) |>
  mutate(sample_name = factor(sample_name, levels = c(c("25C_Rep2",
                                                        "25C_Rep1",
                                                        "18C_Rep2",
                                                        "18C_Rep1")))) |>
  tally() |>
  ggplot(aes(x = n/10^3,
             y = sample_name,
             label = scales::comma_format()(n))) +
  geom_col(color = "grey20",
           fill = "grey80",
           width = 0.8)+
  scale_x_continuous(expand = c(0, 0),
                     name = "Number of nuclei",
                     labels = 
                       scales::label_number(scale = 1, suffix = "k")) +
  geom_label(nudge_x = -0.5,
             label.size = 0.1) +
  scale_y_discrete(name = element_blank(),
                   labels = c("25°C Rep 2",
                              "25°C Rep 1",
                              "18°C Rep 2",
                              "18°C Rep 1")) +
  cowplot::theme_minimal_vgrid() 

# Total number of ATAC reads
p_reads_atac <- dat@meta.data |>
  group_by(sample_name) |>
  mutate(sample_name = factor(sample_name, levels = c(c("25C_Rep2",
                                                        "25C_Rep1",
                                                        "18C_Rep2",
                                                        "18C_Rep1")))) |>
  summarize(n = sum(nCount_ATAC)) |>
  ggplot(aes(x = n/10^6,
             y = sample_name,
             label = scales::label_number(scale = 1, suffix = "M")(round(n/(10^6), digits = 2)))) +
  geom_col(color = "grey20",
           fill = "grey80",
           width = 0.8)+
  scale_x_continuous(expand = c(0, 0),
                     name = "Number of ATAC reads",
                     labels = 
                       scales::label_number(scale = 1, suffix = "M")) +
  geom_label(nudge_x = -3,
             label.size = 0.1) +
  scale_y_discrete(name = element_blank(),
                   labels = c("25°C Rep 2",
                              "25°C Rep 1",
                              "18°C Rep 2",
                              "18°C Rep 1")) +
  cowplot::theme_minimal_vgrid() 

# Total number of ATAC reads
p_reads_rna <- dat@meta.data |>
  group_by(sample_name) |>
  mutate(sample_name = factor(sample_name, levels = c(c("25C_Rep2",
                                                        "25C_Rep1",
                                                        "18C_Rep2",
                                                        "18C_Rep1")))) |>
  summarize(n = sum(nCount_RNA)) |>
  ggplot(aes(x = n/10^6,
             y = sample_name,
             label = scales::label_number(scale = 1, suffix = "M")(round(n/(10^6), digits = 2)))) +
  geom_col(color = "grey20",
           fill = "grey80",
           width = 0.8)+
  scale_x_continuous(expand = c(0, 0),
                     name = "Number of RNA reads",
                     labels = 
                       scales::label_number(scale = 1, suffix = "M")) +
  geom_label(nudge_x = -0.5,
             label.size = 0.1) +
  scale_y_discrete(name = element_blank(),
                   labels = c("25°C Rep 2",
                              "25°C Rep 1",
                              "18°C Rep 2",
                              "18°C Rep 1")) +
  cowplot::theme_minimal_vgrid() 

# Per sample metrics
p_samples <- cowplot::plot_grid(
  p_cells_sample,
  p_reads_atac,
  p_reads_rna,
  labels = c("D", "E", "F"),
  vjust = 1,
  nrow = 3
)

# Second row of plots
p_row_2 <- cowplot::plot_grid(
  p_scatter,
  p_samples,
  labels = c("C", ""),
  rel_widths = c(1.3, 1),
  nrow = 1
)

# All plots together
fs1 <- cowplot::plot_grid(
  p_knees, 
  p_row_2,
  rel_heights = c(1, 1.3),
  nrow = 2
)

# Save
ggsave(here::here(fig_dir, "fig_s1.png"),
       plot = fs1,
       height = 25,
       width = 25,
       units = "cm")
ggsave(here::here(fig_dir, "fig_s1.pdf"),
       plot = fs1,
       height = 25,
       width = 25,
       units = "cm")


# Figure S2 --------------------------------------------------------------------

# Count ATAC
p_c_atac <- dat@meta.data |>
  ggplot(aes(x = nCount_ATAC,
             y = sample_name_fig)) +
  geom_violin(color = "grey20",
              fill = "grey80",
              width = 0.8) +
  geom_boxplot(color = "grey20",
              fill = "grey90",
              width = 0.1,
              outlier.shape = NA) +
  scale_y_discrete(name = element_blank()) +
  scale_x_continuous(expand = c(0, 0),
                     labels = scales::comma_format(),
                     name = "ATAC reads per cell",
                     limits = c(0, 22000),
                     position = "top") +
  cowplot::theme_minimal_grid() 

# Feature ATAC
p_f_atac <- dat@meta.data |>
  ggplot(aes(x = nFeature_ATAC,
             y = sample_name_fig)) +
  geom_violin(color = "grey20",
              fill = "grey80",
              width = 0.8) +
  geom_boxplot(color = "grey20",
               fill = "grey90",
               width = 0.1,
               outlier.shape = NA) +
  scale_y_discrete(name = element_blank()) +
  scale_x_continuous(expand = c(0, 0),
                     labels = scales::comma_format(),
                     name = "ATAC features per cell",
                     limits = c(0, 6500),
                     position = "top") +
  cowplot::theme_minimal_grid() 

# Count RNA
p_c_rna <- dat@meta.data |>
  ggplot(aes(x = nCount_RNA,
             y = sample_name_fig)) +
  geom_violin(color = "grey20",
              fill = "grey80",
              width = 0.8) +
  geom_boxplot(color = "grey20",
               fill = "grey90",
               width = 0.1,
               outlier.shape = NA) +
  scale_y_discrete(name = element_blank()) +
  scale_x_continuous(expand = c(0, 0),
                     labels = scales::comma_format(),
                     name = "RNA reads per cell",
                     limits = c(0, 7000),
                     position = "top") +
  cowplot::theme_minimal_grid() 

# Feature RNA
p_f_rna <- dat@meta.data |>
  ggplot(aes(x = nFeature_RNA,
             y = sample_name_fig)) +
  geom_violin(color = "grey20",
              fill = "grey80",
              width = 0.8) +
  geom_boxplot(color = "grey20",
               fill = "grey90",
               width = 0.1,
               outlier.shape = NA) +
  scale_y_discrete(name = element_blank()) +
  scale_x_continuous(expand = c(0, 0),
                     labels = scales::comma_format(),
                     name = "Expressed genes detected per cell",
                     limits = c(0, 2250),
                     position = "top") +
  cowplot::theme_minimal_grid() 

# Count and feature plots together
p_c_f <- cowplot::plot_grid(
  p_c_rna,
  p_c_atac,
  p_f_rna,
  p_f_atac,
  nrow = 2,
  labels = c("A", "C", "B", "D")
)


# Plot transcription start site enrichment for each sample
p_tss <- TSSPlot(
  dat,
  assay = "ATAC",
  group.by = "sample_name_fig") +
  scale_color_manual(values = rep("grey50", 4)) +
  cowplot::theme_minimal_grid() +
  theme(legend.position = "none",
        plot.title = element_blank(),
        strip.background = element_rect(fill = "grey90"))

# Detailed RNA QC metrics

# Percent mito
p_mt <- dat@meta.data |>
  ggplot(aes(x = percent.mt/100,
             y = sample_name_fig)) +
  geom_violin(color = "grey20",
              fill = "grey80",
              width = 0.8) +
  geom_boxplot(color = "grey20",
               fill = "grey90",
               width = 0.1,
               outlier.shape = NA) +
  scale_y_discrete(name = element_blank()) +
  scale_x_continuous(expand = c(0, 0),
                     name = "RNA reads mapped to mitochondrial genome",
                     limits = c(0, .0525),
                     labels = scales::label_percent(),
                     position = "top") +
  cowplot::theme_minimal_grid() 

# Percent ribo
p_ribo <- dat@meta.data |>
  ggplot(aes(x = percent.ribo/100,
             y = sample_name_fig)) +
  geom_violin(color = "grey20",
              fill = "grey80",
              width = 0.8) +
  geom_boxplot(color = "grey20",
               fill = "grey90",
               width = 0.1,
               outlier.shape = NA) +
  scale_y_discrete(name = element_blank()) +
  scale_x_continuous(expand = c(0, 0),
                     name = "RNA reads from ribosomal genes",
                     limits = c(0, 0.27),
                     labels = scales::label_percent(),
                     position = "top") +
  cowplot::theme_minimal_grid() 

# Count and feature plots together
p_rna_qc <- cowplot::plot_grid(
  p_mt,
  p_ribo,
  nrow = 2,
  labels = c("F", "G"),
  vjust = 1
)

# Detailed ATAC QC metrics

# Nucleosome signal
p_nuc <- dat@meta.data |>
  ggplot(aes(x = nucleosome_signal,
             y = sample_name_fig)) +
  geom_violin(color = "grey20",
              fill = "grey80",
              width = 0.8) +
  geom_boxplot(color = "grey20",
               fill = "grey90",
               width = 0.1,
               outlier.shape = NA) +
  scale_y_discrete(name = element_blank()) +
  scale_x_continuous(expand = c(0, 0),
                     name = "Nucleosome signal",
                     limits = c(0, 1.4),
                     breaks = seq(from = 0, to = 1.25, by = 0.25),
                     position = "top") +
  cowplot::theme_minimal_grid() 

# FRiP 
p_frip <- dat@meta.data |>
  ggplot(aes(x = FRiP,
             y = sample_name_fig)) +
  geom_violin(color = "grey20",
              fill = "grey80",
              width = 0.8) +
  geom_boxplot(color = "grey20",
               fill = "grey90",
               width = 0.1,
               outlier.shape = NA) +
  scale_y_discrete(name = element_blank()) +
  scale_x_continuous(expand = c(0, 0),
                     name = "Fraction of reads in peaks (FRiP)",
                     limits = c(0, 1.05),
                     breaks = seq(from = 0, to = 1, by = 0.25),
                     position = "top") +
  cowplot::theme_minimal_grid() 


# Count and feature plots together
p_atac_qc <- cowplot::plot_grid(
  p_nuc,
  p_frip,
  nrow = 2,
  labels = c("H", "I"),
  vjust = 1
)

p_qc <- cowplot::plot_grid(
  p_rna_qc,
  p_atac_qc,
  nrow = 1
)

# Combine all plots
cowplot::plot_grid(
  p_c_f,
  p_tss,
  p_qc,
  nrow = 3,
  labels = c(NA, "E", NA)
)

# Save
ggsave(here::here(fig_dir, "fig_s2.png"),
       height = 30,
       width = 25,
       units = "cm")
ggsave(here::here(fig_dir, "fig_s2.pdf"),
       height = 30,
       width = 25,
       units = "cm")

# Figure S3 --------------------------------------------------------------------


# ATAC-only dim reduction 
dat_atac <- readRDS(
  here::here("data/processed/seurat_object/07_dat_atac_only_dim.rds")
)
# RNA-only dim reduction
dat_rna <- readRDS(
  here::here("data/processed/seurat_object/07_dat_rna_only_dim.rds")
)

# Look at the different clustering resolutions
library(clustree)
p_clust <- clustree(
  dat,
  prefix = "wknn_res.") +
  theme(legend.key.size = unit(8, "pt"))

# RNA-only UMAP
p_rna_umap <- DimPlot(
  dat_rna,
  group.by = "wknn_res.0.7",
  reduction = "umap") +
  labs(title = element_blank()) +
  scale_color_manual(
    values = scales::hue_pal()(12) |> 
      colorspace::desaturate(amount = 0.5)
  ) +
  theme_void() +
  theme(legend.position = "none") +
  guides(color = guide_legend(byrow = TRUE,
                              nrow = 1,
                              override.aes = list(size = 3)))

p_rna_umap_labeled <- LabelClusters(
  p_rna_umap, 
  id = "wknn_res.0.7",
  fontface = "bold",
  color = "grey20",
  alpha = 0.8,
  box = TRUE,
  repel = TRUE
)

# Acclimation temperature conditions plotted on top of each other
p_atac_umap <- DimPlot(
  dat_atac,
  group.by = "wknn_res.0.7",
  reduction = "umap") +
  labs(title = element_blank()) +
  scale_color_manual(
    values = scales::hue_pal()(12) |> 
      colorspace::desaturate(amount = 0.5)
  ) +
  theme_void() +
  theme(legend.position = "none") +
  guides(color = guide_legend(byrow = TRUE,
                              nrow = 1,
                              override.aes = list(size = 3)))

p_atac_umap_labeled <- LabelClusters(
  p_atac_umap, 
  id = "wknn_res.0.7",
  fontface = "bold",
  color = "grey20",
  alpha = 0.8,
  box = TRUE,
  repel = TRUE
)


p_umaps <- cowplot::plot_grid(
  p_atac_umap_labeled,
  p_rna_umap_labeled,
  labels = c("B", "C"),
  nrow = 2
)

fs3 <- cowplot::plot_grid(
  p_clust,
  p_umaps,
  rel_widths = c(1.2, 1),
  labels = c("A"),
  nrow = 1
)

# Save
ggsave(here::here(fig_dir, "fig_s3.png"),
       fs3,
       height = 20,
       width = 22.5,
       units = "cm")
ggsave(here::here(fig_dir, "fig_s3.pdf"),
       plot = fs3,
       height = 20,
       width = 22.5,
       units = "cm")

# Figure S4 --------------------------------------------------------------------

link_gene_peak <- linked_marker_peaks |> 
  select(gene, peak)

atac_hmap_data <- atac_min_max_norm |> 
  as.data.frame() |> 
  rownames_to_column("peak") |> 
  as_tibble()

exp_hmap_data <- dat_agg@assays$SCT$scale.data |> 
  as.data.frame() |> 
  rownames_to_column("gene") |> 
  as_tibble()

corr_links_df <- link_gene_peak |> 
  left_join(atac_hmap_data, 
            by = "peak")

corr_links_df <- corr_links_df |>   
  left_join(exp_hmap_data, 
            by = "gene",
            suffix = c("_ATAC", "_RNA")) |> 
  pivot_longer(contains("_ATAC"), 
               names_to = "cell_type_ATAC",
               values_to = "ATAC") |> 
  pivot_longer(contains("_RNA"), 
               names_to = "cell_type_RNA",
               values_to = "RNA") |> 
  mutate(cell_type_ATAC = str_remove_all(cell_type_ATAC, "_ATAC")) |> 
  mutate(cell_type_RNA = str_remove_all(cell_type_RNA, "_RNA")) |> 
  filter(cell_type_ATAC == cell_type_RNA) |> 
  mutate(cell_type = cell_type_RNA) |> 
  select(gene, peak, cell_type, ATAC, RNA)

fs4 <- corr_links_df |> 
  ggplot(aes(x = ATAC,
             y = RNA)) +
  geom_point(alpha = 0.2) +
  geom_smooth(method = "lm",
              color = "black",
              se = FALSE,
              linetype = 2) +
  ggpubr::stat_cor(aes(label = after_stat(r.label))) +
  facet_wrap(~cell_type) +
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1),
                     labels = c("0", "0.25", "0.5", "0.75", "1"),
                     name = "ATAC-Seq (normalized min to max)") +
  scale_y_continuous(name = "RNA-Seq (centered and scaled)") +
  theme_minimal() +
  theme(strip.background = element_rect(fill = "grey90",
                                        color = "grey50"))

# Save
ggsave(here::here(fig_dir, "fig_s4.png"),
       fs4,
       height = 20,
       width = 20,
       units = "cm")
ggsave(here::here(fig_dir, "fig_s4.pdf"),
       plot = fs4,
       height = 20,
       width = 20,
       units = "cm")


# Figure S5 --------------------------------------------------------------------

# Marker gene dot plot
markers <- read_csv(here::here("output/markers/cluster_markers_all.csv"))

# Top three markers for each cluster
top_markers <- markers |>
  group_by(cluster) |>
  slice_min(max_pval, n = 3, with_ties = FALSE) |>
  mutate(cluster = factor(cluster, levels = as.character(0:11))) |>
  arrange(cluster)

# Set default assay for fetching data
DefaultAssay(dat) <- "SCT"
Idents(dat) <- "seurat_clusters"

# Percent nuclei with expressed gene
pct_exp <- FetchData(dat,
                     vars = c(top_markers$gene, "seurat_clusters"),
                     layer = "counts") |>
  pivot_longer(cols = -c("seurat_clusters"),
               names_to = "gene",
               values_to = "count") |>
  mutate(exp = count > 0) |>
  group_by(seurat_clusters, gene) |>
  summarize(pct_exp = sum(exp)/n()*100)

# Average expression
avg_exp <- AverageExpression(dat)$SCT |>
  as.data.frame() |>
  rownames_to_column("gene") |>
  filter(gene %in% top_markers$gene) |>
  pivot_longer(cols = -c(gene),
               names_to = "cluster",
               values_to = "expression") |>
  mutate(cluster = str_remove_all(cluster, "g")) |>
  group_by(gene) |>
  mutate(norm_expression = scale(expression))

# Combine the average expression with percent expressed data
df <- full_join(
  avg_exp,
  pct_exp,
  by = c("cluster" = "seurat_clusters", "gene")
)

# Get a pleasing gene order for plotting
gene_order <- df |>
  mutate(cluster = as.numeric(cluster)) |>
  group_by(cluster) |>
  arrange(cluster, desc(norm_expression)) |>
  slice_max(norm_expression, n = 3) |>
  pull(gene)

# Marker gene dot plot
p_mark_dot <- df |>
  mutate(cluster = factor(cluster, levels = as.character(0:11))) |>
  mutate(gene = factor(gene, levels = top_markers$gene)) |>
  ggplot(aes(x = gene,
             y = cluster,
             fill = norm_expression,
             size = pct_exp)) +
  geom_point(shape = 21, stroke = 0.25) +
  scale_fill_gradient(low = "grey95",
                      high = "firebrick",
                      name = "Normalized\nExpression") +
  theme_minimal(base_size = 18) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.98),
        axis.title = element_blank(),
        strip.background = element_rect(fill = "grey90", color = "grey90")) +
  ggh4x::force_panelsizes(cols = c(1, 0.4)) +
  scale_size_continuous(range = c(0.25, 12),
                        labels = scales::label_percent(scale = 1),
                        name = "Nuclei\nexpressing\nthe gene")

# Select genes for individual violin plots
markers_plot <- c(
  "jeb",
  "salm",
  "cnc",
  "scrt",
  "wnd",
  "DNaseII",
  "sns",
  "btl",
  "Blimp-1",
  "otp",
  "srp",
  "meru"
)

# Make sure plotting the correct cluster resolution
Idents(dat) <- "wknn_res.0.7"

# Create all violin plots
p_vlns <- VlnPlot(dat,
        assay = "RNA",
        split.by = "acc_temp",
        features = markers_plot) & 
  scale_fill_manual(values = acc_colors) &
  theme(axis.title.x = element_blank())

fs5 <- cowplot::plot_grid(
  p_mark_dot,
  p_vlns,
  labels = c("A", "B"),
  label_size = 22,
  nrow = 2,
  rel_heights = c(0.75, 1)
)

# Save
ggsave(here::here(fig_dir, "fig_s5.png"),
       fs5,
       height = 40,
       width = 60,
       units = "cm")
ggsave(here::here(fig_dir, "fig_s5.pdf"),
       plot = fs5,
       height = 40,
       width = 60,
       units = "cm")

# Figure S6 --------------------------------------------------------------------

# Set cell types as identity to plot across
Idents(dat) <- "cell_type"
DefaultAssay(dat) <- "peaks"

# Btl coverage plot

# Set gene for coverage plot
gene <- "btl"

# Accessibility peaks
p_cov <- CoveragePlot(
  object = dat,
  region = gene,
  features = gene,
  links = FALSE,
  peaks = FALSE,
  annotation = FALSE,
  extend.upstream = 0,
  extend.downstream = 0) &
  scale_fill_manual(values = cell_type_colors) &
  theme(legend.position = "none",
        strip.background = element_blank(),
        strip.text = element_text(size = 6),
        axis.line.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank())

# Expression plot
p_exp <- ExpressionPlot(
  dat,
  assay = "RNA",
  features = gene) &
  scale_fill_manual(values = cell_type_colors) &
  theme(axis.title.y = element_blank(),
        axis.text.x = element_text(size = 6))

# Gene annotation plot with boxes for gene
p_ann <- AnnotationPlot(
  dat,
  region = gene,
  #size_gene_label = 3,
  extend.upstream = 0,
  extend.downstream = 0)  &
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.line.y = element_blank())

# Combine together
p_btl <- CombineTracks(
  list(p_cov, p_ann),
  expression.plot = p_exp,
  widths = c(8, 2),
  heights = c(5, 1)) &
  theme(axis.text = element_text(size = 6),
        axis.title = element_text(size = 8))

# Mef2 coverage plot 

# Set gene for coverage plot
gene <- "Mef2"

# Accessibility peaks
p_cov <- CoveragePlot(
  object = dat,
  region = gene,
  features = gene,
  links = FALSE,
  peaks = FALSE,
  annotation = FALSE,
  extend.upstream = 0,
  extend.downstream = 0) &
  scale_fill_manual(values = cell_type_colors) &
  theme(legend.position = "none",
        strip.text = element_blank(),
        strip.text.y.left = element_blank(),
        strip.background = element_blank(),
        axis.line.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank())

# Expression plot
p_exp <- ExpressionPlot(
  dat,
  assay = "RNA",
  features = gene) &
  scale_fill_manual(values = cell_type_colors) &
  theme(axis.title.y = element_blank(),
        axis.text.x = element_text(size = 6))

# Gene annotation plot with boxes for gene
p_ann <- AnnotationPlot(
  dat,
  region = gene,
  extend.upstream = 0,
  extend.downstream = 0)  &
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.line.y = element_blank())

# Combine together
p_Mef2 <- CombineTracks(
  list(p_cov, p_ann),
  expression.plot = p_exp,
  widths = c(8, 2),
  heights = c(5, 1)) &
  theme(axis.text = element_text(size = 6),
        axis.title = element_text(size = 8))

# otp coverage plot

# Set gene for coverage plot
gene <- "otp"

# Accessibility peaks
p_cov <- CoveragePlot(
  object = dat,
  region = gene,
  features = gene,
  links = FALSE,
  peaks = FALSE,
  annotation = FALSE,
  extend.upstream = 0,
  extend.downstream = 0) &
  scale_fill_manual(values = cell_type_colors) &
  theme(legend.position = "none",
        strip.text = element_blank(),
        strip.text.y.left = element_blank(),
        strip.background = element_blank(),
        axis.line.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank())

# Expression plot
p_exp <- ExpressionPlot(
  dat,
  assay = "RNA",
  features = gene) &
  scale_fill_manual(values = cell_type_colors) &
  theme(axis.title.y = element_blank(),
        axis.text.x = element_text(size = 6))

# Gene annotation plot with boxes for gene
p_ann <- AnnotationPlot(
  dat,
  region = gene,
  #size_gene_label = 3,
  extend.upstream = 0,
  extend.downstream = 0)  &
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.line.y = element_blank())

# Combine together
p_otp <- CombineTracks(
  list(p_cov, p_ann),
  expression.plot = p_exp,
  widths = c(8, 2),
  heights = c(5, 1)) &
  theme(axis.text = element_text(size = 6),
        axis.title = element_text(size = 8))

# Combine the three tracks
fs5 <- cowplot::plot_grid(
  p_btl,
  p_Mef2,
  p_otp,
  rel_widths = c(1, 0.7, 0.7),
  nrow = 1
)

# Save
ggsave(here::here(fig_dir, "fig_s6.png"),
       plot = fs5,
       height = 10,
       width = 20,
       units = "cm")
ggsave(here::here(fig_dir, "fig_s6.pdf"),
       plot = fs5,
       height = 10,
       width = 20,
       units = "cm")


# Figure S7 --------------------------------------------------------------------

# Bulk volcano plots

# DARs bulk volcano plot
p_dars_volcano <- dars_bulk |> 
  ggplot(aes(y = -log10(p_val_adj), 
             x = -avg_log2FC)) +
  geom_point(aes(fill =  sig),
             color = "grey50",
             shape = 21) +
  geom_hline(aes(yintercept = -log10(0.05)),
             color = "grey20",
             linetype = 2) +
  geom_vline(aes(xintercept = -0.25),
             color = "grey20",
             linetype = 2) +
  geom_vline(aes(xintercept = 0.25),
             color = "grey20",
             linetype = 2) +
  geom_vline(aes(xintercept = 0)) +
  scale_fill_manual(values = c("grey80", "firebrick")) +  
  scale_x_continuous(name = "log2(fold-difference)") +
  scale_y_continuous(name = "-log10(adj. p-value)") +
  guides(fill = FALSE) +
  theme_minimal_grid()

# DEGs bulk volcano plot
p_degs_volcano <- degs_bulk |> 
  as.data.frame() |> 
  rownames_to_column("gene") |> 
  as_tibble() |> 
  ggplot(aes(y = -log10(padj), 
             x = -log2FoldChange)) +
  geom_point(aes(fill = padj < 0.05),
             color = "grey50",
             shape = 21) +
  geom_hline(aes(yintercept = -log10(0.05)),
             color = "grey20",
             linetype = 2) +
  geom_vline(aes(xintercept = 0)) +
  scale_x_continuous(name = "log2(fold-difference)") +
  scale_y_continuous(name = "-log10(adj. p-value)") +
  scale_fill_manual(values = c("grey80", "firebrick")) +  
  guides(fill = FALSE) +
  theme_minimal_grid()


# Both volcano plots
p_volcanos <- cowplot::plot_grid(
  p_dars_volcano, 
  p_degs_volcano,
  labels = c("A", "B"),
  nrow = 1
)


zelda_go <- read_csv(here::here("output/zelda/zelda_go.csv"))

# GO dot plot
p_zelda_go <- zelda_go |> 
  mutate(Description = fct_reorder(Description, Ratio)) |> 
  ggplot() +
  geom_point(aes(y = Description,
                 x = Ratio,
                 size = Size,
                 fill = FDR),
             shape = 21,
             color = "grey30") +
  scale_size(range = c(1, 5),
             breaks = c(30, 50, 100, 150),
             name = "Number of\ngenes") +
  scale_fill_gradient(low = "red", high = "grey80") +
  scale_y_discrete(name = element_blank()) +
  scale_x_continuous(name = "Gene Ratio") +
  theme_minimal()

p_zelda_go_space <- cowplot::plot_grid(
  NULL,
  p_zelda_go,
  labels = c("C", ""),
  rel_widths = c(0.1, 2),
  nrow = 1
)

# Put panels together
fs7 <- cowplot::plot_grid(
  p_volcanos,
  p_zelda_go_space,
  rel_heights = c(0.75, 2),
  nrow = 2
)
# Save
ggsave(here::here(fig_dir, "fig_s7.png"),
       plot = fs7,
       height = 30,
       width = 20,
       units = "cm")
ggsave(here::here(fig_dir, "fig_s7.pdf"),
       plot = fs7,
       height = 30,
       width = 20,
       units = "cm")




# Figure S8 --------------------------------------------------------------------

# See SCENIC+ output in ~/src/04_scenic_plus/00_scenic_plus.ipynb


# Figure S9 --------------------------------------------------------------------


svg(here::here(fig_dir, "fig_s7_dendrogram.svg"))
PlotDendrogram(
  dat, 
  main = "hdWGCNA Dendrogram"
)
dev.off()

# Cell-type specific mean module eigengene score
fs9b_only <- MEs |> 
  full_join(dat@meta.data |> rownames_to_column("cell") |> select(cell, cell_type)) |> 
  group_by(module, cell_type) |> 
  summarize(avg_mod = mean(mod_eigengene)) |> 
  mutate(module = factor(module, 
                         levels = c("red", "brown", "blue", 
                                    "black", "green", "yellow", "turquoise"))) |> 
  full_join(color_cell_type) |> 
  ggplot(aes(y = forcats::fct_rev(cell_type), 
             x = avg_mod,
             alpha = avg_mod > 0,
             fill = colors)) +
  geom_col(color = "grey50") +
  geom_vline(xintercept = 0, 
             color = "grey20") +
  scale_x_continuous(limits = c(-9.1, 9.1),
                     breaks = c(-8, -4, 0 , 4, 8),
                     name = "mean module eigengene score") +
  scale_fill_identity() +
  facet_wrap(~ module, nrow = 2) +
  scale_alpha_discrete(range = c(0.3, 1)) +
  cowplot::theme_minimal_grid() +
  theme(legend.position = "none",
        axis.title.y = element_blank())

# Cell-type specificity index
fs9c_only <- MEs |> 
  full_join(dat@meta.data |> rownames_to_column("cell") |> select(cell, cell_type)) |> 
  group_by(module, cell_type) |> 
  summarize(avg_mod = mean(mod_eigengene)) |> 
  ungroup(cell_type) |> 
  summarize(sd = sd(avg_mod)) |> 
  arrange(sd) |> 
  mutate(module = fct_reorder(module, sd, .desc = TRUE)) |> 
  ggplot(aes(y = module, 
             x = sd)) +
  geom_col(color = "grey20",
           fill = "grey80", 
           width = 0.8) +
  geom_label(aes(label = sprintf("%.2f",round(sd, 2))),
             nudge_x = -0.25,
             size = 3) +
  scale_x_continuous(expand = c(0,0.01),
                     name = "cell-type specificity index") +
  cowplot::theme_minimal_vgrid() +
  theme(axis.title.y = element_blank())


fs9b <- cowplot::plot_grid(
  plotlist = list(NULL, fs9b_only),
  nrow = 1,
  rel_widths = c(0.05, 1)
)


fs9c <- cowplot::plot_grid(
  plotlist = list(NULL, fs9c_only, NULL),
  nrow = 1,
  rel_widths = c(0.2, 1, 0.2)
)


fs9 <- cowplot::plot_grid(
  fs9b,
  fs9c,
  nrow = 2,
  rel_heights = c(1, 0.5)
)

ggsave(here::here(fig_dir, "fig_s9.png"),
       fs9,
       bg = "white",
       height = 20,
       width = 25,
       units = "cm")
ggsave(here::here(fig_dir, "fig_s9.pdf"),
       fs9,
       bg = "white",
       height = 20,
       width = 25,
       units = "cm")
ggsave(here::here(fig_dir, "fig_s9.svg"),
       fs9,
       bg = "white",
       height = 20,
       width = 25,
       units = "cm")



# Figure S10 -------------------------------------------------------------------

# Attach all other go 
blue_go <- read_csv(here::here("output/wgcna/blue_go.csv"))
turquoise_go <- read_csv(here::here("output/wgcna/turquoise_go.csv"))
black_go <- read_csv(here::here("output/wgcna/black_go.csv"))
yellow_go <- read_csv(here::here("output/wgcna/yellow_go.csv"))

other_go <- bind_rows(list("blue" = blue_go,
                           "turquoise" = turquoise_go,
                           "black" = black_go,
                           "yellow" = yellow_go), .id = "group") |> 
  mutate(group = factor(group, levels = c("blue", "turquoise", 
                                          "black", "yellow")))

all_go <- bind_rows(go, other_go) |>
  mutate(group = factor(group, levels = c("red", "brown",
                                          "blue", "turquoise", 
                                          "black", "yellow", "green")))


p_go_red <- all_go |>
  filter(group %in% c("red")) |>
  arrange(Ratio) |>
  mutate(Description = fct_reorder(Description, Ratio)) |>
  ggplot() +
  geom_point(aes(y = Description,
                 x = Ratio,
                 size = Size,
                 fill = FDR),
             shape = 21,
             color = "grey30") +
  scale_size(range = c(1, 5),
             breaks = c(30, 50, 100, 150),
             name = "Number of\ngenes") +
  scale_fill_gradient(low = "red", high = "grey80", limits = c(0, 0.1)) +
  scale_y_discrete(name = element_blank()) +
  scale_x_continuous(name = "Gene Ratio") +
  theme_minimal() +
  facet_wrap(~ group, ncol = 1, scales = "free") +
  theme(legend.position = "none")

p_go_brown <- all_go |>
  filter(group %in% c("brown")) |>
  arrange(Ratio) |>
  mutate(Description = fct_reorder(Description, Ratio)) |>
  ggplot() +
  geom_point(aes(y = Description,
                 x = Ratio,
                 size = Size,
                 fill = FDR),
             shape = 21,
             color = "grey30") +
  scale_size(range = c(1, 5),
             breaks = c(30, 50, 100, 150),
             name = "Number of\ngenes") +
  scale_fill_gradient(low = "red", high = "grey80", limits = c(0, 0.1)) +
  scale_y_discrete(name = element_blank()) +
  scale_x_continuous(name = "Gene Ratio") +
  theme_minimal() +
  facet_wrap(~ group, ncol = 1, scales = "free") + 
  theme(legend.position = "none")

p_go_blue <- all_go |>
  filter(group %in% c("blue")) |>
  arrange(Ratio) |>
  mutate(Description = fct_reorder(Description, Ratio)) |>
  ggplot() +
  geom_point(aes(y = Description,
                 x = Ratio,
                 size = Size,
                 fill = FDR),
             shape = 21,
             color = "grey30") +
  scale_size(range = c(1, 5),
             breaks = c(30, 50, 100, 150),
             name = "Number of\ngenes") +
  scale_fill_gradient(low = "red", high = "grey80", limits = c(0, 0.1)) +
  scale_y_discrete(name = element_blank()) +
  scale_x_continuous(name = "Gene Ratio") +
  theme_minimal() +
  facet_wrap(~ group, ncol = 1, scales = "free") + 
  theme(legend.position = "none")

p_go_yellow <- all_go |>
  filter(group %in% c("yellow")) |>
  arrange(Ratio) |>
  mutate(Description = fct_reorder(Description, Ratio)) |>
  ggplot() +
  geom_point(aes(y = Description,
                 x = Ratio,
                 size = Size,
                 fill = FDR),
             shape = 21,
             color = "grey30") +
  scale_size(range = c(1, 5),
             breaks = c(30, 50, 100, 150),
             name = "Number of\ngenes") +
  scale_fill_gradient(low = "red", high = "grey80", limits = c(0, 0.1)) +
  scale_y_discrete(name = element_blank()) +
  scale_x_continuous(name = "Gene Ratio") +
  theme_minimal() +
  facet_wrap(~ group, ncol = 1, scales = "free") + 
  theme(legend.position = "none")

p_go_turquoise <- all_go |>
  filter(group %in% c("turquoise")) |>
  arrange(Ratio) |>
  mutate(Description = fct_reorder(Description, Ratio)) |>
  ggplot() +
  geom_point(aes(y = Description,
                 x = Ratio,
                 size = Size,
                 fill = FDR),
             shape = 21,
             color = "grey30") +
  scale_size(range = c(1, 5),
             breaks = c(30, 50, 100, 150),
             name = "Number of\ngenes") +
  scale_fill_gradient(low = "red", high = "grey80", limits = c(0, 0.1)) +
  scale_y_discrete(name = element_blank()) +
  scale_x_continuous(name = "Gene Ratio") +
  theme_minimal() +
  facet_wrap(~ group, ncol = 1, scales = "free") + 
  theme(legend.position = "none")

p_go_black <- all_go |>
  filter(group %in% c("black")) |>
  arrange(Ratio) |>
  mutate(Description = fct_reorder(Description, Ratio)) |>
  ggplot() +
  geom_point(aes(y = Description,
                 x = Ratio,
                 size = Size,
                 fill = FDR),
             shape = 21,
             color = "grey30") +
  scale_size(range = c(1, 5),
             breaks = c(30, 50, 100, 150),
             name = "Number of\ngenes") +
  scale_fill_gradient(low = "red", high = "grey80", limits = c(0, 0.1)) +
  scale_y_discrete(name = element_blank()) +
  scale_x_continuous(name = "Gene Ratio") +
  theme_minimal() +
  facet_wrap(~ group, ncol = 1, scales = "free") + 
  theme(legend.position = "none")

p_go_green <- all_go |>
  filter(group %in% c("green")) |>
  arrange(Ratio) |>
  mutate(Description = fct_reorder(Description, Ratio)) |>
  ggplot() +
  geom_point(aes(y = Description,
                 x = Ratio,
                 size = Size,
                 fill = FDR),
             shape = 21,
             color = "grey30") +
  scale_size(range = c(1, 5),
             breaks = c(30, 50, 100, 150),
             name = "Number of\ngenes") +
  scale_fill_gradient(low = "red", high = "grey80", limits = c(0, 0.1)) +
  scale_y_discrete(name = element_blank()) +
  scale_x_continuous(name = "Gene Ratio") +
  theme_minimal() +
  facet_wrap(~ group, ncol = 1, scales = "free")

p_go_left <- ggpubr::ggarrange(p_go_red, 
                               p_go_brown, 
                               p_go_blue, 
                               p_go_yellow,
                               ncol = 1,
                               heights = c(6, 17, 37, 16),
                               align = "v")

p_go_right <- ggpubr::ggarrange(p_go_turquoise, 
                                p_go_black, 
                                p_go_green, 
                                ncol = 1,
                                heights = c(37, 13, 22),
                                align = "v")

p_go_full <- cowplot::plot_grid(NULL, p_go_left, NULL, p_go_right,
                                rel_widths = c(0.05, 1, 0.05, 1.1),
                                ncol = 4)


ggsave(here::here(fig_dir, "fig_s10.png"),
       p_go_full,
       bg = "white",
       height = 35,
       width = 35,
       units = "cm")

ggsave(here::here(fig_dir, "fig_s10.svg"),
       p_go_full,
       bg = "white",
       height = 35,
       width = 35,
       units = "cm")

# Figure S11 --------------------------------------------------------------------

# Representative DEGs and DARs

degs |> 
  filter(sig) |> 
  slice_max(pct.18, n = 10)

top_degs <- degs |> 
  filter(sig) |> 
  filter(pct.18 > 0.25 & pct.25 > 0.25) |> 
  group_by(cell_type) |> 
  slice_min(p_val_adj, n = 5)

top_degs <- top_degs[c(2,4,9,14,20,23,26,31,39), ]


DefaultAssay(dat) <- "RNA"

vln_data <- FetchData(object = dat, 
          vars = top_degs$gene,
          slot = "data") |>
  rownames_to_column("cell") |> 
  pivot_longer(-cell, names_to = "feature", values_to = "data") |> 
  full_join(dat@meta.data |> 
              rownames_to_column("cell") |> 
              select(cell, acc_temp, cell_type)) 

vln_plots <- list()

for (i in 1:nrow(top_degs)) {
  p_vln <- vln_data |> 
    filter(feature == top_degs$gene[i]) |> 
    filter(cell_type == top_degs$cell_type[i]) |> 
    ggplot(aes(
      x = acc_temp,
      y = data,
      fill = acc_temp)) +
    geom_violin(position = position_dodge(width = 0.9)) +
    geom_point(alpha = 0.6, shape = 21,
               position = position_jitterdodge(seed = 1, 
                                               dodge.width = 0.9, 
                                               jitter.width = 0.25)) +
    scale_fill_manual(values = acc_colors) +
    theme_minimal(base_size = 16) +
    labs(title = top_degs$cell_type[i],
         subtitle = top_degs$gene[i],
         y = "log-norm. exp.") +
    theme(axis.title.x = element_blank(),
          legend.position = "none")
  
  vln_plots[[i]] <- p_vln
}

fs11 <- cowplot::plot_grid(
  vln_plots[[1]],
  vln_plots[[2]],
  vln_plots[[3]],
  vln_plots[[4]],
  vln_plots[[5]],
  vln_plots[[6]],
  vln_plots[[7]],
  vln_plots[[8]],
  vln_plots[[9]],
  labels = LETTERS[1:9])


ggsave(here::here(fig_dir, "fig_s11.png"),
       fs11,
       height = 30,
       width = 40,
       units = "cm"
)
ggsave(here::here(fig_dir, "fig_s11.pdf"),
       fs11,
       height = 30,
       width = 40,
       units = "cm"
)
