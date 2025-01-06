# ------------------------------------------------------------------------------
# Plot differentially accessible regions
# TS O'Leary
# ------------------------------------------------------------------------------

# Load libraries
require(tidyverse)
require(Seurat)
require(Signac)
require(cowplot)

# Load plot themes 
source(here::here("src/04_plots/00_plot_themes.R"))

# Output fig dir
fig_dir <- "output/figs/dars"

# Load data
dat <- readRDS(here::here("data/processed/seurat_object/10_dat_linked.rds"))

degs <- readRDS(here::here("output/degs/degs_cell-type.rds")) |> 
  filter(p_val_adj < 0.05)
dars <- readRDS(here::here("output/dars/dars_cell-type.rds")) |> 
  filter(p_val_adj < 0.05) |> 
  separate(region, 
           sep = "-", 
           into = c("chrom", "start", "end"),
           convert = TRUE,
           remove = FALSE) |> 
  arrange(chrom, start)

# Label regions with overlapping genes
dars_df <- dars |> 
  select(region, chrom, start, end)
genes_df <- as_tibble(dat@assays$peaks@annotation) |> 
  dplyr::rename("chrom" = "seqnames") 

# Closest gene
gene_close <- valr::bed_closest(dars_df, genes_df) |> 
  arrange(region.x, abs(.dist)) |> 
  distinct(region.x, .keep_all = TRUE) |> 
  select(region.x, gene_name.y, .dist) |> 
  rename("region"= "region.x",
         "gene" = "gene_name.y",
         "dist" = ".dist") |> 
  filter(abs(dist) < 5000)

dars <- dars |> 
  left_join(gene_close, by = "region") |> 
  dplyr::select(cell_type, region, gene, dist, everything())
    
saveRDS(dars, here::here("output/dars/dars_cell-type_gene.rds"))

dars_degs <- dars |> 
  left_join(degs, by = c("gene", "cell_type"),
            suffix = c(".ATAC", ".RNA")) |> 
  filter(!is.na(padj.RNA))


saveRDS(dars_degs, here::here("output/dars/dars_degs.rds"))


# Volcanoplot
dars_pseudo <- readRDS(here::here("output/dars/dars.rds")) |> 
  rownames_to_column("region")

dp_p_df <- dars_pseudo |> 
  slice_min(p_val_adj, n = 7) |> 
  separate(region, 
           sep = "-", 
           into = c("chrom", "start", "end"),
           convert = TRUE,
           remove = FALSE) |> 
  arrange(chrom, start)

dp_fc_df <- dars_pseudo |> 
  slice_max(abs(avg_log2FC), n = 7) |> 
  separate(region, 
           sep = "-", 
           into = c("chrom", "start", "end"),
           convert = TRUE,
           remove = FALSE) |> 
  arrange(chrom, start)

dp_df <- bind_rows(dp_p_df, dp_fc_df)

dars_label <- valr::bed_closest(dp_df, genes_df) |> 
  arrange(region.x, abs(.dist)) |> 
  distinct(region.x, .keep_all = TRUE) |> 
  select(region.x, gene_name.y, .dist) |> 
  rename("region"= "region.x",
         "gene" = "gene_name.y",
         "dist" = ".dist") |> 
  filter(abs(dist) < 5000)


dars_pseudo |> 
  left_join(dars_label) |> 
  ggplot(aes(label = gene,
             y = -log10(p_val_adj),
             x = avg_log2FC,
             fill = p_val_adj < 0.05,
             size = p_val_adj < 0.05)) +
  geom_vline(xintercept = 0, color = "grey40") +
  geom_hline(yintercept = 0, color = "grey40") +
  geom_hline(yintercept = -log10(0.05),
             color = "grey50",
             linetype = 2) +
  geom_point(color = "grey90",
             stroke = 0.2,
             shape = 21,
             alpha = 0.9) +
  annotate(geom = "segment",
           x = 0.25, 
           xend = .75, 
           y = 35, 
           yend = 35,
           linewidth = 1,
           lineend = "round",
           linejoin = "round",
           color = "#43aa8b",
           alpha = 0.5,
           arrow = arrow(length = unit(0.2, "cm"))) +
  annotate(geom = "segment",
           x = -0.25, 
           xend = -0.75, 
           y = 35, 
           yend = 35,
           linewidth = 1,
           lineend = "round",
           linejoin = "round",
           color = "#f3722c",
           alpha = 0.5,
           arrow = arrow(length = unit(0.2, "cm"))) +
  annotate(geom = "text", 
           x = .5, 
           y = 37, 
           alpha = 0.5,
           label = "Higher in 18°C", 
           color =  "#43aa8b") +
  annotate(geom = "text", 
           x = -.5, 
           y = 37, 
           alpha = 0.5,
           label = "Higher in 25°C", 
           color = "#f3722c") +
  ggrepel::geom_label_repel(
    #data = dars_label,
                            color = "grey20",
                            fill = "grey90",
                            alpha = 0.8,
                            fontface = "bold",
                            size = 3,
                            max.iter = 10000,
                            min.segment.length = 0.1,
                            max.overlaps = 100) +
  scale_size_manual(values = c(1, 2)) +
  scale_fill_manual(values = c("grey80", "firebrick")) +
  scale_y_continuous(expand = c(0, 0.5),
                     limits = c(0, 40)) +
  scale_x_continuous(name = "18°C vs. 25°C\nlog2(fold change)") +
  cowplot::theme_minimal_grid() +
  theme(legend.position = "none") 

# Save volcano
ggsave(here::here(fig_dir, "pseudobulk_volcano.pdf"),
       height = 20,
       width = 20,
       units = "cm")
ggsave(here::here(fig_dir, "pseudobulk_volcano.png"),
       height = 20,
       width = 20,
       units = "cm")

  
# Plot number of DARs per cell-type --------------------------------------------
readRDS(here::here("output/dars/dars_cell-type.rds")) |> 
  filter(p_val_adj < 0.05) |> 
  group_by(cell_type) |> 
  tally() |> 
  full_join(color_cell_type) |> 
  filter(!is.na(n)) |> 
  mutate(colors = colorspace::desaturate(colors, amount = 0.6)) |> 
  mutate(cluster = fct_reorder(cell_type, n, .fun = "identity")) |> 
  ggplot() +
  geom_col(aes(x = n,
               y = cluster,
               fill = colors),
           na.rm = TRUE,
           alpha = 0.5,
           color = "grey20") +
  labs(y = "") +
  scale_x_continuous(position = "top",
                     name = "Number of differentially accessible regions",
                     expand = c(0, 0)) +
  scale_fill_identity() +
  cowplot::theme_minimal_vgrid()

# Save plot
ggsave(here::here(fig_dir, "dars_celltype.pdf"),
       height = 10,
       width = 20,
       units = "cm")
ggsave(here::here(fig_dir, "dars_celltype.png"),
       height = 10,
       width = 20,
       units = "cm")

# Plot percentage of of DARs per cell-type -------------------------------------
readRDS(here::here("output/dars/dars_cell-type.rds")) |> 
  group_by(cell_type) |> 
  add_tally(name = "n_peaks") |> 
  filter(p_val_adj < 0.05) |> 
  group_by(cell_type) |> 
  add_tally(name = "n_dars") |> 
  mutate(proportion_dars = n_dars/n_peaks) |> 
  select(cell_type, proportion_dars) |> 
  distinct(cell_type, .keep_all = TRUE) |> 
  ungroup() |> 
  mutate(cell_type = fct_reorder(cell_type, proportion_dars)) |>
  ggplot() +
  geom_col(aes(x = proportion_dars,
               y = cell_type),
           na.rm = TRUE,
           color = "grey20",
           fill = "grey80") +
  labs(y = "") +
  scale_x_continuous(position = "top",
                     name = "Percentage of differentially accessible regions",
                     labels = scales::percent_format(accuracy = 1),
                     expand = c(0, 0), 
                     limits = c(0, 0.09)) +
  cowplot::theme_minimal_vgrid()

# Save plot
ggsave(here::here(fig_dir, "dars_celltype_percent.pdf"),
       height = 10,
       width = 20,
       units = "cm")
ggsave(here::here(fig_dir, "dars_celltype_percent.png"),
       height = 10,
       width = 20,
       units = "cm")

dars_degs |> 
  ggplot() +
  geom_hline(yintercept = 0, color = "grey20") +
  geom_vline(xintercept = 0, color = "grey20") +
  geom_abline(slope = 1, 
              intercept = 0,
              linetype = 2,
              color = "grey50") +
  geom_point(aes(x = avg_log2FC_18_25.ATAC,
                 y = avg_log2FC_18_25.RNA)) +
  lims(x = c(-1.2, 1.2),
       y = c(-1.2, 1.2)) +
  labs(x = "ATAC peak log2FC",
       y = "RNA expression log2FC") +
  cowplot::theme_minimal_grid()

# Save plot
ggsave(here::here(fig_dir, "dars_degs_scatter.pdf"),
       height = 20,
       width = 20,
       units = "cm")
ggsave(here::here(fig_dir, "dars_degs_scatter.png"),
       height = 20,
       width = 20,
       units = "cm")


# all_dars <- readRDS(here::here("output/dars/dars_cell-type.rds")) |> 
#   separate(region, 
#            sep = "-", 
#            into = c("chrom", "start", "end"),
#            convert = TRUE,
#            remove = FALSE) |> 
#   arrange(chrom, start)|> 
#   select(region, chrom, start, end)
# all_dars <- valr::bed_closest(all_dars, genes_df) |>
#   arrange(region.x, abs(.dist)) |>
#   distinct(region.x, .keep_all = TRUE) |>
#   select(region.x, gene_name.y, .dist) |>
#   rename("region"= "region.x",
#          "gene" = "gene_name.y",
#          "dist" = ".dist") |>
#   filter(abs(dist) < 5000)

left_join(
  dars,
  readRDS(here::here("output/degs/degs_cell-type.rds")),
  by = c("cell_type", "gene")
)

dars |> 
  left_join(readRDS(here::here("output/degs/degs_cell-type.rds")), 
            by = c("gene", "cell_type"),
            suffix = c(".ATAC", ".RNA")) |> 
  filter(!is.na(padj.RNA))


# Plot Coverage Plot -----------------------------------------------------------
dat@meta.data <- dat@meta.data |> 
  mutate(celltype_acc = paste(cell_type, acc_temp, sep = "_"))
Idents(dat) <- "celltype_acc"

for (i in 1:nrow(dars)){
  # Coverage plot comparing 18°C and 25°C for a specific cell type
  cov_plot <- CoveragePlot(
    object = dat |> 
      subset(cell_type == dars$cell_type[i]),
    region = dars$region[i],
    annotation = FALSE,
    peaks = FALSE,
    extend.upstream = 5000,
    extend.downstream = 5000) &
    labs(title = str_to_sentence(dars$cell_type[i])) &
    scale_fill_manual(values = acc_colors) & 
    scale_y_continuous(expand = c(0, 0.05),
                       name = "Normalized signal") &
    theme_minimal_hgrid() &
    theme(legend.position = "none",
          strip.background = element_blank(),
          strip.text = element_blank())
  
  peak_rect <- annotate("rect",
    xmin = as.numeric(str_remove_all(
      str_extract(dars$region[i], "-[:digit:]+-"), "-")),
    xmax = as.numeric(str_remove_all(
      str_extract(dars$region[i], "-[:digit:]+$"), "-")),
    ymin = -Inf,
    ymax = Inf,
    fill = "grey75",
    alpha = 0.5)
  
  cov_plot$layers <- c(peak_rect, cov_plot$layers)
  
  gene_plot <- AnnotationPlot(
    object = dat |> 
      subset(cell_type == dars$cell_type[i]),
    region = dars$region[i],
    extend.upstream = 5000,
    extend.downstream = 5000) & 
    theme_cowplot() &
    theme(axis.line.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.y = element_blank(),
          axis.title.y = element_blank())
  
  peak_plot <- PeakPlot(
    object = dat |> 
      subset(cell_type == dars$cell_type[i]),
    region = dars$region[i],
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
    heights = c(8,0.1,2)
  )
  
  # Save plot
  ggsave(here::here(fig_dir, "coverage_plots", paste0(dars$region[i], ".png")),
         height = 20,
         width = 20,
         units = "cm")
  ggsave(here::here(fig_dir, "coverage_plots", paste0(dars$region[i], ".pdf")),
         height = 20,
         width = 20,
         units = "cm")
}


# Plot Coverage & Expression ---------------------------------------------------
dat@meta.data <- dat@meta.data |> 
  mutate(celltype_acc = paste(cell_type, acc_temp, sep = "_"))
Idents(dat) <- "celltype_acc"

for (i in 1:nrow(dars_degs)){
  # Coverage plot comparing 18°C and 25°C for a specific cell type
  cov_plot <- CoveragePlot(
    object = dat |> 
      subset(cell_type == dars_degs$cell_type[i]),
    region = dars_degs$region[i],
    annotation = FALSE,
    peaks = FALSE,
    extend.upstream = 5000,
    extend.downstream = 5000) &
    labs(title = str_to_sentence(dars_degs$cell_type[i])) &
    scale_fill_manual(values = acc_colors) & 
    scale_y_continuous(expand = c(0, 0.05),
                       name = "Normalized signal") &
    theme_minimal_hgrid() &
    theme(legend.position = "none",
          strip.background = element_blank(),
          strip.text = element_blank())
  
  peak_rect <- annotate("rect",
                        xmin = as.numeric(str_remove_all(
                          str_extract(dars_degs$region[i], "-[:digit:]+-"), "-")),
                        xmax = as.numeric(str_remove_all(
                          str_extract(dars_degs$region[i], "-[:digit:]+$"), "-")),
                        ymin = -Inf,
                        ymax = Inf,
                        fill = "grey75",
                        alpha = 0.5)
  
  cov_plot$layers <- c(peak_rect, cov_plot$layers)
  
  gene_plot <- AnnotationPlot(
    object = dat |> 
      subset(cell_type == dars_degs$cell_type[i]),
    region = dars_degs$region[i],
    extend.upstream = 5000,
    extend.downstream = 5000) & 
    theme_cowplot() &
    theme(axis.line.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.y = element_blank(),
          axis.title.y = element_blank())
  
  peak_plot <- PeakPlot(
    object = dat |> 
      subset(cell_type == dars_degs$cell_type[i]),
    region = dars_degs$region[i],
    assay = "peaks",
    extend.upstream = 5000,
    extend.downstream = 5000) &
    theme_cowplot() &
    theme(axis.line.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.y = element_blank(),
          axis.title.y = element_blank(),
          
          legend.position = "none")
  
  exp_plot <- ExpressionPlot(
    object = dat |> 
      subset(cell_type == dars_degs$cell_type[i]),
    features = dars_degs$gene[i],
    assay = "SCT") & 
    scale_fill_manual(values = acc_colors) &
    theme_cowplot() &
    theme(legend.position = "none",
          axis.text.x = element_text(size = 10),
          axis.line.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.y = element_blank(),
          strip.background = element_blank(),
          strip.text = element_blank(),
          axis.title.y = element_blank())
  
  CombineTracks(
    plotlist = list(
      cov_plot,
      peak_plot,
      gene_plot),
    expression.plot = exp_plot,
    heights = c(8, 0.1, 2),
    widths = c(10, 2)
  )
  
  # Save plot
  ggsave(here::here(fig_dir, "cov_exp", paste0(dars_degs$gene[i], ".png")),
         height = 20,
         width = 20,
         units = "cm")
  ggsave(here::here(fig_dir, "cov_exp", paste0(dars_degs$gene[i], ".pdf")),
         height = 20,
         width = 20,
         units = "cm")
}


