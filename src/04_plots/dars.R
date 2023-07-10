# ------------------------------------------------------------------------------
# Plot differentially accessible regions
# July 07, 2023
# TS O'Leary
# ------------------------------------------------------------------------------

# Load libraries
require(tidyverse)
require(Seurat)
require(Signac)
require(cowplot)

# Set output fig_dir
fig_dir <- "output/figs/dars"

# Load data
dat <- readRDS(here::here("data/processed/seurat_object/10_dat_linked.rds"))
degs <- readRDS(here::here("output/degs/degs_cell-type.rds")) |> 
  filter(p_val_adj < 0.05)
dars <- readRDS(here::here("output/dars/dars_cell-type.rds")) |> 
  filter(padj < 0.05) |> 
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
  rename("seqnames" = "chrom")

# Closest gene
gene_close <- valr::bed_closest(dars_df, genes_df) |> 
  arrange(region.x, abs(.dist)) |> 
  distinct(region.x, .keep_all = TRUE) |> 
  select(region.x, gene_name.y, .dist) |> 
  rename("region.x" = "region",
         "gene_name.y" = "gene",
         ".dist" = "dist") |> 
  filter(abs(dist) < 5000)

dars <- dars |> 
  left_join(gene_close, by = "region") |> 
  dplyr::select(cell_type, region, gene, dist, everything())

saveRDS(dars, here::here("output/dars/dars_cell-type_gene.rds"))

dars_degs <- dars |> 
  left_join(degs, by = c("gene", "cell_type"),
            suffix = c(".ATAC", ".RNA")) |> 
  filter(!is.na(padj.RNA))

  
# Plot number of DARs per cell-type --------------------------------------------
dars |> 
  group_by(cell_type) |> 
  tally() |> 
  mutate(cluster = fct_reorder(cell_type, n, .fun = "identity")) |> 
  ggplot() +
  geom_col(aes(x = n,
               y = cluster),
           na.rm = TRUE,
           color = "grey20",
           fill = "grey80") +
  labs(y = "") +
  scale_x_continuous(position = "top",
                     name = "Number of differentially accessible regions",
                     expand = c(0, 0)) +
  cowplot::theme_cowplot()

# Save plot
ggsave(here::here(fig_dir, "dars_celltype.pdf"),
       height = 10,
       width = 20,
       units = "cm")
ggsave(here::here(fig_dir, "dars_celltype.png"),
       height = 10,
       width = 20,
       units = "cm")


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
    scale_fill_manual(values = c("#43aa8b", "#f3722c")) & 
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
    scale_fill_manual(values = c("#43aa8b", "#f3722c")) & 
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
    scale_fill_manual(values = c("#43aa8b", "#f3722c")) &
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


