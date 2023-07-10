# ------------------------------------------------------------------------------
# Find differentially expressed genes between conditions
# June 23, 2023
# TS O'Leary
# ------------------------------------------------------------------------------

# Load libraries
require(tidyverse)
require(Seurat)
require(Signac)

# Set output fig_dir
fig_dir <- "output/figs/degs"

# Load data
dat <- readRDS(here::here("data/processed/seurat_object/10_dat_linked.rds"))
degs <- readRDS(here::here("output/degs/degs_cell-type.rds")) |> 
  filter(p_val_adj < 0.05)

# Number of degs per cell type
degs_cluster |> 
  group_by(cell_type) |> 
  tally() |> 
  mutate(cell_type = fct_reorder(cell_type, n, .fun = "identity")) |> 
  ggplot() +
  geom_col(aes(x = n,
               y = cell_type),
           na.rm = TRUE,
           color = "grey20",
           fill = "grey80") +
  labs(y = "") +
  scale_x_continuous(position = "top",
                     name = "Number of differentially expressed genes",
                     expand = c(0, 0)) +
  cowplot::theme_cowplot()

ggsave(here::here(fig_dir, "degs_count.pdf"),
       height = 15,
       width = 30,
       units = "cm")
ggsave(here::here(fig_dir, "degs_count.png"),
       height = 15,
       width = 30,
       units = "cm")

# Set default assay
DefaultAssay(dat) <- "RNA"

# Reset primary idents to the acclimation temperatures
Idents(dat) <- "acc_temp"

# Average expression of RNA across all cells
avg_exp <- log10(AverageExpression(dat, verbose = FALSE)$RNA)
avg_exp <- avg_exp %>%
  as_tibble(rownames = "gene")

# Scatter plot of the average expression for all cells -------------------------
p <- avg_exp %>%
  ggplot() +
  geom_point(aes(x = `25°C`,
                 y = `18°C`,
                 label = gene)) +
  labs(x = "18°C",
       y = "25°C") +
  theme_minimal()

ggsave(here::here(fig_dir, "scatter.pdf"),
       plot = p,
       width = 20,
       height = 20,
       units = "cm")
ggsave(here::here(fig_dir, "scatter.png"),
       plot = p,
       width = 20,
       height = 20,
       units = "cm")

# Quick plotly with gene labels
plotly::ggplotly(p)

# degs_cluster |>
#   ggplot(aes(label = gene)) +
#   geom_hline(yintercept = -log10(0.05),
#              color = "grey50",
#              linetype = 2) +
#   geom_point(aes(y = -log10(padj),
#                  x = avg_log2FC_18_25,
#                  fill = padj < 0.05,
#                  size = padj < 0.05),
#              color = "grey90",
#              stroke = 0.2,
#              shape = 21) +
#   scale_size_manual(values = c(1, 2)) +
#   scale_fill_manual(values = c("grey80", "firebrick")) +
#   cowplot::theme_cowplot() +
#   theme(legend.position = "none") +
#   facet_wrap(~cell_type)
#
# ggsave(here::here(fig_dir, "celltype_volcano.pdf"),
#        height = 30,
#        width = 30,
#        units = "cm")
# ggsave(here::here(fig_dir, "celltype_volcano.png"),
#        height = 30,
#        width = 30,
#        units = "cm")



# Coverage plot example -----
Idents(dat) <- "cell_type"
gene <- "CG13427"

p1 <- CoveragePlot(
  object = dat |> subset(acc_temp == "18°C"),
  region = gene,
  features = gene,
  expression.assay = "SCT",
  extend.upstream = 500,
  extend.downstream = 1000
)

p2 <- CoveragePlot(
  object = dat |> subset(acc_temp == "25°C"),
  region = gene,
  features = gene,
  expression.assay = "SCT",
  extend.upstream = 500,
  extend.downstream = 1000
) + cowplot::theme_cowplot(font_size = 24)

# Plot the coverage
p <- cowplot::plot_grid(p1, p2,
                        labels = c("A. 18°C", "B. 25°C"))

# Create title for the plot
title <- cowplot::ggdraw() +
  cowplot::draw_label(
    paste0("Chromatin accessibility and expression of ", gene),
    fontface = 'bold',
    x = 0,
    hjust = 0
  ) +
  theme(
    plot.margin = margin(0, 0, 0, 7)
  )

# Put plot together with the title
cowplot::plot_grid(
  title, p,
  ncol = 1,
  rel_heights = c(0.1, 1)
)

ggsave(here::here(fig_dir, gene, "coverage.png"),
       height = 25,
       width = 40,
       units = "cm")

# dat |> 
#   subset(cell_type == "mesoderm prim.") |> 
#   CoveragePlot(
#     region = gene,
#     features = gene,
#     split.by = "acc_temp",
#     expression.assay = "SCT",
#     extend.upstream = 500,
#     extend.downstream = 1000
#   )



# Violin Plot
VlnPlot(dat,
        features = gene,
        split.by = "acc_temp",
        pt.size = 0,
        assay = "RNA") +
  scale_y_continuous(expand = c(0, 0.05),
                     position = "right") +
  scale_x_discrete(name = element_blank()) +
  scale_fill_manual(values = c("#43aa8b", "#f3722c")) +
  coord_flip() +
  cowplot::theme_cowplot()


ggsave(here::here(fig_dir, gene, "vln.png"),
       height = 25,
       width = 20,
       units = "cm")


