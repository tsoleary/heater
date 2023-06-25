# ------------------------------------------------------------------------------
# Find differentially expressed genes between conditions
# June 23, 2023
# TS O'Leary
# ------------------------------------------------------------------------------

# Load libraries
require(tidyverse)

# Load data
dat <- readRDS(here::here("data/processed/seurat_object/10_dat_linked.rds"))

degs_cluster <- readRDS(here::here("output/degs/degs_clusters.rds")) |> 
  filter(p_val_adj < 0.05)


write_tsv(
  degs_cluster |> 
    filter(cluster == "mesoderm primordium") |> 
    select(gene, avg_log2FC_18_25), 
  here::here("output/degs/degs_meso_lfc.tsv")
)

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

ggsave(here::here("output/figs/degs/scatter.pdf"),
       plot = p,
       width = 20,
       height = 20,
       units = "cm")
ggsave(here::here("output/figs/degs/scatter.png"),
       plot = p,
       width = 20,
       height = 20,
       units = "cm")

# Quick plotly with gene labels
plotly::ggplotly(p)
