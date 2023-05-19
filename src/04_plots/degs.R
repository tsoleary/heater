# ------------------------------------------------------------------------------
# Find differentially expressed genes between conditions
# May 17, 2023
# TS O'Leary
# ------------------------------------------------------------------------------

# Load libraries
require(tidyverse)

# Load data
dat <- readRDS(here::here("data/processed/seurat_object/03_dat_clustered.rds"))

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
