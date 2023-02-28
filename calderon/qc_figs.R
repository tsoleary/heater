# ------------------------------------------------------------------------------
# Quality control descriptive figures
# November 04, 2022
# TS O'Leary
# ------------------------------------------------------------------------------

# Load libraries
require(tidyverse)

# Load data

# Analyze data

#### SET UP ABOVE!1!!!!!! =---------------------asdfa0dsf0qw4-ri9afid

dat_comb@meta.data %>%
  ggplot(aes(x = nCount_RNA,
             y = nFeature_RNA)) +
  geom_point()

dat_comb@meta.data %>%
  ggplot(aes(x = nCount_RNA)) +
  geom_histogram(color = "grey20",
                 fill = "grey70",
                 bins = 100)

dat_comb@meta.data %>%
  ggplot(aes(x = nFeature_RNA)) +
  geom_histogram(color = "grey20",
                 fill = "grey70",
                 bins = 100)

# Plot Counts & Total features -------------------------------------------------
p1 <- dat_comb@meta.data %>%
  ggplot(aes(y = nCount_RNA,
             x = seurat_clusters)) +
  geom_violin(alpha = 0.8) +
  geom_boxplot(width = 0.1) +
  theme_bw() +
  theme(panel.border = element_blank(),
        axis.line = element_line(color = "grey20"))

p2 <- dat_comb@meta.data %>%
  ggplot(aes(y = nFeature_RNA,
             x = seurat_clusters)) +
  geom_violin(alpha = 0.8) +
  geom_boxplot(width = 0.1) +
  theme_bw() +
  theme(panel.border = element_blank(),
        axis.line = element_line(color = "grey20"))

cowplot::plot_grid(p1, 
                   p2,
                   nrow = 2)
