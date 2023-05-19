# ------------------------------------------------------------------------------
# Knee plot figures
# March 27, 2023
# TS O'Leary
# ------------------------------------------------------------------------------

# Description -----
# Basic quality control metrics and filtering out non-cell barcodes.

# Load libraries
library(tidyverse)
library(Seurat)
library(Signac)

# Load  data
dat <- readRDS(here::here("data/processed/seurat_object/01_dat_raw.rds"))

# Load 10X filtered data
dat_10x <- readRDS(
  here::here("data/processed/seurat_object/02_dat_10x_cells.rds")
)

# Define count cutoffs
low_ATAC <- 800
low_RNA <- 350

# Knee-plot for ATAC -----------------------------------------------------------
# Create ranked data.frame
knee_atac_df <- tibble(
  nCount_ATAC = sort(dat$nCount_ATAC, 
                     decreasing = TRUE)
) %>%
  rownames_to_column("rank") %>%
  mutate(rank = as.numeric(rank)) %>%
  head(100000)

# Plot the knee_plot
p1 <- knee_atac_df %>%
  ggplot() +
  geom_point(aes(x = rank,
                 y = nCount_ATAC),
             color = "grey50",
             shape = 21) +
  geom_vline(aes(xintercept = 16000),
             linetype = 2,
             color = "orange",
             alpha = 0.8,
             linewidth = 1.1) +
  geom_hline(aes(yintercept = low_ATAC),
             linetype = 2,
             color = "red",
             alpha = 0.8,
             linewidth = 1.1) +
  scale_y_continuous(trans = "log10", 
                     breaks = c(100, 1000, 10000, 100000),
                     labels = c("100", "1k", "10k", "100k"),
                     name = "ATAC counts") +
  scale_x_continuous(trans = "log10",
                     breaks = c(1, 10, 100, 1000, 10000, 100000),
                     labels = c("1", "10", "100", "1k", "10k", "100k"),
                     name = "Barcodes rank ordered by total number of ATAC fragments") +
  cowplot::theme_minimal_grid()

# Knee-plot for RNA ------------------------------------------------------------
# Create ranked data.frame
knee_rna_df <- tibble(
  nCount_RNA = sort(dat$nCount_RNA, 
                    decreasing = TRUE)
) %>%
  rownames_to_column("rank") %>%
  mutate(rank = as.numeric(rank)) %>%
  head(100000)

# Plot the knee_plot
p2 <- knee_rna_df %>%
  ggplot() +
  geom_point(aes(x = rank,
                 y = nCount_RNA),
             color = "grey50",
             shape = 21) +
  geom_vline(aes(xintercept = 16000),
             linetype = 2,
             color = "orange",
             alpha = 0.8,
             linewidth = 1.1) +
  geom_hline(aes(yintercept = low_RNA),
             linetype = 2,
             color = "red",
             alpha = 0.8,
             linewidth = 1.1) +
  scale_y_continuous(trans = "log10", 
                     breaks = c(100, 1000, 10000, 100000),
                     labels = c("100", "1k", "10k", "100k"),
                     name = "RNA counts") +
  scale_x_continuous(trans = "log10",
                     breaks = c(1, 10, 100, 1000, 10000, 100000),
                     labels = c("1", "10", "100", "1k", "10k", "100k"),
                     name = "Barcodes rank ordered by total number of RNA reads") +
  cowplot::theme_minimal_grid()

# Knee plots together
cowplot::plot_grid(p1, p2,
                   nrow = 2)

# Save knee plot
ggsave(here::here("output/figs/qc/knee_plots.pdf"),
       height = 20,
       width = 20,
       units = "cm")
ggsave(here::here("output/figs/qc/knee_plots.png"),
       height = 20,
       width = 20,
       units = "cm")

# Scatter plot of ATAC & RNA counts --------------------------------------------

# Cells called by Cell Ranger ARC
dat@meta.data %>%
  filter(nCount_RNA >= 1 &
           nCount_ATAC >= 1) %>%
  rownames_to_column("cells") %>%
  mutate(cells_10x = ifelse(cells %in% Cells(dat_10x), 
                            "Cell Ranger ARC cells", 
                            "Non-cell barcodes")) %>%
  mutate(cells_10x = factor(cells_10x,
                            levels = c("Non-cell barcodes",
                                       "Cell Ranger ARC cells"))) %>%
  ggplot(aes(x = nCount_RNA,
             y = nCount_ATAC)) +
  geom_point(aes(color = cells_10x),
             shape = 21) +
  scale_x_continuous(trans = "log10", 
                     breaks = c(1, 10, 100, 1000, 10000, 100000),
                     labels = c("1", "10", "100", "1k", "10k", "100k"),
                     name = "RNA counts") +
  scale_y_continuous(trans = "log10", 
                     breaks = c(1, 10, 100, 1000, 10000, 100000),
                     labels = c("1", "10", "100", "1k", "10k", "100k"),
                     name = "ATAC counts") +
  scale_color_manual(values = c("grey90", "#91AB73"),
                     name = element_blank()) +
  cowplot::theme_minimal_grid() +
  theme(legend.position = "bottom")

# Save scatter plot
ggsave(here::here("output/figs/qc/scatter_rna_atac_10x_cells.pdf"),
       height = 20,
       width = 20,
       units = "cm")
ggsave(here::here("output/figs/qc/scatter_rna_atac_10x_cells.png"),
       height = 20,
       width = 20,
       units = "cm")

# Violin plot of QC metrics ----------------------------------------------------
p1 <- dat_10x@meta.data %>%
  mutate(project = "heater") %>%
  ggplot() +
  geom_violin(aes(y = nCount_RNA, 
                  x = project),
              fill = "grey80",
              color = "grey20") +
  scale_x_discrete(label = "RNA Counts",
                   name = element_blank()) +
  scale_y_continuous(labels = scales::label_number(suffix = "k", scale = 1e-3),
                     name = element_blank()) +
  theme_classic()

p2 <- dat_10x@meta.data %>%
  mutate(project = "heater") %>%
  ggplot() +
  geom_violin(aes(y = nCount_ATAC, 
                  x = project),
              fill = "grey80",
              color = "grey20") +
  scale_x_discrete(label = "ATAC Counts",
                   name = element_blank()) +
  scale_y_continuous(labels = scales::label_number(suffix = "k", scale = 1e-3),
                     name = element_blank()) +
  theme_classic()

p3 <- dat_10x@meta.data %>%
  mutate(project = "heater") %>%
  ggplot() +
  geom_violin(aes(y = TSS.enrichment, 
                  x = project),
              fill = "grey80",
              color = "grey20") +
  scale_x_discrete(label = "TSS Enrichment",
                   name = element_blank()) +
  ylab(element_blank()) +
  theme_classic()

p4 <- dat_10x@meta.data %>%
  mutate(project = "heater") %>%
  ggplot() +
  geom_violin(aes(y = nucleosome_signal, 
                  x = project),
              fill = "grey80",
              color = "grey20") +
  scale_x_discrete(label = "Nucleosome signal",
                   name = element_blank()) +
  ylab(element_blank()) +
  theme_classic()

p5 <- dat_10x@meta.data %>%
  mutate(project = "heater") %>%
  ggplot() +
  geom_violin(aes(y = percent.mt/100, 
                  x = project),
              fill = "grey80",
              color = "grey20") +
  scale_x_discrete(label = "Percent mitochondrial reads",
                   name = element_blank()) +
  scale_y_continuous(name = element_blank(),
                     labels = scales::percent_format(accuracy = 1)) +
  theme_classic()

# Combine all plots together
cowplot::plot_grid(p2, p3, p4, p1, p5,
                   nrow = 2)

# Save plot
ggsave(here::here("output/figs/qc/violin_qc_metrics.pdf"),
       height = 15,
       width = 30,
       units = "cm")
ggsave(here::here("output/figs/qc/violin_qc_metrics.png"),
       height = 15,
       width = 30,
       units = "cm")
