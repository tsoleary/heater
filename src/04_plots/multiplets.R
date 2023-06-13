# ------------------------------------------------------------------------------
# Plot marked multiplets onto the data
# May 23, 2023
# TS O'Leary
# ------------------------------------------------------------------------------

# Load libraries
require(tidyverse)

# Load data
dat <- readRDS(
  here::here("data/processed/seurat_object/05_dat_multiplets.rds")
  ) 

dat@meta.data <- dat@meta.data |>
  mutate(doublet = ifelse(doublet_amulet == "Multiplet" | 
                            doublet_finder == "Multiplet",
                          "Multiplet",
                          "Singlet"))

# Visualize multiplets on the UMAP projection ----------------------------------
DimPlot(dat, 
        group.by = "doublet") +
  labs(title = element_blank()) +
  scale_color_manual(name = element_blank(),
                     values = c("firebrick", "grey80")) +
  theme_void() +
  theme(legend.position = "bottom")

# Save the plot
ggsave(here::here("output/figs/qc/multiplet_umap.pdf"),
       width = 15,
       height = 15,
       units = "cm")
ggsave(here::here("output/figs/qc/multiplet_umap.png"),
       width = 15,
       height = 15,
       units = "cm")

# Counts mapped on to UMAP projection

# ATAC counts
p1 <- FeaturePlot(dat, 
        feature = "nCount_ATAC") +
  labs(title = element_blank()) +
  scale_color_gradient(name = "ATAC Counts",
                       low = "grey90", 
                       high = "darkgreen") +
  theme_void()

# RNA counts
p2 <- FeaturePlot(dat, 
            feature = "nCount_RNA") +
  scale_color_gradient(name = "RNA Counts",
                       low = "grey90", 
                       high = "darkgreen") +
  labs(title = element_blank()) +
  theme_void()

# Combine plots
cowplot::plot_grid(p1, p2, 
                   nrow = 1)

# Save plot
ggsave(here::here("output/figs/qc/counts_umap.pdf"),
       width = 25,
       height = 15,
       units = "cm")
ggsave(here::here("output/figs/qc/counts_umap.png"),
       width = 20,
       height = 15,
       units = "cm")

# Histogram out counts for RNA and ATAC libraries with multiplets --------------

# RNA count multiplet histogram
p1 <- dat@meta.data |>
  mutate(doublet = factor(doublet,
                          levels = c("Singlet", "Multiplet"))) |>
  ggplot() +
  geom_histogram(aes(x = nCount_RNA,
                     fill = doublet),
                 color = "grey30",
                 position = "dodge",
                 bins = 40) +
  scale_fill_manual(values = c("grey50", "firebrick"),
                    name = element_blank()) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(y = "Count", x = "RNA counts per barcode") +
  cowplot::theme_cowplot()

# ATAC count multiplet histogram
p2 <- dat@meta.data |>
  mutate(doublet = factor(doublet,
                          levels = c("Singlet", "Multiplet"))) |>
  ggplot() +
  geom_histogram(aes(x = nCount_ATAC,
                     fill = doublet),
                 color = "grey30",
                 position = "dodge",
                 bins = 40) +
  scale_fill_manual(values = c("grey50", "firebrick"),
                    name = element_blank()) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(y = "Count", x = "ATAC counts per barcode") +
  cowplot::theme_cowplot()


# Extract the legend from one of the plots
legend <- cowplot::get_legend(
  # Create some space to the left of the legend
  p1 + theme(legend.box.margin = margin(0, 0, 0, 12))
)

# Combine two plots without their legends
plots <- cowplot::plot_grid(p1 + theme(legend.position = "none"), 
                            p2 + theme(legend.position = "none"),
                            nrow = 2,
                            labels = c("A", "B"))

# Add the legend to the plots we made earlier
cowplot::plot_grid(plots, legend, rel_widths = c(2, .5))

# Save the plot
ggsave(here::here("output/figs/qc/multiplets_hist.pdf"),
       width = 20,
       height = 15,
       units = "cm")
ggsave(here::here("output/figs/qc/multiplets_hist.png"),
       width = 20,
       height = 15,
       units = "cm")

# Histogram out counts for RNA and ATAC libraries without multiplets -----------

# RNA count multiplet histogram
p1 <- dat@meta.data |>
  mutate(doublet = factor(doublet,
                          levels = c("Singlet", "Multiplet"))) |>
  filter(doublet == "Singlet" & nCount_RNA < 10000) |>
  ggplot() +
  geom_histogram(aes(x = nCount_RNA),
                 fill = "grey50",
                 color = "grey30",
                 position = "dodge",
                 bins = 30) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(y = "Count", x = "RNA counts per barcode") +
  cowplot::theme_cowplot()

# ATAC count multiplet histogram
p2 <- dat@meta.data |>
  mutate(doublet = factor(doublet,
                          levels = c("Singlet", "Multiplet"))) |>
  filter(doublet == "Singlet" & nCount_RNA < 10000) |>
  ggplot() +
  geom_histogram(aes(x = nCount_ATAC),
                 fill = "grey50",
                 color = "grey30",
                 position = "dodge",
                 bins = 30) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(y = "Count", x = "ATAC counts per barcode") +
  cowplot::theme_cowplot() 

# Combine two plots without their legends
cowplot::plot_grid(p1, 
                   p2,
                   nrow = 2,
                   labels = c("A", "B"))

# Save the plot
ggsave(here::here("output/figs/qc/singlets_only_hist.pdf"),
       width = 20,
       height = 15,
       units = "cm")
ggsave(here::here("output/figs/qc/singlets_only_hist.png"),
       width = 20,
       height = 15,
       units = "cm")
