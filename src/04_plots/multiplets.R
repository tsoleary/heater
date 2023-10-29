# ------------------------------------------------------------------------------
# Plot marked multiplets onto the data
# TS O'Leary
# ------------------------------------------------------------------------------

# Load libraries
require(Seurat)
require(Signac)
require(tidyverse)

# Output fig dir
fig_dir <- "output/figs/qc/multiplets"

# Load data
dat <- readRDS(
  here::here("data/processed/seurat_object/05_dat_multiplets.rds")
)

# Visualize multiplets on the UMAP projection ----------------------------------
DimPlot(dat, 
        group.by = "multiplet_rm") +
  labs(title = element_blank()) +
  scale_color_manual(name = element_blank(),
                     values = c("firebrick", "grey80")) +
  theme_void() +
  theme(legend.position = "bottom")

# Save the plot
ggsave(here::here(fig_dir, "multiplet_umap.pdf"),
       width = 15,
       height = 15,
       units = "cm")
ggsave(here::here(fig_dir, "multiplet_umap.png"),
       width = 15,
       height = 15,
       units = "cm")

# Counts mapped on to UMAP projection

# ATAC counts
p1 <- FeaturePlot(dat, 
                  feature = "nCount_ATAC") +
  labs(title = "ATAC Counts") +
  scale_color_distiller(palette = 2,
                        direction = 1,
                        trans = "log",
                        values = c(0, 0.7, 1),
                        breaks = c(0, 1000, 10000, 100000),
                        labels = c("0", "1k", "10k", "100k")) +
  theme_void() +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        legend.key.width = unit(1.25, "cm"),
        plot.title = element_text(hjust = 0.5))

p1 <- FeaturePlot(dat, 
                  feature = "nCount_ATAC") +
  labs(title = "ATAC Counts") +
  scale_color_gradient(low = "grey90",
                       high = "darkgreen",
                       trans = "log",
                       breaks = c(0, 1000, 10000, 100000),
                       labels = c("0", "1k", "10k", "100k")) +
  theme_void() +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        legend.key.width = unit(1.25, "cm"),
        plot.title = element_text(hjust = 0.5))

# RNA counts
p2 <- FeaturePlot(dat, 
            feature = "nCount_RNA") +
  scale_color_gradient(low = "grey90",
                       high = "darkgreen",
                       trans = "log",
                       breaks = c(0, 1000, 10000, 100000),
                       labels = c("0", "1k", "10k", "100k")) +
  labs(title = "RNA Counts") +
  theme_void() +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        legend.key.width = unit(1.25, "cm"),
        plot.title = element_text(hjust = 0.5))

# Combine plots
cowplot::plot_grid(p1, p2, 
                   nrow = 1)

# Save plot
ggsave(here::here(fig_dir, "counts_umap.pdf"),
       width = 25,
       height = 12,
       units = "cm")
ggsave(here::here(fig_dir, "counts_umap.png"),
       width = 20,
       height = 12,
       units = "cm")

# Histogram out counts for RNA and ATAC libraries with multiplets --------------

# RNA count multiplet histogram
p1 <- dat@meta.data |>
  mutate(doublet = factor(multiplet_rm,
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
  mutate(doublet = factor(multiplet_rm,
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
ggsave(here::here(fig_dir, "multiplets_hist.pdf"),
       width = 20,
       height = 15,
       units = "cm")
ggsave(here::here(fig_dir, "multiplets_hist.png"),
       width = 20,
       height = 15,
       units = "cm")

# Histogram out counts for RNA and ATAC libraries without multiplets -----------

# RNA count multiplet histogram
p1 <- dat@meta.data |>
  mutate(doublet = factor(multiplet_rm,
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
  mutate(doublet = factor(multiplet_rm,
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
ggsave(here::here(fig_dir, "singlets_only_hist.pdf"),
       width = 20,
       height = 15,
       units = "cm")
ggsave(here::here(fig_dir, "singlets_only_hist.png"),
       width = 20,
       height = 15,
       units = "cm")
