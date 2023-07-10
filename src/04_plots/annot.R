# ------------------------------------------------------------------------------
# Cluster and expression plots
# June 23, 2023
# TS O'Leary
# ------------------------------------------------------------------------------

# Load libraries
library(tidyverse)
library(Seurat)
library(Signac)

# Load data
dat <- readRDS(
  here::here("data/processed/seurat_object/09_dat_annot.rds")
)

# Set output fig_dir
fig_dir <- "output/figs/annot"

# Set cell-type as the primary identity
Idents(dat) <- "cell_type"

# 
plot <- DimPlot(dat, 
        split.by = "acc_temp",
        repel = TRUE,
        label.color = "grey20",
        pt.size = 1.5) +
  theme_void() +
  theme(legend.position = "none",
        strip.text = element_text(size = 16, face = "bold"),
        strip.background = element_rect(fill = "grey95", color = "grey95")) +
  scale_color_manual(
    values = c("#ADD9F4",
               "#57A4B2",
               "#D39C9E",
               "#FEF29A",
               "#F9DCEE",
               "grey90",
               "#819FC5",
               "#A7BF9B",
               "#F9FADC",
               "grey95",
               "pink")
    )

LabelClusters(plot = plot, 
              id = "ident",
              box = TRUE,
              alpha = 0.75,
              color = c("grey10"))

# Save plot
ggsave(here::here(fig_dir, "umap_18_25_split.pdf"),
       height = 15,
       width = 30,
       units = "cm")
ggsave(here::here(fig_dir, "umap_18_25_split.png"),
       height = 15,
       width = 30,
       units = "cm")

# Quick bar plot counting the number of cells for cell type
dat@meta.data |>
  mutate(cell_type = factor(cell_type,
                            levels = c("unknown",
                                       "amnioserosa",
                                       "yolk nuclei",
                                       "tracheal prim.",
                                       "ventral nerve cord prim.",
                                       "endoderm prim.",
                                       "foregut/hindgut prim.",
                                       "peripheral nervous system prim.",
                                       "ectoderm prim.",
                                       "mesoderm prim.",
                                       "ubiquitous"))) |> 
  ggplot(aes(y = cell_type,
             fill = acc_temp)) +
  geom_bar(position = "dodge",
           color = "grey50") +
  labs(y = "",
       x = "Number of cells") +
  scale_fill_manual(name = element_blank(),
                    values = c("#43aa8b", "#f3722c")) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.05)),
                     position = "top") +
  cowplot::theme_minimal_vgrid()

# Save plot
ggsave(here::here(fig_dir, "cells_per_celltype.pdf"),
       height = 15,
       width = 30,
       units = "cm")
ggsave(here::here(fig_dir, "cells_per_celltype.png"),
       height = 15,
       width = 30,
       units = "cm")


# Quick bar plot counting the number of cells for cell type
dat@meta.data |>
  mutate(cell_type = factor(cell_type,
                            levels = c("unknown",
                                       "amnioserosa",
                                       "yolk nuclei",
                                       "tracheal prim.",
                                       "ventral nerve cord prim.",
                                       "endoderm prim.",
                                       "foregut/hindgut prim.",
                                       "peripheral nervous system prim.",
                                       "ectoderm prim.",
                                       "mesoderm prim.",
                                       "ubiquitous"))) |> 
  ggplot(aes(y = cell_type)) +
  geom_bar(fill = "grey80",
           color = "grey20") +
  labs(y = "",
       x = "Number of cells") +
  scale_fill_manual(name = element_blank(),
                    values = c("#43aa8b", "#f3722c")) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.05)),
                     position = "top") +
  cowplot::theme_minimal_vgrid()

# Save plot
ggsave(here::here(fig_dir, "cells_per_celltype_combined.pdf"),
       height = 15,
       width = 30,
       units = "cm")
ggsave(here::here(fig_dir, "cells_per_celltype_combined.png"),
       height = 15,
       width = 30,
       units = "cm")


# Quick bar plot counting the number of cells for cell type

################################################################################
# Hard code fix later-----
################################################################################
dat@meta.data |>
  mutate(cell_type = factor(cell_type,
                            levels = c("unknown",
                                       "amnioserosa",
                                       "yolk nuclei",
                                       "tracheal prim.",
                                       "ventral nerve cord prim.",
                                       "endoderm prim.",
                                       "foregut/hindgut prim.",
                                       "peripheral nervous system prim.",
                                       "ectoderm prim.",
                                       "mesoderm prim.",
                                       "ubiquitous"))) |> 
  group_by(cell_type,
           acc_temp) |> 
  tally() |> 
  mutate(prop = ifelse(acc_temp == "18°C", n/6284,
                       n/5756)) |> 
  ggplot(aes(x = prop,
             y = cell_type)) +
  geom_col(aes(fill = acc_temp),
           color = "grey80",
           position = "dodge") +
  labs(y = "",
       x = "") +
  scale_fill_manual(name = element_blank(),
                    values = c("#43aa8b", "#f3722c")) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.05)),
                     labels = scales::percent,
                     position = "top") +
  cowplot::theme_minimal_vgrid()

# Save plot
ggsave(here::here(fig_dir, "cells_per_celltype_prop.pdf"),
       height = 15,
       width = 30,
       units = "cm")
ggsave(here::here(fig_dir, "cells_per_celltype_prop.png"),
       height = 15,
       width = 30,
       units = "cm")

# Number of genes detected per cell-type
dat@meta.data |> 
  ggplot(aes(y = cell_type,
             x = nFeature_SCT)) +
  geom_violin(color = "grey20",
              fill = "grey80") +
  geom_boxplot(color = "grey40",
               fill = "grey95",
               width = 0.1,
               outlier.shape = NA) +
  scale_x_continuous(name = "Number of expressed genes detected per cell",
                     limits = c(0, 1650),
                     expand = c(0, 0.05)) +
  scale_y_discrete(name = element_blank()) +
  theme_minimal_vgrid()

# Save plot
ggsave(here::here(fig_dir, "genes_per_celltype.png"),
       height = 15,
       width = 25,
       units = "cm")
ggsave(here::here(fig_dir, "genes_per_celltype.pdf"),
       height = 15,
       width = 25,
       units = "cm")

# Number of peaks detected per cell-type
dat@meta.data |> 
  ggplot(aes(y = cell_type,
             x = nFeature_ATAC)) +
  geom_violin(color = "grey20",
              fill = "grey80") +
  geom_boxplot(color = "grey40",
               fill = "grey95",
               width = 0.1,
               outlier.shape = NA) +
  scale_x_continuous(name = "Number of peaks detected per cell",
                     limits = c(0, 6200),
                     expand = c(0, 0.05),) +
  scale_y_discrete(name = element_blank()) +
  theme_minimal_vgrid()

# Save plot
ggsave(here::here(fig_dir, "peaks_per_celltype.png"),
       height = 15,
       width = 25,
       units = "cm")
ggsave(here::here(fig_dir, "peaks_per_celltype.pdf"),
       height = 15,
       width = 25,
       units = "cm")

