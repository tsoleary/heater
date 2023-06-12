# ------------------------------------------------------------------------------
# Plots for Quality Control metrics after cell calling
# May 25, 2023
# TS O'Leary
# ------------------------------------------------------------------------------

# Load libraries
require(tidyverse)

# Load libraries
library(tidyverse)
library(Seurat)
library(Signac)

# Load raw data
dat_raw <- readRDS(
  here::here("data/processed/seurat_object/00_dat_raw.rds")
)

# Load filtered cells - mutiplets removed
dat <- readRDS(
  here::here("data/processed/seurat_object/06_dat_qc.rds")
)

# Define count cutoffs
low_ATAC <- 800
low_RNA <- 200

# Cells called by Cell Ranger ARC & low count threshold
dat_raw@meta.data |>
  filter(nCount_RNA >= 1 &
           nCount_ATAC >= 1) |>
  rownames_to_column("cells") |>
  mutate(cells_10x = ifelse(cells %in% Cells(dat) &
                              nCount_RNA >= low_RNA &
                              nCount_ATAC >= low_ATAC, 
                            "Cells", 
                            "Non-cell barcodes")) |>
  mutate(cells_10x = factor(cells_10x,
                            levels = c("Non-cell barcodes",
                                       "Cells"))) |>
  ggplot(aes(x = nCount_RNA,
             y = nCount_ATAC)) +
  geom_point(aes(color = cells_10x),
             shape = 21) +
  geom_segment(aes(x = low_RNA, xend = low_RNA,
                   y = low_ATAC, yend = 150000),
               linetype = 2,
               color = "grey50",
               alpha = 0.8,
               linewidth = 1.1) +
  geom_segment(aes(y = low_ATAC, yend = low_ATAC,
                   x = low_RNA, xend = 50000),
               linetype = 2,
               color = "grey50",
               alpha = 0.8,
               linewidth = 1.1) +
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
ggsave(here::here("output/figs/qc/scatter_rna_atac_final_filter.pdf"),
       height = 20,
       width = 20,
       units = "cm")
ggsave(here::here("output/figs/qc/scatter_rna_atac_final_filter.png"),
       height = 20,
       width = 20,
       units = "cm")

# Violin plot after filtering
p1 <- dat@meta.data |>
  mutate(project = "heater") |>
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

p2 <- dat@meta.data |>
  mutate(project = "heater") |>
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

p3 <- dat@meta.data |>
  mutate(project = "heater") |>
  ggplot() +
  geom_violin(aes(y = TSS.enrichment, 
                  x = project),
              fill = "grey80",
              color = "grey20") +
  scale_x_discrete(label = "TSS Enrichment",
                   name = element_blank()) +
  ylab(element_blank()) +
  theme_classic()

p4 <- dat@meta.data |>
  mutate(project = "heater") |>
  ggplot() +
  geom_violin(aes(y = nucleosome_signal, 
                  x = project),
              fill = "grey80",
              color = "grey20") +
  scale_x_discrete(label = "Nucleosome signal",
                   name = element_blank()) +
  ylab(element_blank()) +
  theme_classic()

p5 <- dat@meta.data |>
  mutate(project = "heater") |>
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
ggsave(here::here("output/figs/qc/violin_qc_metrics_post_filter.pdf"),
       height = 15,
       width = 30,
       units = "cm")
ggsave(here::here("output/figs/qc/violin_qc_metrics_post_filter.png"),
       height = 15,
       width = 30,
       units = "cm")

# TSS Enrichment plot for each sample ------------------------------------------
TSSPlot(dat, 
        group.by = "sample_name") +
  scale_color_manual(values = rep("grey50",4)) +
  cowplot::theme_cowplot() +
  theme(legend.position = "none")


# Save plot
ggsave(here::here("output/figs/qc/tss_enrichment_samples.pdf"),
       height = 15,
       width = 25,
       units = "cm")
ggsave(here::here("output/figs/qc/tss_enrichment_samples.png"),
       height = 15,
       width = 25,
       units = "cm")

# Fragment histogram for each sample -------------------------------------------
FragmentHistogram(dat, 
                  group.by = "sample_name")

# Save plot
ggsave(here::here("output/figs/qc/fragment_hist.pdf"),
       height = 15,
       width = 25,
       units = "cm")
ggsave(here::here("output/figs/qc/fragment_hist.png"),
       height = 15,
       width = 25,
       units = "cm")

# FRiP
dat@meta.data |>
  group_by(sample_name) |>
  summarise(FRiP = median(FRiP))  |>
  mutate(sample_name = factor(sample_name, levels = c(c("25C_Rep2",
                                                        "25C_Rep1",
                                                        "18C_Rep2",
                                                        "18C_Rep1")))) |>
  ggplot(aes(x = FRiP,
             y = sample_name,
             label = round(FRiP, digits = 3))) +
  geom_col(color = "grey20",
           fill = "grey80",
           width = 0.8)+
  scale_x_continuous(expand = c(0, 0),
                     name = "Median per cell fraction of reads in peaks (FRiP)") +
  geom_label(nudge_x = -0.035) +
  scale_y_discrete(name = element_blank(),
                   labels = c("25°C Rep 2",
                              "25°C Rep 1",
                              "18°C Rep 2",
                              "18°C Rep 1")) +
  cowplot::theme_cowplot() 
 
# Save plot
ggsave(here::here("output/figs/qc/frip_samples.pdf"),
       height = 10,
       width = 20,
       units = "cm")
ggsave(here::here("output/figs/qc/frip_samples.png"),
       height = 10,
       width = 20,
       units = "cm")

# FRiT
dat@meta.data |>
  group_by(sample_name) |>
  summarise(FRiT = median(FRiT))  |>
  mutate(sample_name = factor(sample_name, levels = c(c("25C_Rep2",
                                                        "25C_Rep1",
                                                        "18C_Rep2",
                                                        "18C_Rep1")))) |>
  ggplot(aes(x = FRiT,
             y = sample_name,
             label = round(FRiT, digits = 3))) +
  geom_col(color = "grey20",
           fill = "grey80",
           width = 0.8)+
  scale_x_continuous(expand = c(0, 0),
                     name = "Median per cell fraction of reads in TSS (FRiT)") +
  geom_label(nudge_x = -0.05) +
  scale_y_discrete(name = element_blank(),
                   labels = c("25°C Rep 2",
                              "25°C Rep 1",
                              "18°C Rep 2",
                              "18°C Rep 1")) +
  cowplot::theme_cowplot() 

# Save plot
ggsave(here::here("output/figs/qc/frit_samples.pdf"),
       height = 10,
       width = 20,
       units = "cm")
ggsave(here::here("output/figs/qc/frit_samples.png"),
       height = 10,
       width = 20,
       units = "cm")
  
# Cells per sample
dat@meta.data |>
  group_by(sample_name) |>
  mutate(sample_name = factor(sample_name, levels = c(c("25C_Rep2",
                                                        "25C_Rep1",
                                                        "18C_Rep2",
                                                        "18C_Rep1")))) |>
  tally() |>
  ggplot(aes(x = n,
             y = sample_name,
             label = n)) +
  geom_col(color = "grey20",
           fill = "grey80",
           width = 0.8)+
  scale_x_continuous(expand = c(0, 0),
                     name = "Number of cells") +
  geom_label(nudge_x = -200) +
  scale_y_discrete(name = element_blank(),
                   labels = c("25°C Rep 2",
                              "25°C Rep 1",
                              "18°C Rep 2",
                              "18°C Rep 1")) +
  cowplot::theme_cowplot() 

# Save plot
ggsave(here::here("output/figs/qc/cells_sample.pdf"),
       height = 10,
       width = 20,
       units = "cm")
ggsave(here::here("output/figs/qc/cells_sample.png"),
       height = 10,
       width = 20,
       units = "cm")

# Median RNA and ATAC counts
dat@meta.data |>
  group_by(sample_name) |>
  summarise(nCount_RNA_median = median(nCount_RNA),
            nCount_ATAC_median = median(nCount_ATAC))  |>
  mutate(sample_name = factor(sample_name, levels = c(c("25C_Rep2",
                                                        "25C_Rep1",
                                                        "18C_Rep2",
                                                        "18C_Rep1")))) |>
  ggplot(aes(x = nCount_RNA_median,
             y = sample_name,
             label = nCount_RNA_median)) +
  geom_col(color = "grey20",
           fill = "grey80",
           width = 0.8)+
  scale_x_continuous(expand = c(0, 0),
                     name = "Median per cell RNA count") +
  geom_label(nudge_x = -75) +
  scale_y_discrete(name = element_blank(),
                   labels = c("25°C Rep 2",
                              "25°C Rep 1",
                               "18°C Rep 2",
                              "18°C Rep 1")) +
  cowplot::theme_cowplot() 

# Save plot
ggsave(here::here("output/figs/qc/median_rna_sample.pdf"),
       height = 10,
       width = 20,
       units = "cm")
ggsave(here::here("output/figs/qc/median_rna_sample.png"),
       height = 10,
       width = 20,
       units = "cm")

# Median RNA and ATAC counts
dat@meta.data |>
  group_by(sample_name) |>
  summarise(nCount_RNA_median = median(nCount_RNA),
            nCount_ATAC_median = median(nCount_ATAC))  |>
  mutate(sample_name = factor(sample_name, levels = c(c("25C_Rep2",
                                                        "25C_Rep1",
                                                        "18C_Rep2",
                                                        "18C_Rep1")))) |>
  ggplot(aes(x = nCount_ATAC_median,
             y = sample_name,
             label = nCount_ATAC_median)) +
  geom_col(color = "grey20",
           fill = "grey80",
           width = 0.8)+
  scale_x_continuous(expand = c(0, 0),
                     name = "Median per cell ATAC count") +
  geom_label(nudge_x = -300) +
  scale_y_discrete(name = element_blank(),
                   labels = c("25°C Rep 2",
                              "25°C Rep 1",
                              "18°C Rep 2",
                              "18°C Rep 1")) +
  cowplot::theme_cowplot() 

# Save plot
ggsave(here::here("output/figs/qc/median_atac_sample.pdf"),
       height = 10,
       width = 20,
       units = "cm")
ggsave(here::here("output/figs/qc/median_atac_sample.png"),
       height = 10,
       width = 20,
       units = "cm")


