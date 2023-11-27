# ------------------------------------------------------------------------------
# Final figures
# TS O'Leary
# ------------------------------------------------------------------------------

# Load libraries
library(tidyverse)
library(Seurat)
library(Signac)
library(cowplot)

# Output fig dir
fig_dir <- "output/figs/final"

# Phenotype data
pheno <- read_csv(
  here::here("data/raw/pheno/acc_hs_survival.csv")
)

# Load raw data
dat_raw <- readRDS(
  here::here("data/processed/seurat_object/00_dat_raw.rds")
)

# Load final seurat object
dat <- readRDS(
  here::here("data/processed/seurat_object/10_dat_linked.rds")
)

# Recode sample name for prettier plotting
dat@meta.data <- dat@meta.data |> 
  mutate(sample_name_fig = case_when(sample_name == "18C_Rep1" ~ "18°C Rep 1",
                                     sample_name == "18C_Rep2" ~ "18°C Rep 2",
                                     sample_name == "25C_Rep1" ~ "25°C Rep 1",
                                     sample_name == "25C_Rep2" ~ "25°C Rep 2")) |> 
  mutate(sample_name_fig = factor(sample_name_fig, levels = c(c("25°C Rep 2",
                                                        "25°C Rep 1",
                                                        "18°C Rep 2",
                                                        "18°C Rep 1"))))

# ------------------------------------------------------------------------------
# Manuscript figures -----------------------------------------------------------
# ------------------------------------------------------------------------------

# Figure 1 ---------------------------------------------------------------------

p_exp <- ggdraw() + 
  draw_image(here::here("output/figs/exp_design/acc_design.png"))


# Plot Canton S hatching success after heat shock
p_pheno <- pheno |>
  mutate(acc_temp = paste0(acc_temp, "°C")) |>
  ggplot(aes(x = acc_temp,
             y = n_hatched/n_eggs*100)) +
  geom_boxplot(width = 0.3,
               color = "grey20",
               outlier.shape = NA) +
  ggbeeswarm::geom_beeswarm(size = 3,
                            cex = 2,
                            shape = 21,
                            fill = "grey50",
                            color = "grey20",
                            alpha = 0.9) +
  xlab("Acclimation temperature") +
  scale_y_continuous(name = "Egg hatching success after acute heat shock",
                     limits = c(0, 100),
                     expand = expansion(mult = c(0, 0.05)),
                     labels = function(x) paste0(x, "%")) + 
  cowplot::theme_minimal_hgrid()

# Knee plots together
cowplot::plot_grid(
  p_exp, 
  p_pheno,
  labels = c("A", "B"),
  rel_heights = c(1, 2),
  nrow = 2
)

# Save
ggsave(here::here(fig_dir, "fig_1.png"),
       height = 20,
       width = 20,
       units = "cm")
ggsave(here::here(fig_dir, "fig_1.pdf"),
       height = 20,
       width = 20,
       units = "cm")

# ------------------------------------------------------------------------------
# Supplementary figures --------------------------------------------------------
# ------------------------------------------------------------------------------

# Figure S1 --------------------------------------------------------------------


# Define count cutoffs
low_ATAC <- 800
low_RNA <- 200

# Knee-plot for ATAC 
p_knee_atac <- tibble(nCount_ATAC = sort(dat_raw$nCount_ATAC,
                                         decreasing = TRUE)) |>
  rownames_to_column("rank") |>
  mutate(rank = as.numeric(rank)) |>
  head(100000) |>
  ggplot() +
  geom_point(aes(x = rank,
                 y = nCount_ATAC),
             color = "grey50",
             shape = 21) +
  geom_hline(aes(yintercept = low_ATAC),
             linetype = 2,
             color = "grey20",
             alpha = 0.8,
             linewidth = 1.1) +
  scale_y_continuous(trans = "log10", 
                     breaks = c(100, 1000, 10000, 100000),
                     labels = c("100", "1k", "10k", "100k"),
                     name = "ATAC counts") +
  scale_x_continuous(trans = "log10",
                     breaks = c(1, 10, 100, 1000, 10000, 100000),
                     labels = c("1", "10", "100", "1k", "10k", "100k"),
                     name = "Barcodes in rank order") +
  theme_minimal_grid()

# Knee-plot for RNA 
p_knee_rna <- tibble(nCount_RNA = sort(dat_raw$nCount_RNA,
                                       decreasing = TRUE)) |>
  rownames_to_column("rank") |>
  mutate(rank = as.numeric(rank)) |>
  head(100000) |>
  ggplot() +
  geom_point(aes(x = rank,
                 y = nCount_RNA),
             color = "grey50",
             shape = 21) +
  geom_hline(aes(yintercept = low_RNA),
             linetype = 2,
             color = "grey20",
             alpha = 0.8,
             linewidth = 1.1) +
  scale_y_continuous(trans = "log10", 
                     breaks = c(100, 1000, 10000, 100000),
                     labels = c("100", "1k", "10k", "100k"),
                     name = "RNA counts") +
  scale_x_continuous(trans = "log10",
                     breaks = c(1, 10, 100, 1000, 10000, 100000),
                     labels = c("1", "10", "100", "1k", "10k", "100k"),
                     name = "Barcodes in rank order") +
  theme_minimal_grid()


# Knee plots together
p_knees <- cowplot::plot_grid(
  p_knee_atac, 
  p_knee_rna,
  labels = c("A", "B"),
  nrow = 1
)

# Scatter plot
p_scatter <- dat_raw@meta.data |>
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
  arrange(cells_10x) |> 
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
  scale_color_manual(values = c("grey80", "#91AB73"),
                     name = element_blank()) +
  theme_minimal_grid() +
  theme(legend.position = "bottom")

# Number of cells per sample
p_cells_sample <- dat@meta.data |>
  group_by(sample_name) |>
  mutate(sample_name = factor(sample_name, levels = c(c("25C_Rep2",
                                                        "25C_Rep1",
                                                        "18C_Rep2",
                                                        "18C_Rep1")))) |>
  tally() |>
  ggplot(aes(x = n/10^3,
             y = sample_name,
             label = scales::comma_format()(n))) +
  geom_col(color = "grey20",
           fill = "grey80",
           width = 0.8)+
  scale_x_continuous(expand = c(0, 0),
                     name = "Number of cells",
                     labels = 
                       scales::label_number(scale = 1, suffix = "k")) +
  geom_label(nudge_x = -0.5,
             label.size = 0.1) +
  scale_y_discrete(name = element_blank(),
                   labels = c("25°C Rep 2",
                              "25°C Rep 1",
                              "18°C Rep 2",
                              "18°C Rep 1")) +
  cowplot::theme_minimal_vgrid() 

# Total number of ATAC reads
p_reads_atac <- dat@meta.data |>
  group_by(sample_name) |>
  mutate(sample_name = factor(sample_name, levels = c(c("25C_Rep2",
                                                        "25C_Rep1",
                                                        "18C_Rep2",
                                                        "18C_Rep1")))) |>
  summarize(n = sum(nCount_ATAC)) |>
  ggplot(aes(x = n/10^6,
             y = sample_name,
             label = scales::comma_format()(n))) +
  geom_col(color = "grey20",
           fill = "grey80",
           width = 0.8)+
  scale_x_continuous(expand = c(0, 0),
                     name = "Number of ATAC reads",
                     labels = 
                       scales::label_number(scale = 1, suffix = "M")) +
  geom_label(nudge_x = -5,
             label.size = 0.1) +
  scale_y_discrete(name = element_blank(),
                   labels = c("25°C Rep 2",
                              "25°C Rep 1",
                              "18°C Rep 2",
                              "18°C Rep 1")) +
  cowplot::theme_minimal_vgrid() 

# Total number of ATAC reads
p_reads_rna <- dat@meta.data |>
  group_by(sample_name) |>
  mutate(sample_name = factor(sample_name, levels = c(c("25C_Rep2",
                                                        "25C_Rep1",
                                                        "18C_Rep2",
                                                        "18C_Rep1")))) |>
  summarize(n = sum(nCount_RNA)) |>
  ggplot(aes(x = n/10^6,
             y = sample_name,
             label = scales::comma_format()(n))) +
  geom_col(color = "grey20",
           fill = "grey80",
           width = 0.8)+
  scale_x_continuous(expand = c(0, 0),
                     name = "Number of RNA reads",
                     labels = 
                       scales::label_number(scale = 1, suffix = "M")) +
  geom_label(nudge_x = -1,
             label.size = 0.1) +
  scale_y_discrete(name = element_blank(),
                   labels = c("25°C Rep 2",
                              "25°C Rep 1",
                              "18°C Rep 2",
                              "18°C Rep 1")) +
  cowplot::theme_minimal_vgrid() 

# Per sample metrics
p_samples <- cowplot::plot_grid(
  p_cells_sample,
  p_reads_atac,
  p_reads_rna,
  labels = c("D", "E", "F"),
  vjust = 1,
  nrow = 3
)

# Second row of plots
p_row_2 <- cowplot::plot_grid(
  p_scatter,
  p_samples,
  labels = c("C", ""),
  rel_widths = c(1.3, 1),
  nrow = 1
)

# All plots together
cowplot::plot_grid(
  p_knees, 
  p_row_2,
  rel_heights = c(1, 1.3),
  nrow = 2
)

# Save
ggsave(here::here(fig_dir, "fig_s1.png"),
       height = 25,
       width = 25,
       units = "cm")
ggsave(here::here(fig_dir, "fig_s1.pdf"),
       height = 25,
       width = 25,
       units = "cm")


# Figure S2 --------------------------------------------------------------------

# Count ATAC
p_c_atac <- dat@meta.data |>
  ggplot(aes(x = nCount_ATAC,
             y = sample_name_fig)) +
  geom_violin(color = "grey20",
              fill = "grey80",
              width = 0.8) +
  geom_boxplot(color = "grey20",
              fill = "grey90",
              width = 0.1,
              outlier.shape = NA) +
  scale_y_discrete(name = element_blank()) +
  scale_x_continuous(expand = c(0, 0),
                     labels = scales::comma_format(),
                     name = "ATAC reads per cell",
                     limits = c(0, 22000),
                     position = "top") +
  cowplot::theme_minimal_grid() 

# Feature ATAC
p_f_atac <- dat@meta.data |>
  ggplot(aes(x = nFeature_ATAC,
             y = sample_name_fig)) +
  geom_violin(color = "grey20",
              fill = "grey80",
              width = 0.8) +
  geom_boxplot(color = "grey20",
               fill = "grey90",
               width = 0.1,
               outlier.shape = NA) +
  scale_y_discrete(name = element_blank()) +
  scale_x_continuous(expand = c(0, 0),
                     labels = scales::comma_format(),
                     name = "ATAC features per cell",
                     limits = c(0, 6500),
                     position = "top") +
  cowplot::theme_minimal_grid() 

# Count RNA
p_c_rna <- dat@meta.data |>
  ggplot(aes(x = nCount_RNA,
             y = sample_name_fig)) +
  geom_violin(color = "grey20",
              fill = "grey80",
              width = 0.8) +
  geom_boxplot(color = "grey20",
               fill = "grey90",
               width = 0.1,
               outlier.shape = NA) +
  scale_y_discrete(name = element_blank()) +
  scale_x_continuous(expand = c(0, 0),
                     labels = scales::comma_format(),
                     name = "RNA reads per cell",
                     limits = c(0, 7000),
                     position = "top") +
  cowplot::theme_minimal_grid() 

# Feature RNA
p_f_rna <- dat@meta.data |>
  ggplot(aes(x = nFeature_RNA,
             y = sample_name_fig)) +
  geom_violin(color = "grey20",
              fill = "grey80",
              width = 0.8) +
  geom_boxplot(color = "grey20",
               fill = "grey90",
               width = 0.1,
               outlier.shape = NA) +
  scale_y_discrete(name = element_blank()) +
  scale_x_continuous(expand = c(0, 0),
                     labels = scales::comma_format(),
                     name = "Expressed genes detected per cell",
                     limits = c(0, 2250),
                     position = "top") +
  cowplot::theme_minimal_grid() 

# Count and feature plots together
p_c_f <- cowplot::plot_grid(
  p_c_rna,
  p_c_atac,
  p_f_rna,
  p_f_atac,
  nrow = 2,
  labels = c("A", "C", "B", "D")
)


# Plot transcription start site enrichment for each sample
p_tss <- TSSPlot(
  dat,
  assay = "ATAC",
  group.by = "sample_name_fig") +
  scale_color_manual(values = rep("grey50", 4)) +
  cowplot::theme_minimal_grid() +
  theme(legend.position = "none",
        plot.title = element_blank(),
        strip.background = element_rect(fill = "grey90"))

# Detailed RNA QC metrics -----

# Percent mito
p_mt <- dat@meta.data |>
  ggplot(aes(x = percent.mt/100,
             y = sample_name_fig)) +
  geom_violin(color = "grey20",
              fill = "grey80",
              width = 0.8) +
  geom_boxplot(color = "grey20",
               fill = "grey90",
               width = 0.1,
               outlier.shape = NA) +
  scale_y_discrete(name = element_blank()) +
  scale_x_continuous(expand = c(0, 0),
                     name = "RNA reads mapped to mitochondrial genome",
                     limits = c(0, .0525),
                     labels = scales::label_percent(),
                     position = "top") +
  cowplot::theme_minimal_grid() 

# Percent ribo
p_ribo <- dat@meta.data |>
  ggplot(aes(x = percent.ribo/100,
             y = sample_name_fig)) +
  geom_violin(color = "grey20",
              fill = "grey80",
              width = 0.8) +
  geom_boxplot(color = "grey20",
               fill = "grey90",
               width = 0.1,
               outlier.shape = NA) +
  scale_y_discrete(name = element_blank()) +
  scale_x_continuous(expand = c(0, 0),
                     name = "RNA reads from ribosomal genes",
                     limits = c(0, 0.27),
                     labels = scales::label_percent(),
                     position = "top") +
  cowplot::theme_minimal_grid() 

# Count and feature plots together
p_rna_qc <- cowplot::plot_grid(
  p_mt,
  p_ribo,
  nrow = 2,
  labels = c("F", "G"),
  vjust = 1
)

# Detailed ATAC QC metrics -----

# Nucleosome signal
p_nuc <- dat@meta.data |>
  ggplot(aes(x = nucleosome_signal,
             y = sample_name_fig)) +
  geom_violin(color = "grey20",
              fill = "grey80",
              width = 0.8) +
  geom_boxplot(color = "grey20",
               fill = "grey90",
               width = 0.1,
               outlier.shape = NA) +
  scale_y_discrete(name = element_blank()) +
  scale_x_continuous(expand = c(0, 0),
                     name = "Nucleosome signal",
                     limits = c(0, 1.4),
                     breaks = seq(from = 0, to = 1.25, by = 0.25),
                     position = "top") +
  cowplot::theme_minimal_grid() 

# FRiP 
p_frip <- dat@meta.data |>
  ggplot(aes(x = FRiP,
             y = sample_name_fig)) +
  geom_violin(color = "grey20",
              fill = "grey80",
              width = 0.8) +
  geom_boxplot(color = "grey20",
               fill = "grey90",
               width = 0.1,
               outlier.shape = NA) +
  scale_y_discrete(name = element_blank()) +
  scale_x_continuous(expand = c(0, 0),
                     name = "Fraction of reads in peaks (FRiP)",
                     limits = c(0, 1.05),
                     breaks = seq(from = 0, to = 1, by = 0.25),
                     position = "top") +
  cowplot::theme_minimal_grid() 


# Count and feature plots together
p_atac_qc <- cowplot::plot_grid(
  p_nuc,
  p_frip,
  nrow = 2,
  labels = c("H", "I"),
  vjust = 1
)

p_qc <- cowplot::plot_grid(
  p_rna_qc,
  p_atac_qc,
  nrow = 1
)

# Combine all plots
cowplot::plot_grid(
  p_c_f,
  p_tss,
  p_qc,
  nrow = 3,
  labels = c(NA, "E", NA)
)

# Save
ggsave(here::here(fig_dir, "fig_s2.png"),
       height = 30,
       width = 25,
       units = "cm")
ggsave(here::here(fig_dir, "fig_s2.pdf"),
       height = 30,
       width = 25,
       units = "cm")

