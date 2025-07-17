# ------------------------------------------------------------------------------
# Compare DARs
# TS O'Leary
# ------------------------------------------------------------------------------

# Load libraries
library(tidyverse)

# Load data
dars_bulk <- readRDS(here::here(out_dir, "dars_bulk_MAST_exclusion.rds"))
dars_all <- readRDS(here::here(out_dir, "dars_cell-type_MAST_exclusion.rds"))

# Pseudobulk comparison --------------------------------------------------------


# Count --
# Number of detected regions per cell type
dars_t <- dars_bulk |> 
  filter(pct.18.all > 0 | pct.25.all > 0) |> 
  tally(name = "n_total") 

# Number of dars
dars_n <- dars_bulk |> 
  filter(p_val_adj.all < 0.05 & shared) |>
  tally(name = "n_dars") |> 
  arrange(desc(n_dars))

dars_02 <- dars_bulk |> 
  filter(pct.18.all > 0.02 | pct.25.all > 0.02) |> 
  tally(name = "n_02_pct")

dars_10 <- dars_bulk |> 
  filter(pct.18.all > 0.1 | pct.25.all > 0.1) |> 
  tally(name = "n_02_pct")

# Compile numbers
dars_n # 176
dars_02 # 90228
dars_10 # 53744
dars_t # 101344
dars_n/dars_t # 0.001736659 --- 0.17% dars 
dars_n/dars_02 # 0.001950614 --- 0.19% dars 

# Plot
dars_bulk |> 
  ggplot(aes(x = avg_log2FC_18_25.18.1,
             y = avg_log2FC_18_25.18.2)) +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  geom_point(shape = 21, 
             fill = "grey90", 
             color = "grey70",
             alpha = 0.8, size = 0.5) +
  geom_abline(slope = 1, intercept = 0,
              linetype = 2) +
  geom_smooth(data = dars_bulk |> 
                filter(p_val_adj.all < 0.05 & shared),
              method = "lm",
              color = "tomato3",
              fullrange = TRUE,
              se = FALSE,
              linewidth = 0.75,
              linetype = 3) +
  geom_point(data = dars_bulk |> 
               filter(p_val_adj.all < 0.05 & shared),
             color = "grey20",
             fill = "tomato3",
             stroke = 0.75,
             shape = 21) +
  ggpmisc::stat_poly_eq(data = dars_bulk |> 
                          filter(p_val_adj.all < 0.05 & shared),
                        formula = y ~ x,
                        aes(label = ..rr.label..),
                        parse = TRUE,
                        geom = "label",
                        label.x = 5.5,
                        label.y = 5.5,
                        alpha = 0.8,
                        size = 3,
                        color = "tomato3",
                        fill = "white") +
  lims(x = c(-lim_size, lim_size),
       y = c(-lim_size, lim_size)) +
  labs(x = expression("Log"[2]*"(fold-difference) 18°C Rep 1 only"),
       y = expression("Log"[2]*"(fold-difference) 18°C Rep 2 only")) + 
  theme_bw()

# Save 
ggsave(
  here::here("output/figs/sample_exclusion/dars_bulk_sample_exclusion.pdf"),
  width = 5,
  height = 5,
  units = "in"
)

# Cell-type specific comparison ------------------------------------------------

# Create limits for square plots
lim_size <- 3.5

# Scatter plot of all points and DARs
dars_all |> 
  ggplot(aes(x = avg_log2FC_18_25.18.1,
             y = avg_log2FC_18_25.18.2)) +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  geom_point(shape = 21, 
             fill = "grey90", 
             color = "grey70",
             alpha = 0.8, size = 0.5) +
  geom_abline(slope = 1, intercept = 0,
              linetype = 2) +
  geom_smooth(data = dars_bulk |> 
                filter(p_val_adj.all < 0.05 & shared),
              method = "lm",
              color = "tomato3",
              fullrange = TRUE,
              se = FALSE,
              linewidth = 0.75,
              linetype = 3) +
  geom_point(data = dars_all |> 
               filter(p_val_adj.all < 0.05 & shared),
             color = "grey20",
             fill = "tomato3",
             stroke = 0.25,
             shape = 21) +
  ggpmisc::stat_poly_eq(data = dars_bulk |> 
                          filter(p_val_adj.all < 0.05 & shared),
                        formula = y ~ x,
                        aes(label = ..rr.label..),
                        parse = TRUE,
                        geom = "label",
                        label.x = 2.5,
                        label.y = 2.5,
                        alpha = 0.8,
                        size = 3,
                        color = "tomato3",
                        fill = "white") +
  lims(x = c(-lim_size, lim_size),
       y = c(-lim_size, lim_size)) +
  labs(x = expression("Log"[2]*"(fold-difference) 18°C Rep 1 only"),
       y = expression("Log"[2]*"(fold-difference) 18°C Rep 2 only")) + 
  facet_wrap(~cell_type) +
  theme_bw()

# Save 
ggsave(
  here::here("output/figs/sample_exclusion/dars_sample_exclusion_all.pdf"),
  width = 6,
  height = 6,
  units = "in"
)

# Count the shared and non-shared 
p_only <- dars_all |> 
  filter(p_val_adj.all < 0.05) |> 
  mutate(quad_1_3 = case_when((avg_log2FC_18_25.18.1 > 0 & avg_log2FC_18_25 > 0) | 
                                (avg_log2FC_18_25.18.1 < 0 & avg_log2FC_18_25 < 0) ~ "shared",
                              .default = "off")) |> 
  group_by(quad_1_3) |> 
  tally() 


p_and_nuc <- dars_all |> 
  filter(p_val_adj.all < 0.05) |>
  filter(pct.18.all > 0.1 | pct.25.all > 0.1) |> 
  mutate(quad_1_3 = case_when((avg_log2FC_18_25.18.1 > 0 & avg_log2FC_18_25 > 0) | 
                                (avg_log2FC_18_25.18.1 < 0 & avg_log2FC_18_25 < 0) ~ "shared",
                              .default = "off")) |> 
  group_by(quad_1_3) |> 
  tally()

p_nuc_lfc <- dars_all |> 
  filter(p_val_adj.all < 0.05) |>
  filter(pct.18.all > 0.1 | pct.25.all > 0.1) |> 
  filter(abs(avg_log2FC_18_25.all) > 0.25) |> 
  mutate(quad_1_3 = case_when((avg_log2FC_18_25.18.1 > 0 & avg_log2FC_18_25 > 0) | 
                                (avg_log2FC_18_25.18.1 < 0 & avg_log2FC_18_25 < 0) ~ "shared",
                              .default = "off")) |> 
  group_by(quad_1_3) |> 
  tally() 

# Tally the number of DARs --- all shared response the thresholds only change 
#  the number of DARs
bind_rows(list("p-value only" = p_only, 
               "p-value & nuclei detection" = p_and_nuc, 
               "p-value, nuclei detection, & lfc" = p_nuc_lfc),
          .id = "filter")



