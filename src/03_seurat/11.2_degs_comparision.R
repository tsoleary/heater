# ------------------------------------------------------------------------------
# Compare degs
# TS O'Leary
# ------------------------------------------------------------------------------

# Load libraries
library(tidyverse)

# Load data
degs_bulk <- readRDS(here::here(out_dir, "degs_bulk_MAST_exclusion.rds"))
degs_all <- readRDS(here::here(out_dir, "degs_cell-type_MAST_exclusion.rds"))

# Pseudobulk comparison --------------------------------------------------------

# Count --
# DEGs
# Number of detected genes per cell type
degs_t <- degs_bulk |> 
  filter(pct.1.all > 0 | pct.2.all > 0) |> 
  tally(name = "n_total") 

# Number of DEGs
degs_n <- degs_bulk |> 
  filter(p_val_adj.all < 0.05 & shared) |>
  tally(name = "n_degs") |> 
  arrange(desc(n_degs))

degs_02 <- degs_bulk |> 
  filter(pct.1.all > 0.02 | pct.2.all > 0.02) |> 
  tally(name = "n_02_pct")

# Compile numbers
degs_n # 508
degs_02 # 4555
degs_t # 8982
degs_n/degs_t # 0.05655756 --- 5.7% DEGs 
degs_n/degs_02 # 0.05655756 --- 11.2% DEGs 

# Plot
degs_bulk |> 
  ggplot(aes(x = avg_log2FC.18.1,
             y = avg_log2FC.18.2)) +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  geom_point(shape = 21, 
             fill = "grey90", 
             color = "grey70",
             alpha = 0.8, size = 0.5) +
  geom_abline(slope = 1, intercept = 0,
              linetype = 2) +
  geom_smooth(data = degs_bulk |> 
                filter(p_val_adj.all < 0.05 & shared),
              method = "lm",
              color = "tomato3",
              fullrange = TRUE,
              se = FALSE,
              linewidth = 0.75,
              linetype = 3) +
  geom_point(data = degs_bulk |> 
               filter(p_val_adj.all < 0.05 & shared),
             color = "grey20",
             fill = "tomato3",
             stroke = 0.75,
             shape = 21) +
  ggpmisc::stat_poly_eq(data = degs_all |> 
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
  here::here("output/figs/sample_exclusion/degs_bulk_sample_exclusion.pdf"),
  width = 5,
  height = 5,
  units = "in"
)



# Cell-type specific comparison ------------------------------------------------

# Combine all for comparison


# DEGs
# Number of detected genes per cell type
degs_t <- degs_all |> 
  filter(pct.18.all > 0 | pct.25.all > 0) |> 
  group_by(cell_type) |> 
  tally(name = "n_total") 

# Number of DEGs
degs_n <- degs_all |> 
  filter(p_val_adj.all < 0.05 & shared) |>
  group_by(cell_type) |> 
  tally(name = "n_degs") |> 
  arrange(desc(n_degs))

degs_02 <- degs_all |> 
  filter(pct.18.all > 0.02 | pct.25.all > 0.02) |> 
  group_by(cell_type) |> 
  tally(name = "n_02_pct")

# Full join
degs_tally <- full_join(degs_n, degs_02) |> 
  full_join(degs_t) |> 
  mutate(p_degs = n_degs/n_total) |> 
  mutate(p_degs_tested = n_degs/n_02_pct)

write_csv(degs_tally, here::here("output/degs/degs_tally.csv"))

# Plots ------------------------------------------------------------------------

# Create limits for square plots
lim_size <- 6.6

# Scatter plot of all points and DEGs
degs_all |> 
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
  geom_smooth(data = degs_all |> 
                filter(p_val_adj.all < 0.05 & shared),
              method = "lm",
              color = "tomato3",
              fullrange = TRUE,
              se = FALSE,
              linewidth = 0.75,
              linetype = 3) +
  geom_point(data = degs_all |> 
               filter(p_val_adj.all < 0.05 & shared),
             color = "grey20",
             fill = "tomato3",
             stroke = 0.75,
             shape = 21) +
  ggpmisc::stat_poly_eq(data = degs_all |> 
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
  labs(x = expression("Log"[2]*"(fold-difference) – 18°C Rep 1 only"),
       y = expression("Log"[2]*"(fold-difference) – 18°C Rep 2 only")) + 
  facet_wrap(~cell_type) +
  theme_bw()

# Save 
ggsave(
  here::here("output/figs/sample_exclusion/degs_sample_exclusion_all.pdf"),
  width = 6,
  height = 6,
  units = "in"
)

# # Count the shared and non-shared 
# p_only <- degs_all |> 
#   filter(p_val_adj.all < 0.05) |> 
#   mutate(quad_1_3 = case_when((avg_log2FC_18_25.18.1 > 0 & avg_log2FC_18_25 > 0) | 
#                                 (avg_log2FC_18_25.18.1 < 0 & avg_log2FC_18_25 < 0) ~ "shared",
#                               .default = "off")) |> 
#   group_by(quad_1_3) |> 
#   tally() |> 
#   pivot_wider(names_from = quad_1_3, values_from = n) |> 
#   mutate(total_degs = shared + off) |> 
#   mutate(prop_off = off/total_degs)
# 
# 
# p_and_nuc <- degs_all |> 
#   filter(p_val_adj.all < 0.05) |>
#   filter(pct.18.all > 0.1 | pct.25.all > 0.1) |> 
#   mutate(quad_1_3 = case_when((avg_log2FC_18_25.18.1 > 0 & avg_log2FC_18_25 > 0) | 
#                                 (avg_log2FC_18_25.18.1 < 0 & avg_log2FC_18_25 < 0) ~ "shared",
#                               .default = "off")) |> 
#   group_by(quad_1_3) |> 
#   tally() |> 
#   pivot_wider(names_from = quad_1_3, values_from = n) |> 
#   mutate(total_degs = shared + off) |> 
#   mutate(prop_off = off/total_degs)
# 
# p_nuc_lfc <- degs_all |> 
#   filter(p_val_adj.all < 0.05) |>
#   filter(pct.18.all > 0.1 | pct.25.all > 0.1) |> 
#   filter(abs(avg_log2FC_18_25.all) > 0.25) |> 
#   mutate(quad_1_3 = case_when((avg_log2FC_18_25.18.1 > 0 & avg_log2FC_18_25 > 0) | 
#                                 (avg_log2FC_18_25.18.1 < 0 & avg_log2FC_18_25 < 0) ~ "shared",
#                               .default = "off")) |> 
#   group_by(quad_1_3) |> 
#   tally() |> 
#   pivot_wider(names_from = quad_1_3, values_from = n) |> 
#   mutate(total_degs = shared + off) |> 
#   mutate(prop_off = off/total_degs)
# 
# 
# bind_rows(list("p-value only" = p_only, 
#                "p-value & nuclei detection" = p_and_nuc, 
#                "p-value, nuclei detection, & lfc" = p_nuc_lfc),
#           .id = "filter") |> 
#   select(filter, total_degs, prop_off, shared, off)
