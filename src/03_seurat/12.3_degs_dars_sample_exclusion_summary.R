# ------------------------------------------------------------------------------
# Sample exclusion combine
# TS O'Leary
# ------------------------------------------------------------------------------

# Load libraries
library(tidyverse)

# Cell-type specific data
degs <- readRDS(here::here("output/degs/degs_cell-type_MAST.rds"))
degs_18_Rep1 <- readRDS(
  here::here("output/degs/sample_exclusion/degs_cell-type_MAST_18_Rep1.rds")
)
degs_18_Rep2 <- readRDS(
  here::here("output/degs/sample_exclusion/degs_cell-type_MAST_18_Rep2.rds")
)
dars <- readRDS(here::here("output/dars/dars_cell-type_MAST.rds"))
dars_18_Rep1 <- readRDS(
  here::here("output/dars/sample_exclusion/dars_cell-type_MAST_18_Rep1.rds")
)
dars_18_Rep2 <- readRDS(
  here::here("output/dars/sample_exclusion/dars_cell-type_MAST_18_Rep2.rds")
)

# Combine all for comparison
degs_all <- full_join(
  degs, degs_18_Rep1,
  by = c("cell_type", "gene"),
  suffix = c("", ".18.1")) |>
  full_join(
    degs_18_Rep2,
    by = c("cell_type", "gene"),
    suffix = c(".all", ".18.2")) |>
  mutate(shared = case_when(
    (avg_log2FC_18_25.18.1 > 0 & avg_log2FC_18_25.18.2 > 0) |
      (avg_log2FC_18_25.18.1 < 0 & avg_log2FC_18_25.18.2 < 0) ~ TRUE,
    .default = FALSE)) |> 
  mutate(tested = ifelse((pct.18.all >= 0.1 | pct.25.all >= 0.1) ,
                         TRUE, FALSE)) |> 
  mutate(sig = ifelse(p_val_adj.all < 0.05 & 
                        shared & tested & abs(avg_log2FC_18_25.all) > 0.25, 
                      TRUE, FALSE))

# Save 
saveRDS(
  degs_all,
  here::here("output/degs/degs_all.rds")
)


# Count the number of genes tested and DEGs detected in each cell_type
degs_tally <- degs_all |> 
  group_by(cell_type) |> 
  filter(tested) |> 
  add_tally(name = "n_tested_genes") |> 
  filter(sig) |> 
  add_tally(name = "n_degs") |> 
  select(cell_type, n_tested_genes, n_degs) |> 
  distinct(cell_type, .keep_all = TRUE) |> 
  mutate(percent_degs = (n_degs/n_tested_genes)*100) |> 
  arrange(desc(percent_degs))

# Save 
write_csv(
  degs_tally,
  here::here("output/degs/degs_tally.csv")
)

# Combine all for comparison
dars_all <- full_join(
  dars, dars_18_Rep1,
  by = c("cell_type", "region"),
  suffix = c("", ".18.1")) |>
  full_join(
    dars_18_Rep2,
    by = c("cell_type", "region"),
    suffix = c(".all", ".18.2")) |>
  mutate(shared = case_when(
    (avg_log2FC_18_25.18.1 > 0 & avg_log2FC_18_25.18.2 > 0) |
      (avg_log2FC_18_25.18.1 < 0 & avg_log2FC_18_25.18.2 < 0) ~ TRUE,
    .default = FALSE)) |> 
  mutate(tested = ifelse((pct.18.all >= 0.1 | pct.25.all >= 0.1) ,
                         TRUE, FALSE)) |> 
  mutate(sig = ifelse(p_val_adj.all < 0.05 & 
                        shared & tested & abs(avg_log2FC_18_25.all) > 0.25, 
                      TRUE, FALSE))

# Save 
saveRDS(
  dars_all,
  here::here("output/dars/dars_all.rds")
)

# Count the number of regions tested and dars detected in each cell_type
dars_tally <- dars_all |> 
  group_by(cell_type) |> 
  filter(tested) |> 
  add_tally(name = "n_tested_regions") |> 
  group_by(sig, cell_type) |> 
  add_tally(name = "n_dars") |> 
  filter(sig | cell_type == "amnioserosa") |>
  mutate(n_dars = ifelse(sig, n_dars, 0)) |> 
  ungroup(sig) |> 
  select(-sig) |> 
  select(cell_type, n_tested_regions, n_dars) |> 
  distinct(cell_type, .keep_all = TRUE) |> 
  mutate(percent_dars = (n_dars/n_tested_regions)*100) |> 
  arrange(desc(percent_dars))

# Save 
write_csv(
  dars_tally,
  here::here("output/dars/dars_tally.csv")
)
