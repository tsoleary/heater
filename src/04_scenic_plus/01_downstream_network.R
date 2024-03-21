# ------------------------------------------------------------------------------
# SCENIC+ downstream crude network analysis
# TS O'Leary
# ------------------------------------------------------------------------------

# Load libraries
library(tidyverse)

# Load all cell-type specific acclimation comparisons
degs <- readRDS(here::here("output/degs/degs_cell-type_MAST.rds")) |> 
  mutate(sig = ifelse(abs(avg_log2FC_18_25) >= 0.25 & p_val_adj < 0.05 &
                        (pct.18 >= 0.10 | pct.25 >= 0.10), 
                      TRUE, FALSE))
dars <- readRDS(here::here("output/dars/dars_cell-type_MAST.rds")) |> 
  mutate(sig = ifelse(abs(avg_log2FC_18_25) >= 0.25 & p_val_adj < 0.05 &
                        (pct.18 >= 0.10 | pct.25 >= 0.10), 
                      TRUE, FALSE)) 

degs <- readRDS(here::here("output/degs/degs_cell-type.rds")) |> 
  mutate(sig = ifelse(abs(avg_log2FC_18_25) >= 0.25 & p_val_adj < 0.05 &
                        (pct.18 >= 0.10 | pct.25 >= 0.10), 
                      TRUE, FALSE))
dars <- readRDS(here::here("output/dars/dars_cell-type.rds")) |> 
  mutate(sig = ifelse(abs(avg_log2FC_18_25) >= 0.25 & p_val_adj < 0.05 &
                        (pct.18 >= 0.10 | pct.25 >= 0.10), 
                      TRUE, FALSE)) 


# eRegulon data
eReg <- read_csv(here::here("data/processed/scenic/eRegulon_data.csv")) |> 
  mutate(eRegulon_type = factor(
    str_extract(Region_signature_name, "(?<=_)\\+|-(?=_)"), levels = c("+", "-"))) |> 
  mutate(eRegulon_type = ifelse(eRegulon_type == "+", "activator", "repressor"))


# Part I --- Differentially expressed TFs --------------------------------------
deTFs <- degs |> 
  filter(sig) |> 
  filter(gene %in% unique(eReg$TF)) |> 
  mutate(tf_cell = paste(gene, cell_type, sep = " | "))

# Part II --- Target peak accessibility LFC distribution -----------------------
target_peaks <- eReg |> 
  filter(TF %in% deTFs$gene) |> 
  select(TF, Region_signature_name, eRegulon_type, Region) |> 
  distinct() |> 
  mutate(Region = str_remove_all(Region, "chr")) |> 
  mutate(Region = str_replace_all(Region, ":", "-")) 

# Create target peak data frame
target_peak_df <- dars |> 
  left_join(target_peaks,
            by = c("region" = "Region"),
            relationship = "many-to-many") |> 
  filter(!is.na(TF)) |> 
  mutate(tf_cell = paste(TF, cell_type, sep = " | ")) |> 
  filter(tf_cell %in% deTFs$tf_cell) |> 
  left_join(deTFs, by = c("tf_cell"), suffix = c("", ".TF"))

# One-sample t-test if the peak LFC distribution for each is non-zero
target_peak_t <- target_peak_df |>
  group_by(cell_type, Region_signature_name, TF) |> 
  filter(pct.18 >= 0.1 | pct.25 >= 0.1) |> 
  mutate(
    alternative = 
      case_when(avg_log2FC_18_25.TF > 0 & eRegulon_type == "activator" ~ "greater",
                avg_log2FC_18_25.TF < 0 & eRegulon_type == "activator" ~ "less",
                avg_log2FC_18_25.TF > 0 & eRegulon_type == "repressor" ~ "less",
                avg_log2FC_18_25.TF < 0 & eRegulon_type == "repressor" ~ "greater")) |> 
  nest() |> 
  mutate(t_test = map(data,
                      ~t.test(.$avg_log2FC_18_25,
                              alternative = unique(.$alternative)))) |> 
  mutate(t_test = map(t_test, broom::tidy)) |> 
  unnest(t_test) |> 
  mutate(tf_cell = paste(TF, cell_type, sep = " | "))

# Part III - Target gene expression LFC distribution ---------------------------
target_genes <- eReg |> 
  filter(TF %in% deTFs$gene) |> 
  select(TF, Region_signature_name, eRegulon_type, Gene) |> 
  distinct()

# Create target gene data frame
target_gene_df <- degs |> 
  left_join(target_genes,
            by = c("gene" = "Gene"),
            relationship = "many-to-many") |> 
  filter(!is.na(TF)) |> 
  mutate(tf_cell = paste(TF, cell_type, sep = " | ")) |>
  filter(tf_cell %in% deTFs$tf_cell) |> 
  left_join(deTFs, by = c("tf_cell"), suffix = c("", ".TF"))

# One-sample t-test if the gene LFC distribution for each is non-zero
target_gene_t <- target_gene_df |>
  group_by(cell_type, Region_signature_name, TF) |> 
  filter(pct.18 >= 0.1 | pct.25 >= 0.1) |> 
  mutate(
    alternative = 
      case_when(avg_log2FC_18_25.TF > 0 & eRegulon_type == "activator" ~ 
                  "greater",
                avg_log2FC_18_25.TF < 0 & eRegulon_type == "activator" ~ "less",
                avg_log2FC_18_25.TF > 0 & eRegulon_type == "repressor" ~ "less",
                avg_log2FC_18_25.TF < 0 & eRegulon_type == "repressor" ~ 
                  "greater")) |> 
  nest() |> 
  mutate(t_test = map(data,
                      ~t.test(.$avg_log2FC_18_25,
                              alternative = unique(.$alternative)))) |> 
  mutate(t_test = map(t_test, broom::tidy)) |> 
  unnest(t_test) |> 
  mutate(tf_cell = paste(TF, cell_type, sep = " | "))


# The GRNs that pass all three tests -------------------------------------------
# Compile stats from tests into one data frame and adjust p-values
GRNs <- target_gene_t |> 
  full_join(target_peak_t, by = c("cell_type", "TF", "Region_signature_name"),
            suffix = c(".gene", ".peak")) |> 
  full_join(deTFs, by = c("cell_type", "TF" = "gene")) |> 
  mutate(p.value.adj.peak = p.adjust(p.value.peak, method = "BH", n = 26),
         p.value.adj.gene = p.adjust(p.value.gene, method = "BH", n = 26))

# Print out GRNs that pass three filter
GRNs |> 
  select(Region_signature_name, avg_log2FC_18_25, estimate.peak, 
         estimate.gene, p_val_adj, p.value.adj.peak, p.value.adj.gene) |> 
  filter(Region_signature_name %in% unique(dotplot_dat$Region_signature_name)) |> 
  filter(p_val_adj < 0.05 & p.value.adj.peak < 0.05 & p.value.adj.gene < 0.05)

diffGRN <- GRNs |> 
  select(Region_signature_name, avg_log2FC_18_25, estimate.peak, 
         estimate.gene, p_val_adj, p.value.peak, p.value.gene) |> 
  filter(p_val_adj < 0.05 & p.value.peak < 0.05 & p.value.gene < 0.05) |> 
  filter(Region_signature_name %in% unique(dotplot_dat$Region_signature_name)) |> 
  mutate(cell_type_eGRN = paste(Region_signature_name, cell_type, sep = "\n")) |> 
  pull(cell_type_eGRN)



# Print out GRNs that pass three filter
GRNs |> 
  select(Region_signature_name, avg_log2FC_18_25, estimate.peak, 
         estimate.gene, p_val_adj, p.value.peak, p.value.gene) |> 
  filter(p_val_adj < 0.05 & p.value.peak < 0.05 & p.value.gene < 0.05)

# Sara's question about TF not reaching target peaks
# TF is DiffExp, but target peaks are not or are in the other direction
GRNs |> 
  select(Region_signature_name, avg_log2FC_18_25, estimate.peak, 
         estimate.gene, p_val_adj, p.value.adj.peak, p.value.adj.gene) |> 
  filter(p_val_adj < 0.05 & p.value.adj.peak > 0.05)

# Opposite direction
GRNs |> 
  select(Region_signature_name, avg_log2FC_18_25, estimate.peak, 
         estimate.gene, p_val_adj, p.value.adj.peak, p.value.adj.gene) |> 
  filter(p_val_adj < 0.05 & p.value.adj.peak > 0.05) |> 
  filter(avg_log2FC_18_25*estimate.peak < 0)
# ^ This actually seems more common in the 25°C to have a high gene expresion, 
# but target motif is not differentially accessible...

# Sara's question extended – limitation of transcription of target genes
GRNs |> 
  select(Region_signature_name, avg_log2FC_18_25, estimate.peak, 
         estimate.gene, p_val_adj, p.value.adj.peak, p.value.adj.gene) |> 
  filter(p_val_adj < 0.05 & p.value.adj.peak < 0.05 & p.value.adj.gene > 0.05)
# And for this: the distribution of LFC of gene expression matches TF and peaks 
# So this doesn't really match the case of incomplete compensation/limited
# transcriptional apparatus