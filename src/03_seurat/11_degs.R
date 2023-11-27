# ------------------------------------------------------------------------------
# Find differentially expressed genes between conditions
# TS O'Leary
# ------------------------------------------------------------------------------

# Load libraries
library(tidyverse)
library(Seurat)

# Load data
dat <- readRDS(here::here("data/processed/seurat_object/09_dat_annot.rds"))

# Output dir
out_dir <- "output/degs"

# Set default assay
DefaultAssay(dat) <- "SCT"

# Differential expression bulk between two conditions regardless of cluster ----
Idents(dat) <- "acc_temp"
degs <- FindMarkers(dat, 
                    ident.1 = "18°C",
                    ident.2 = "25°C")

# Filter only padj < 0.05
degs |>
  filter(p_val_adj < 0.05)

# Save degs regardless of cluster
saveRDS(degs, here::here(out_dir, "degs.rds"))

# Differential expression between conditions within clusters -------------------

# Create combined groups of clusters and acclimation temperatures
dat$celltype_acc <- paste(dat$cell_type, dat$acc_temp, sep = "_")

# Set that new variable as the ident
Idents(dat) <- "celltype_acc"

# Find markers genes between different clusters 
degs <- list()

# Loop through each cluster
for (i in unique(dat$cell_type)) {
  degs[[i]] <- FindMarkers(
    dat, 
    ident.1 = paste0(i, "_18°C"),
    ident.2 = paste0(i, "_25°C"),
    logfc.threshold = 0,
    min.pct = 0.02) |>
    rownames_to_column("gene")
}

# Combine all clusters together
degs <- bind_rows(degs, .id = "cell_type")

# Rename cols for later and adjust p-vals again to account for all clusters 
# being tested at the same time groups
degs <- degs |> 
  dplyr::rename("avg_log2FC_18_25" = "avg_log2FC",
         "pct.18" = "pct.1",
         "pct.25" = "pct.2") |> 
  mutate(padj = p.adjust(p_val_adj, method = "BH"))

# Save the degs cell types
saveRDS(degs, here::here(out_dir, "degs_cell-type_lfc_0_min.pct_0.02.rds"))

# Total number of DEGs
degs |>
  filter(padj < 0.05) |>
  tally()

# Number of DEGs per cell-type -----
degs |>
  group_by(cell_type) |> 
  add_tally(name = "n_genes") |> 
  filter(padj < 0.05) |>
  group_by(cell_type) |>
  add_tally() |> 
  mutate(prop_deg = n/n_genes) |> 
  distinct(cell_type, .keep_all = TRUE) |> 
  select(cell_type, n , n_genes, prop_deg) |> 
  arrange(desc(prop_deg))

# Pseudobulk analysis with DESeq2 ----------------------------------------------

# Load libraries
library(DESeq2)

# Aggregate samples data into counts matrix
counts <- AggregateExpression(dat, 
                              group.by = "sample_name",
                              assays = "RNA",
                              slot = "counts",
                              return.seurat = FALSE)[[1]]

# Create metadata for counts matrix
metadata <- tibble(
  samples = c("18C_Rep1","18C_Rep2","25C_Rep1","25C_Rep1"),
  acc_temp = factor(c("18C", "18C", "25C", "25C"), 
                    levels = c("25C", "18C"))
)

# Create DESeq Data Set
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = metadata,
                              design = ~ acc_temp)

# Filter out genes with low counts across the four samples
keep_genes <- rowSums(counts(dds)) >= 10
dds <- dds[keep_genes,]

# Run DESeq2
dds <- DESeq(dds)

# Generate results object
res <- results(dds, name = "acc_temp_18C_vs_25C")
deg_res <- res |>
  as_tibble(rownames = "gene")
degs <- res |>
  as_tibble(rownames = "gene") |>
  dplyr::filter(padj < 0.05) |>
  arrange(padj)

# Save results
saveRDS(res, here::here(out_dir, "pseudobulk_DESeq_res.rds"))
saveRDS(degs, here::here(out_dir, "pseudobulk_DESeq.rds"))

# Save csv with only gene names for GO analysis
write_csv(
  degs |> 
    select(gene), 
  here::here(out_dir, "pseudobulk_DESeq_gene_names.csv")
)


