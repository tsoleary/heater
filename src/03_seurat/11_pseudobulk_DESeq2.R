# ------------------------------------------------------------------------------
# Pseudobulk analysis
# June 08, 2023
# TS O'Leary
# ------------------------------------------------------------------------------

# Load libraries
require(tidyverse)
require(Seurat)
require(Signac)
require(DESeq2)

# Load data
dat <- readRDS(
  here::here("data/processed/seurat_object/06_dat_qc.rds")
)

# Aggregate samples data into counts matrix
counts <- AggregateExpression(dat, 
                              group.by = "sample_name",
                              assays = "RNA",
                              slot = "counts",
                              return.seurat = FALSE)[[1]]

# Create metadata for counts matrix
metadata <- tibble(
  samples = c("18C_Rep1","18C_Rep2","25C_Rep1","25C_Rep1"),
<<<<<<< HEAD
  acc_temp = factor(c("18C", "18C", "25C", "25C"), 
                    levels = c("25C", "18C"))
=======
  acc_temp = as.factor(c("18C", "18C", "25C", "25C"))
>>>>>>> a8de26aef8a500cdbd82fa277045d3e778c75e21
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
<<<<<<< HEAD
res <- results(dds, name = "acc_temp_18C_vs_25C")
=======
res <- results(dds, name = "acc_temp_25C_vs_18C")
>>>>>>> a8de26aef8a500cdbd82fa277045d3e778c75e21
deg_res <- res |>
  as_tibble(rownames = "gene")
degs <- res |>
  as_tibble(rownames = "gene") |>
  dplyr::filter(padj< 0.05) |>
  arrange(padj)

# Save results
saveRDS(res, here::here("output/degs/pseudobulk_DESeq_res.rds"))
saveRDS(degs, here::here("output/degs/pseudobulk_DESeq.rds"))
write_csv(degs |> select(gene), here::here("output/degs/pseudobulk_DESeq_gene_names.csv"))
