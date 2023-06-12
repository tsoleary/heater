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
  acc_temp = as.factor(c("18C", "18C", "25C", "25C"))
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
res <- results(dds, name = "acc_temp_25C_vs_18C")
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


# Quick volcano plot
p <- res |>
  as_tibble(rownames = "gene") |>
  ggplot(aes(label = gene)) +
  geom_hline(yintercept = -log10(0.05),
             color = "grey50",
             linetype = 2) +
  geom_point(aes(y = -log10(padj),
                 x = log2FoldChange,
                 fill = padj < 0.05,
                 size = padj < 0.05),
             color = "grey90",
             stroke = 0.2,
             shape = 21) +
  scale_size_manual(values = c(1, 2)) +
  scale_fill_manual(values = c("grey80", "firebrick")) +
  scale_y_continuous(expand = c(0, 0.5),
                     breaks = c(0, 5, 10, 15),
                     labels = c("0", "5", "10", "15")) +
  #ggh4x::coord_axes_inside(labels_inside = TRUE) +
  cowplot::theme_cowplot() +
  theme(legend.position = "none",
        axis.title.x = element_blank())

plotly::ggplotly(p, tooltip = c("gene"))
