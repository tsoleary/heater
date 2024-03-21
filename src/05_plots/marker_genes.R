# ------------------------------------------------------------------------------
# Create plots for marker genes
# TS O'Leary
# ------------------------------------------------------------------------------

# Load libraries
library(tidyverse)
library(Seurat)
library(Signac)

# Load plot themes 
source(here::here("src/04_plots/00_plot_themes.R"))

# Output fig dir
fig_dir <- "output/figs/markers"

# Load data
dat <- readRDS(here::here("data/processed/seurat_object/10_dat_linked.rds"))
markers <- read_csv(here::here("output/markers/cluster_markers_all.csv"))

# Set resolution that was used for marker gene detection
Idents(dat) <- "wknn_res.0.7"

DefaultAssay(dat) <- "SCT"

# Top three markers for each cluster
top_markers <- markers |> 
  group_by(cluster) |> 
  slice_min(max_pval, n = 3, with_ties = FALSE) |> 
  mutate(cluster = factor(cluster, levels = as.character(0:15))) |> 
  arrange(cluster)

# Output fig dir
fig_dir <- "output/figs/marker_genes"

# Percent nuclei with expressed gene
pct_exp <- FetchData(dat, 
          vars = c(top_markers$gene, "seurat_clusters"), 
          layer = "counts") |> 
  pivot_longer(cols = -c("seurat_clusters"),
               names_to = "gene",
               values_to = "count") |> 
  mutate(exp = count > 0) |> 
  group_by(seurat_clusters, gene) |> 
  summarize(pct_exp = sum(exp)/n()*100)

# Average expression 
avg_exp <- AverageExpression(dat)$SCT |> 
  as.data.frame() |> 
  rownames_to_column("gene") |> 
  filter(gene %in% top_markers$gene) |> 
  pivot_longer(cols = -c(gene),
               names_to = "cluster",
               values_to = "expression") |> 
  group_by(gene) |> 
  mutate(norm_expression = scale(expression))

# Combine the average expression with percent expressed data
df <- full_join(
  avg_exp, 
  pct_exp,
  by = c("cluster" = "seurat_clusters", "gene")
)

# Get a pleasing gene order for plotting
gene_order <- df |> 
  mutate(cluster = as.numeric(cluster)) |> 
  group_by(cluster) |> 
  arrange(cluster, desc(norm_expression)) |> 
  slice_max(norm_expression, n = 3) |> 
  pull(gene)

# Marker gene dot plot
df |> 
  mutate(cluster = factor(cluster, levels = as.character(0:15))) |> 
  mutate(gene = factor(gene, levels = top_markers$gene)) |> 
  ggplot(aes(x = gene, 
             y = cluster, 
             fill = norm_expression,
             size = pct_exp)) +
  geom_point(shape = 21, stroke = 0.25) +
  scale_fill_gradient(low = "grey95",
                      high = "firebrick",
                      name = "Normalized\nExpression") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 0.98),
        axis.title = element_blank(),
        strip.background = element_rect(fill = "grey90", color = "grey90")) +
  force_panelsizes(cols = c(1, 0.4)) +
  scale_size_continuous(range = c(0.25, 6),
                        labels = scales::label_percent(scale = 1),
                        name = "Nuclei\nexpressing\nthe gene")

# Gene example -----------------------------------------------------------------

# Set assay for expression data and clusters
DefaultAssay(dat) <- "SCT"

# Loop through all top genes
for (i in 1:nrow(top_markers)) {
  gene <- top_markers$gene[i]
  
  # Violin Plot
  p1 <- VlnPlot(dat,
                top_markers$gene[i],
                slot = "counts",
                cols = acc_colors,
                split.by = "acc_temp", 
                pt.size = 0) +
    labs(
      title = 
        paste0(top_markers$gene[i], 
               ": Cluster ", 
               top_markers$cluster[i],
               " marker")) +
    scale_y_continuous(expand = c(0, 0)) +
    scale_x_discrete(name = "Cluster") +
    cowplot::theme_minimal_hgrid()
  
  # Expression level on top of the UMAP projection
  p2 <- FeaturePlot(dat, 
                    top_markers$gene[i],
                    cols = c("grey85", "darkgreen"),
                    reduction = "umap",
                    keep.scale = "all",
                    split.by = "acc_temp") &
    theme_void() &
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5, face = "bold")) 
  
  # Combine the two plots
  cowplot::plot_grid(
    rel_heights = c(1, 1.25),
    p1, p2,
    nrow = 2
  )
  
  #Save each plot
  ggsave(
    here::here(
      fig_dir,
      paste0(top_markers$cluster[i],
             "_cluster_", 
             top_markers$gene[i],
             ".png")),
         height = 30,
         width = 25,
         units = "cm")
}
