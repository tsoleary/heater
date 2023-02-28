# ------------------------------------------------------------------------------
# Differential expression b/w clusters in conditions
# November 04, 2022
# TS O'Leary
# ------------------------------------------------------------------------------

# Load libraries
require(Seurat)
require(tidyverse)

# Load data
dat_comb <- readRDS(here::here("calderon/data/dat_comb.rds"))

# Create a Cluster-Condition group
dat_comb$clust_cond <- paste("cluster",
                             Idents(dat_comb), 
                             dat_comb$time, 
                             "hrs",
                             sep = "_")

# Set the above groups as the idents for comparison
Idents(dat_comb) <- "clust_cond"


# Parameters for differential expression testing
min.pct <- 0.02
logfc.threshold <- 0.05


# Initialize list to save DEGs for each cluster treatment comparison
degs_all <- list()

# Loop through all clusters to find markers in each cluster
for (i in levels(dat_comb$seurat_clusters)) {
  
  # Print in the beginning
  print(paste("Begin loop", i))
  
  # Set up comparison groups
  group_5 <- paste("cluster", i, "5", "hrs", sep = "_")
  group_10 <- paste("cluster", i, "10", "hrs", sep = "_")
  
  # Find the conserved markers between 5 hr and 10 hr embryos
  marker <- FindMarkers(dat_comb,
                        ident.1 = group_5,
                        ident.2 = group_10,
                        group.by = "clust_cond",
                        min.pct = min.pct,
                        logfc.threshold = logfc.threshold)
  
  # Save markers in the list
  degs_all[[i]] <- marker %>% 
    rownames_to_column("gene")
  
  # Print in the end
  print(paste("End loop", i))
}

# Combine into a single data_frame
dat_comb_cluster_deg <- bind_rows(degs_all, 
                                  .id = "cluster")

saveRDS(dat_comb_cluster_deg, here::here("calderon/data/dat_comb_cluster_deg.rds"))

# Tally total number of total_cells in each cluster
cell_tally <- dat_comb@meta.data %>%
  group_by(seurat_clusters) %>%
  tally() %>%
  dplyr::rename(total_cells = n,
                cluster = seurat_clusters)

cell_time_tally <- dat_comb@meta.data %>%
  group_by(seurat_clusters, time) %>%
  tally() %>%
  pivot_wider(values_from = n,
              names_from = time) %>%
  dplyr::rename(cells_5 = `5`,
                cells_10 = `10`,
                cluster = seurat_clusters)

# Tally total number of genes tested in each category
gene_tally <- dat_comb_cluster_deg %>%
  group_by(cluster) %>%
  tally() %>%
  dplyr::rename(genes = n)

# Only with adjusted p-vals less than 0.05
deg_tally <- dat_comb_cluster_deg %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  group_by(cluster) %>%
  tally() %>%
  dplyr::rename(degs = n)

# Join these all together
df <- full_join(cell_tally, 
                full_join(gene_tally, deg_tally))

df <- full_join(cell_time_tally, df)


df %>%
  ggplot(aes(x = total_cells, y = genes)) +
  geom_point()

df %>%
  ggplot() +
  geom_point(aes(x = total_cells,
                 y = genes))

df %>%
  ggplot() +
  geom_point(aes(x = total_cells,
                 y = degs))


df %>%
  ggplot(aes(x = total_cells, y = degs/genes,
             size = genes)) +
  geom_point()

# Summary of nFeature & nCount per cell
nCounts <- dat_comb@meta.data %>%
  group_by(seurat_clusters) %>%
  summarise(nFeature_RNA = median(nFeature_RNA),
            nCount_RNA = median(nCount_RNA)) %>%
  dplyr::rename(cluster = seurat_clusters)


df <- full_join(df, nCounts)

# Set the theme for plotting after
theme_set(theme_bw() +
            theme(panel.border = element_blank(),
                  axis.line = element_line(color = "grey20")))


df %>%
  ggplot(aes(x = nFeature_RNA, 
             y = genes)) +
  geom_point()

# Average feature and average counts per cell in each cluster
df %>%
  ggplot(aes(x = nFeature_RNA, 
             y = nCount_RNA)) +
  geom_point() +
  lims(x = c(0, 1200),
       y = c(0, 2600))

# Number of genes in test versus the nCount_RNA per cell average
df %>%
  ggplot(aes(x = nCount_RNA, 
             y = genes)) +
  geom_point() +
  lims(x = c(0, 2600),
       y = c(0, 500))

# Number of genes in test versus the nFeature_RNA per cell average
df %>%
  ggplot(aes(x = nFeature_RNA,
             y = genes)) +
  geom_point() +
  lims(x = c(0, 2600),
       y = c(0, 500))
