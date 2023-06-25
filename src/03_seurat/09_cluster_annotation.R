# ------------------------------------------------------------------------------
# Annotating Seurat Clusters using Fisher's enrichment of BDGP in situ data
# June 12, 2023
# TS O'Leary
# ------------------------------------------------------------------------------

# Load libraries
library(tidyverse)

# Load data --------------------------------------------------------------------

# Berkeley Drosophila Genome Project in situ database
# https://insitu.fruitfly.org/insitu-mysql-dump/insitu_annot.csv.gz
insitu_annot <- read_csv(here::here("data/raw/annot/insitu_annot.csv"))

# Cluster annotations from Calderon data at 4 to 6 hours 
cluster_annot_calderon <- read_csv(
  here::here("data/raw/annot/calderon_cluster_annot.csv")) |>
  filter(`time window (modeled)` == "4-6") |>
  dplyr::rename(cluster = `cluster #`) |>
  dplyr::select(cluster, `selected annotation (manual)`) |> 
  distinct(`selected annotation (manual)`)

# Marker genes from the Seurat clusters -- all markers with pval < 0.05
markers <- readRDS(here::here("data/processed/genes/markers.rds")) |> 
  group_by(cluster) |> 
  filter(max_pval < 0.05)

# Fisher's Exact Test for enrichment of cell-type specific marker genes --------

# Reference panel of a specific level of the in situ data
ref_panel <- insitu_annot
cluster_annot <- NULL

pb <- progress::progress_bar$new(
  total = length(unique(markers$cluster)) * length(unique(ref_panel$annot)))
# Loop through each cluster present in the query markers data
for (cluster_i in unique(markers$cluster)) {
  # Unique markers for a cluster
  query_markers <- markers |> 
    filter(cluster == cluster_i) |> 
    ungroup() |> 
    distinct(gene) |> 
    deframe()
  
  # Loop through each annotation type in the reference
  for (ref_annot_j in unique(ref_panel$annot)) {
    
    # Distinct reference markers for this annotation
    ref_markers <- ref_panel |> 
      filter(annot == ref_annot_j) |> 
      distinct(gene) |> 
      deframe()
    
    # Collect numbers for the Fisher's exact test contingency table ------------
    
    # Number of markers in common
    n_common <- length(intersect(query_markers, ref_markers))
    # Number unique to the query
    n_query_unique <- length(setdiff(query_markers, ref_markers))
    # Number unique to the reference
    n_ref_unique <- length(setdiff(ref_markers, query_markers))		
    # Number remaining
    n_remain <- length(unique(c(unique(markers$gene), rownames(ref_panel)))) - 
      n_common - n_query_unique - n_ref_unique
    
    # Proportion of query markers that are shared with the reference
    prop_overlap <- n_common/length(query_markers)
    
    # Contingency table
    con_table <- matrix(
      c(n_common, 
        n_query_unique, 
        n_ref_unique, 
        n_remain), 
      ncol = 2)
    
    # Run the Fisher's exact test and save the p-value
    fet_pval <- fisher.test(con_table)$p.value
    
    # Store results
    df <- tibble(
      cluster = cluster_i,
      annot = ref_annot_j,
      prop_overlap = prop_overlap,
      pval = fet_pval
    )
    
    # Bind results together into a dataframe
    cluster_annot <- bind_rows(df, cluster_annot)
    
    # Print progress bar
    pb$tick()
  }
}

pb$terminate()

# Adjusted pval for FDR correction
cluster_annot <- cluster_annot |> 
  mutate(padj = p.adjust(pval, method = "BH"))

x <- cluster_annot |> 
  group_by(cluster) |> 
  filter(padj < 0.05) |> 
  tally()

# Save the full annotation results
saveRDS(cluster_annot, here::here("data/processed/annot/cluster_annot_all.rds"))

# Add in the consensus results
cluster_annot_top <- cluster_annot |> 
  filter(padj < 0.05) |> 
  mutate(cluster = as.numeric(cluster)) |> 
  group_by(cluster) |> 
  arrange(padj)
  slice_min(padj, n = 3)
  
# Four doesn't have an annotation that passes FDR 
# so just looking at top annotations here
cluster_annot |> 
  filter(cluster == 4) |> 
  arrange(padj) |> 
  slice_min(padj, n = 10)

cluster_annot |> 
  filter(cluster == 13) |> 
  arrange(padj) |> 
  slice_min(padj, n = 10)

cluster_annot |> 
  filter(cluster == 27) |> 
  arrange(padj) |> 
  slice_min(padj, n = 10)

# Add in the consensus results
cluster_annot_together <- cluster_annot_top |> 
  summarise(top_3 = paste0(annot, collapse = "; "))

# Manual annotations by TSO using the 11 annotations present in at 4 to 6 hours
# of development in the Calderon data
cluster_annot_manual <- list(
  "ubiquitous" = c(0, 2, 6, 16, 23, 27),
  "ectoderm primordium" = c(10, 11, 15, 18, 21),
  "endoderm primordium" = c(12, 13),
  "mesoderm primordium" = c(1, 5, 13, 14, 20),
  "ventral nerve cord primordium" = c(9),
  "peripheral nervous system primordium" = c(4, 8),
  "tracheal primordium" = c(7),
  "foregut/hindgut primordium" = c(3, 19),
  "amnioserosa" = c(22),
  "yolk" = c(17, 26),
  "unknown" = c(24, 25)
)

# Save the consensus results
saveRDS(
  cluster_annot_manual, 
  here::here("data/processed/annot/cluster_annot_manual.rds")
)


# Recode metadata with cluster annotations -----

# Load data
dat <- readRDS(here::here("data/processed/seurat_object/07_dat_cluster.rds"))

# Recode
dat@meta.data <- dat@meta.data |> 
  mutate(
    cell_type = 
      case_match(
        seurat_clusters,
        c("0", "2", "6", "16", "23", "27") ~ "ubiquitous",
        c("10", "11", "15", "18", "21") ~ "ectoderm prim.",
        c("12", "13") ~ "endoderm prim.",
        c("1", "5", "13", "14", "20") ~ "mesoderm prim.",
        c("9") ~ "ventral nerve cord prim.",
        c("4", "8") ~ "peripheral nervous system prim.",
        c("7") ~ "tracheal prim.",
        c("3", "19") ~ "foregut/hindgut prim.",
        c("22") ~ "amnioserosa",
        c("17", "26") ~ "yolk nuclei",
        c("24", "25") ~ "unknown")
    )

# Save data
saveRDS(dat, here::here("data/processed/seurat_object/09_dat_annot.rds"))
