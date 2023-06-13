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
cluster_annot <- cluster_annot %>%
  mutate(padj = p.adjust(pval, method = "BH"))

# Save the full annotation results
saveRDS(cluster_annot, here::here("data/processed/annot/cluster_annot_all.rds"))

# Manual annotations by TSO using the 11 annotations present in at 4 to 6 hours
# of development in the Calderon data
cluster_annot_manual <- c(
  "0" = "ectoderm primordium",
  "1" = "mesoderm primordium",
  "2" = "ventral nerve cord primordium",
  "3" = "endoderm primordium",
  "4" = "ectoderm primordium",
  "5" = "muscle system primordium",
  "6" = "tracheal primordium",
  "7" = "ventral nerve cord primordium",
  "8" = "plasmatocytes anlage",
  "9" = "amnioserosa",
  "10" = "unknown"
)

# Add in the consensus results
cluster_annot_manual <- cluster_annot |> 
  filter(padj < 0.05) |> 
  mutate(cluster = as.numeric(cluster)) |> 
  group_by(cluster) |> 
  slice_min(padj, n = 3) |> 
  summarise(top_3 = paste0(annot, collapse = "; ")) |> 
  mutate(annot_manual = cluster_annot_manual) |> 
  select(cluster, annot_manual, top_3)

# Save the consensus results
saveRDS(
  cluster_annot_manual, 
  here::here("data/processed/annot/cluster_annot_manual.rds")
)
