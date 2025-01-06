# ------------------------------------------------------------------------------
# Annotating Seurat Clusters using Fisher's enrichment of BDGP in situ data
# TS O'Leary
# ------------------------------------------------------------------------------

# Load libraries
library(tidyverse)
library(Seurat)
library(Signac)

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
markers <- readRDS(here::here("output/markers/cluster_markers.rds")) |> 
  group_by(cluster) |> 
  filter(max_pval < 0.05)

# Output dir 
out_dir <- "data/processed/seurat_object"

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

cluster_annot |> 
  mutate(cluster = as.numeric(cluster)) |> 
  group_by(cluster) |> 
  filter(padj < 0.05) |> 
  tally()

# Save the full annotation results
saveRDS(cluster_annot, here::here("data/processed/annot/cluster_annot_all.rds"))

# Top results
cluster_annot_top <- cluster_annot |> 
  mutate(cluster = as.numeric(cluster)) |> 
  filter(padj < 0.05 | prop_overlap > 0.01) |> 
  group_by(cluster) |> 
  slice_min(padj, n = 5) |> 
  arrange(cluster, padj) 

# Recode metadata with manual cluster annotations ------------------------------

# Load data
dat <- readRDS(here::here("data/processed/seurat_object/07_dat_cluster.rds"))
dat$seurat_clusters <- dat@meta.data$wknn_res.0.7

# Recode
dat@meta.data <- dat@meta.data |> 
  mutate(
    cell_type = 
      case_match(
        seurat_clusters,
        c("2") ~ "ectoderm prim.",
        c("3") ~ "endoderm prim.",
        c("0", "6") ~ "mesoderm prim.",
        c("1") ~ "ventral nerve cord prim.",
        c("4", "11") ~ "peripheral nervous system prim.",
        c("7") ~ "tracheal prim.",
        c("8", "9") ~ "foregut & hindgut prim.",
        c("10") ~ "amnioserosa",
        c("5") ~ "germ cell")
    )

# Save data with cell-type annotations added
saveRDS(dat, here::here(out_dir, "09_dat_annot.rds"))

# Save summary of the manual annotations ---------------------------------------
annot <- dat@meta.data |> 
  rownames_to_column("gene") |> 
  select(seurat_clusters, cell_type) |> 
  distinct(seurat_clusters, .keep_all = TRUE) |> 
  arrange(seurat_clusters)

saveRDS(annot, here::here("data/processed/annot", "annot.rds"))

# Analyze the number of cells in each cluster and annotation per replicate -----

# Plot the percentage of cells in each cell type
dat@meta.data |> 
  group_by(cell_type, acc_temp) |> 
  tally() |> 
  ungroup(cell_type) |> 
  group_by(acc_temp) |> 
  mutate(percent = (n/sum(n))) |> 
  mutate(cell_type = fct_reorder(cell_type, percent)) |> 
  ggplot() +
  geom_col(aes(y = cell_type, 
               x = percent,
               fill = acc_temp),
           color = "grey80",
           position = "dodge") +
  scale_fill_manual(values = c("#8698C0", "#D09F7B"),
                    name = element_blank()) +
  cowplot::theme_minimal_vgrid() +
  scale_x_continuous(expand = c(0, 0), 
                     labels = scales::percent,
                     name = "Cell-type composition",
                     position = "top") +
  theme(legend.position = "bottom",
        axis.title.y = element_blank())


# Plot the percentage of well
dat@meta.data |> 
  group_by(cell_type, sample_name) |> 
  tally() |> 
  ungroup(cell_type) |> 
  group_by(sample_name) |> 
  mutate(percent = (n/sum(n))) |> 
  mutate(cell_type = fct_reorder(cell_type, percent)) |> 
  ggplot() +
  geom_col(aes(y = cell_type, 
               x = percent,
               fill = sample_name),
           color = "grey80",
           position = "dodge") +
  scale_fill_manual(values = c("#668CDE","#8698C0","#D08956", "#D09F7B"),
                    name = element_blank()) +
  cowplot::theme_minimal_vgrid() +
  scale_x_continuous(expand = c(0, 0), 
                     labels = scales::percent,
                     name = "Cell-type composition",
                     position = "top") +
  theme(legend.position = "bottom",
        axis.title.y = element_blank())

# Create contingency table
t <- dat@meta.data |> 
  group_by(cell_type, sample_name) |> 
  tally() |> 
  pivot_wider(names_from = sample_name,
              values_from = n) |> 
  column_to_rownames("cell_type")

# Run a chi-square test
chisq <- chisq.test(t)
print(chisq)


# Create contingency table
t <- dat@meta.data |> 
  group_by(cell_type, acc_temp) |> 
  tally() |> 
  pivot_wider(names_from = acc_temp,
              values_from = n) |> 
  column_to_rownames("cell_type")

# Run a chi-square test
chisq <- chisq.test(t)
print(chisq)

# Plot out the residuals ----
chisq$residuals |> 
  as.data.frame() |> 
  rownames_to_column("sample_name") |> 
  pivot_longer(cols = amnioserosa:`ventral nerve cord prim.`,
               values_to = "residuals",
               names_to = "cell_type") |> 
  ggplot() +
  geom_point(
    aes(x = sample_name,
        y = cell_type,
        fill = residuals,
        size = abs(residuals)),
    color = "grey70",
    shape = 21) +
  scale_fill_gradient2(low = scales::muted("blue"), 
                       high = scales::muted("red")) +
  cowplot::theme_minimal_grid() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  guides(size = "none")

# Alternative method directly accounting for replication
t <- dat@meta.data |> 
  group_by(cell_type, sample_name, acc_temp) |> 
  tally() |> 
  mutate(replicate = str_extract(sample_name, "\\d$"))

model <- lme4::glmer(
  n ~ acc_temp * cell_type + (1 | cell_type/replicate), 
  data = t, 
  family = poisson)

summary(model)


# Chi-square individually 

