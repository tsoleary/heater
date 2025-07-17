# ------------------------------------------------------------------------------
# High dimensional weighted gene co-expression analysis (hdWGCNA) of RNA and 
#   ATAC data 
# TS O'Leary
# ------------------------------------------------------------------------------

# Load libraries
library(hdWGCNA)
library(tidyverse)
library(Seurat)
library(Signac)
library(cowplot)
library(patchwork)

# Output dir
out_dir <- "output/wgcna"

# Load data
dat <- readRDS(here::here("data/processed/seurat_object/10_dat_linked.rds"))
DefaultAssay(dat) <- "SCT"

# Set up the data -- exclude genes expressed in less than 5% of cells
dat <- SetupForWGCNA(
  dat,
  gene_select = "fraction",
  fraction = 0.05,
  wgcna_name = "wgcna"
)

# Construct metacells  in each group
dat <- MetacellsByGroups(
  seurat_obj = dat,
  group.by = c("cell_type", "acc_temp"), 
  reduction = "umap",
  k = 25,
  max_shared = 10,
  ident.group = "cell_type"
)

# Normalize meta cell expression matrix
dat <- NormalizeMetacells(
  dat
)

# Set dat expr
dat <- SetDatExpr(
  dat,
  group_name = "18°C",
  assay = "SCT",
  group.by = "acc_temp",
  slot = "data"
)

# Test different soft powers
dat <- TestSoftPowers(
  dat,
  networkType = "signed"
)

# Construct co-expression network:
dat <- ConstructNetwork(
  dat, 
  soft_power = 8,
  setDatExpr = FALSE,
  tom_outdir = out_dir,
  tom_name = "wgcna",
  overwrite_tom = TRUE
)

# Scale data first
dat <- ScaleData(
  dat,
  features = VariableFeatures(dat)
)

# Compute all MEs in the full single-cell dataset
dat <- ModuleEigengenes(
  dat,
  group.by.vars = "acc_temp"
)

# Find differential module eigengenes
DMEs <- FindDMEs(
  dat,
  barcodes1 = dat@meta.data |> filter(acc_temp == "18°C") |> rownames(),
  barcodes2 = dat@meta.data |> filter(acc_temp == "25°C") |> rownames(),
  test.use = "wilcox",
  wgcna_name = "wgcna"
)

# Save DMEs
saveRDS(DMEs, here::here(out_dir, "DMEs.rds"))

# Module connectivity
# compute eigengene-based connectivity (kME):
dat <- ModuleConnectivity(
  dat,
  group.by = "acc_temp", 
  group_name = "25°C"
)

PlotKMEs(
  dat,
  ncol = 3,
  plot_widths = c(3, 1),
  text_size = 3,
  n_hubs = 8
)

# Save
ggsave(here::here(fig_dir, "kMEs.png"),
       height = 20,
       width = 25,
       units = "cm")

# Gene module assignment
modules <- GetModules(dat)

# Get genes from significant modules
brown <- modules |> 
  filter(module == "brown") |> 
  pull(gene_name)

# Get genes from significant modules
red <- modules |> 
  filter(module == "red") |> 
  pull(gene_name)

# Get genes from green
green <- modules |> 
  filter(module == "green") |> 
  pull(gene_name)
blue <- modules |> 
  filter(module == "blue") |> 
  pull(gene_name)
turquoise <- modules |> 
  filter(module == "turquoise") |> 
  pull(gene_name)
black <- modules |> 
  filter(module == "black") |> 
  pull(gene_name)
yellow <- modules |> 
  filter(module == "yellow") |> 
  pull(gene_name)


# Save modules
saveRDS(modules, here::here(out_dir, "modules.rds"))
write_csv(brown |> as_tibble(), here::here(out_dir, "brown.csv"))
write_csv(red |> as_tibble(), here::here(out_dir, "red.csv"))
write_csv(green |> as_tibble(), here::here(out_dir, "green.csv"))
write_csv(blue |> as_tibble(), here::here(out_dir, "blue.csv"))
write_csv(turquoise |> as_tibble(), here::here(out_dir, "turquoise.csv"))
write_csv(black |> as_tibble(), here::here(out_dir, "black.csv"))
write_csv(yellow |> as_tibble(), here::here(out_dir, "yellow.csv"))


# Save dat with hdWGCNA in it for plotting
saveRDS(dat, here::here("data/processed/seurat_object/13_dat_wgcna.rds"))

# # ATAC -------------------------------------------------------------------------
# 
# # Output dir
# out_dir <- "output/wgcna/atac"
# 
# # Load data
# dat <- readRDS(here::here("data/processed/seurat_object/10_dat_linked.rds"))
# DefaultAssay(dat) <- "peaks"
# 
# # Set up the data -- exclude genes expressed in less than 5% of cells
# dat <- SetupForWGCNA(
#   dat,
#   gene_select = "fraction",
#   fraction = 0.05,
#   wgcna_name = "wgcna"
# )
# 
# # Construct metacells  in each group
# dat <- MetacellsByGroups(
#   seurat_obj = dat,
#   group.by = c("cell_type", "acc_temp"), 
#   reduction = "umap",
#   k = 25,
#   max_shared = 10,
#   ident.group = "cell_type"
# )
# 
# # Normalize meta cell expression matrix
# dat <- NormalizeMetacells(
#   dat
# )
# 
# # Set dat expr
# dat <- SetDatExpr(
#   dat,
#   group_name = "18°C",
#   assay = "peaks",
#   group.by = "acc_temp",
#   slot = "data"
# )
# 
# # Test different soft powers
# dat <- TestSoftPowers(
#   dat,
#   networkType = "signed"
# )
# 
# # Construct co-expression network
# dat <- ConstructNetwork(
#   dat, 
#   soft_power = 5,
#   setDatExpr = FALSE,
#   tom_outdir = out_dir,
#   tom_name = "wgcna_atac",
#   overwrite_tom = TRUE
# )
# 
# # Scale data first
# dat <- ScaleData(
#   dat,
#   features = VariableFeatures(dat)
# )
# 
# # Compute all MEs in the full single-cell dataset
# dat <- ModuleEigengenes(
#   dat,
#   assay = "peaks",
#   group.by.vars = "acc_temp"
# )
# 
# # Find differential module eigengenes
# DMEs <- FindDMEs(
#   dat,
#   barcodes1 = dat@meta.data |> filter(acc_temp == "18°C") |> rownames(),
#   barcodes2 = dat@meta.data |> filter(acc_temp == "25°C") |> rownames(),
#   test.use = "wilcox",
#   wgcna_name = "wgcna"
# )
# 
# # Save DMEs
# saveRDS(DMEs, here::here(out_dir, "DMEs.rds"))
# 
# # Module connectivity
# # compute eigengene-based connectivity (kME):
# dat <- ModuleConnectivity(
#   dat,
#   group.by = "acc_temp", 
#   group_name = "25°C"
# )
# 
# PlotKMEs(
#   dat,
#   ncol = 3,
#   plot_widths = c(3, 1),
#   text_size = 3,
#   n_hubs = 8
# )
# 
# # Gene module assignment
# modules <- GetModules(dat)
# 
# # Get regions from significant modules
# brown <- modules |> 
#   filter(module == "brown") |> 
#   pull(gene_name)
# 
# # Get regions from significant modules
# blue <- modules |> 
#   filter(module == "blue") |> 
#   pull(gene_name)
# 
# # Save modules
# saveRDS(modules, here::here(out_dir, "modules.rds"))
# write_csv(brown |> as_tibble(), here::here(out_dir, "brown.csv"))
# write_csv(blue |> as_tibble(), here::here(out_dir, "blue.csv"))
# 
# # Get genes linked to the blue peaks
# blue_genes <- GetLinkedGenes(
#   dat, 
#   features = blue,
#   min.abs.score = 0.1
# )
# 
# # Get genes linked to the brown peaks
# brown_genes <- GetLinkedGenes(
#   dat, 
#   features = brown,
#   min.abs.score = 0.1
# )
# 
# # Save the list of linked peaks
# write_csv(brown_genes |> as_tibble(), here::here(out_dir, "brown_genes_0.1.csv"))
# write_csv(blue_genes |> as_tibble(), here::here(out_dir, "blue_genes_0.1.csv"))
# 
# # Save dat with hdWGCNA in it for plotting
# saveRDS(dat, here::here("data/processed/seurat_object/13_dat_wgcna_atac.rds"))
# 
# # Threhold for DME testing ----
# 
# MEs <- GetMEs(dat) |> 
#   select(-grey) |> 
#   rownames_to_column("cell") |> 
#   as_tibble() |> 
#   pivot_longer(-cell, 
#                names_to = "module",
#                values_to = "mod_eigengene") |> 
#   group_by(module) |> 
#   mutate(threshold = quantile(mod_eigengene, probs = 0)) |> 
#   left_join(as.data.frame(dat@meta.data) |> 
#               rownames_to_column("cell") |> 
#               dplyr::select(cell, acc_temp))
# 
# 
# MEs |> 
#   ggplot() +
#   geom_histogram(aes(x = mod_eigengene),
#                  color = "grey20",
#                  fill = "grey80",
#                  bins = 50) +
#   geom_vline(aes(xintercept = threshold), 
#              color = "red",
#              linetype = 2) +
#   facet_wrap(~module,
#              scales = "free_x",
#              nrow = 2) +
#   theme_minimal() +
#   theme(strip.background = element_rect(color = "grey20", fill = "grey95"))
# 
# MEs |> 
#   ggplot() +
#   geom_histogram(aes(x = mod_eigengene,
#                      fill = acc_temp),
#                  position = "dodge",
#                  color = NA,
#                  bins = 50) +
#   geom_vline(aes(xintercept = threshold), 
#              color = "red",
#              linetype = 2) +
#   facet_wrap(~module,
#              scales = "free_x",
#              nrow = 2) +
#   theme_minimal() +
#   theme(strip.background = element_rect(color = "grey20", fill = "grey95"))
# 
# MEs |> 
#   ggplot() +
#   geom_histogram(aes(x = mod_eigengene),
#                  color = "grey20",
#                  fill = "grey80",
#                  bins = 50) +
#   geom_vline(aes(xintercept = threshold), 
#              color = "red",
#              linetype = 2) +
#   facet_wrap(~module,
#              scales = "free_x",
#              nrow = 2) +
#   theme_minimal() +
#   theme(strip.background = element_rect(color = "grey20", fill = "grey95"))
# 
# 
# # Filter to only cells that make it to a certain threshold
# MEs_filered <- MEs |> 
#   filter(mod_eigengene > threshold)  |> 
#   group_split() 
#   
# # Create empty dataframe
# DMEs <- tibble()
# 
# for (i in 1:length(MEs_filered)) {
#   
#   # Test for DMEs with only top 50% of cells in each module
#   DME <- FindDMEs(
#     dat,
#     barcodes1 = MEs_filered[[i]] |> filter(acc_temp == "18°C") |> pull(cell),
#     barcodes2 = MEs_filered[[i]] |> filter(acc_temp == "25°C") |> pull(cell),
#     #test.use = "MAST",
#     wgcna_name = "wgcna") |> 
#     filter(module == unique(MEs_filered[[i]]$module)) |> 
#     as_tibble()
#   
#   # Append results to the specific module
#   DMEs <- DMEs |> 
#     bind_rows(DME)
# }
# 
# DMEs |> arrange(avg_log2FC)

# Something about which modules are in which cell-types
MEs |> 
  full_join(dat@meta.data |> rownames_to_column("cell") |> select(cell, cell_type)) |> 
  group_by(module, cell_type) |> 
  summarize(avg_mod = mean(mod_eigengene)) |> 
  ggplot(aes(x = module, 
             y = avg_mod, 
             label = cell_type)) +
  geom_point() +
  ggrepel::geom_label_repel() + 
  scale_color_manual() +
  theme_minimal_grid()


MEs |> 
  full_join(dat@meta.data |> rownames_to_column("cell") |> select(cell, cell_type)) |> 
  group_by(module, cell_type) |> 
  summarize(avg_mod = mean(mod_eigengene)) |> 
  mutate(module = factor(module, 
                         levels = c("red", "brown", "blue", 
                                    "black", "green", "yellow", "turquoise"))) |> 
  full_join(color_cell_type) |> 
  ggplot(aes(y = forcats::fct_rev(cell_type), 
             x = avg_mod,
             alpha = avg_mod > 0,
             fill = colors)) +
  geom_col(color = "grey50") +
  geom_vline(xintercept = 0, 
             color = "grey20") +
  scale_x_continuous(limits = c(-9.1, 9.1),
                     breaks = c(-8, -4, 0 , 4, 8),
                     name = "mean module eigengene score per cell type") +
  scale_fill_identity() +
  facet_wrap(~ module, nrow = 2) +
  scale_alpha_discrete(range = c(0.3, 1)) +
  cowplot::theme_minimal_grid() +
  theme(legend.position = "none",
        axis.title.y = element_blank())

MEs |> 
  full_join(dat@meta.data |> rownames_to_column("cell") |> select(cell, cell_type)) |> 
  group_by(module, cell_type) |> 
  summarize(avg_mod = mean(mod_eigengene)) |> 
  ungroup(cell_type) |> 
  summarize(sd = sd(avg_mod)) |> 
  arrange(sd)

MEs |> 
  full_join(dat@meta.data |> rownames_to_column("cell") |> select(cell, cell_type)) |> 
  group_by(module, cell_type) |> 
  summarize(avg_mod = mean(mod_eigengene)) |> 
  ungroup(cell_type) |> 
  summarize(sd = sd(avg_mod)) |> 
  arrange(sd) |> 
  mutate(module = fct_reorder(module, sd, .desc = TRUE)) |> 
  ggplot(aes(y = module, 
             x = sd)) +
  geom_col(color = "grey20",
           fill = "grey80", 
           width = 0.8) +
  geom_label(aes(label = sprintf("%.2f",round(sd, 2))),
             nudge_x = -0.25) +
  scale_x_continuous(expand = c(0,0.01)) +
  cowplot::theme_minimal_vgrid()
  


dMEs |> 
  full_join(dat@meta.data |> rownames_to_column("cell") |> select(cell, cell_type)) |> 
  group_by(module, cell_type) |> 
  summarize(mean_mod = mean(mod_eigengene)) |> 
  ungroup(cell_type) |> 
  summarize(sd = sd(mean_mod)) |> 
  arrange(sd)

MEs |> 
  full_join(dat@meta.data |> rownames_to_column("cell") |> select(cell, cell_type)) |> 
  group_by(module, cell_type) |> 
  summarize(mean_mod = mean(mod_eigengene)) |> 
  filter(module == "green") |> arrange(mean_mod)
