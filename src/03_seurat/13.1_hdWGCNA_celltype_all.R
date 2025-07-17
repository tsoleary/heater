# ------------------------------------------------------------------------------
# hdWGCNA on cell type subsets of the data
# High dimensional weighted gene co-expression analysis (hdWGCNA) of RNA and 
#   ATAC data 
# TS O'Leary
# ------------------------------------------------------------------------------

# Load libraries
library(tidyverse)
library(Seurat)
library(Signac)
library(hdWGCNA)
library(cowplot)

# Load data
dat_all <- readRDS(here::here("data/processed/seurat_object/10_dat_linked.rds"))

cell_types <- unique(dat_all$cell_type)

# Create functions for easy repetition on each cell type -----------------------

# Function: setup_hdWGCNA -----
# Description: 
# Inputs: dat, and cell type
# Outputs: soft power plots

setup_hdWGCNA <- function(dat, ct) {
  
  # Define output directory and create clean name for cell type file
  out_dir <- "output/wgcna/cell_type/soft_power"
  ct_clean <- ct |> str_remove("\\.") |> str_replace_all(" ", "_")
  
  # Subset for cell type
  dat <- subset(dat_all, subset = cell_type == ct)
  
  # Default assay to SCT
  DefaultAssay(dat) <- "SCT"
  
  # Set up for WGCNA
  dat <- SetupForWGCNA(
    dat,
    gene_select = "fraction", 
    fraction = 0.05, 
    wgcna_name = "hdWGCNA" 
  )
  
  # Construct metacells in each group
  dat <- MetacellsByGroups(
    seurat_obj = dat,
    group.by = c("acc_temp", "cell_type"),
    reduction = "umap",
    k = 25, 
    max_shared = 10,
    ident.group = "acc_temp" 
  )
  
  # Normalize metacell expression matrix
  dat <- NormalizeMetacells(dat)
  
  # Set the data for the 
  dat <- SetDatExpr(
    dat,
    group_name = "18°C",
    group.by = "acc_temp", 
    assay = "SCT",
    layer = "data"
  )
  
  # Figure out what 
  dat <- TestSoftPowers(
    dat, 
    networkType = "signed"
  )
  
  # Create soft power plot
  p_sp <- patchwork::wrap_plots(PlotSoftPowers(dat), ncol = 2) +
    patchwork::plot_annotation(title = paste("hdWGCNA soft power for", ct))
  
  # Save the plot for soft power selection
  ggsave(plot = p_sp, 
         filename = here::here(out_dir, paste0(ct_clean, ".pdf")),
         height = 8,
         width = 12,
         units = "in")
  
  # Return dat for later processing
  return(dat)
  
} ## End function -----

# Function: run_hdWGCNA -----
# Description: A function to run the hdWGCNA on cell-type subsets
# Inputs: setup_WGCNA prepared data
# Outputs: DMEs, modules, and processed dat

# Required packages

# Example call

run_hdWGCNA <- function(dat_all, soft_power, ct) {
  
  # Define output directory and create clean name for cell type file
  out_dir <- "output/wgcna/cell_type"
  ct_clean <- ct |> str_remove("\\.") |> str_replace_all(" ", "_")
  
  dat <- setup_hdWGCNA(dat_all, ct = ct)
  
  # Construct network 
  dat <- ConstructNetwork(
    dat, 
    soft_power = soft_power,
    tom_outdir = here::here(out_dir, ct |> 
                              str_remove("\\.") |> 
                              str_replace_all(" ", "_")),
    tom_name = "hdWGCNA",
    overwrite_tom = TRUE
  )
  
  # Scale data using variable features
  dat <- ScaleData(dat, features = VariableFeatures(dat))
  
  # Compute module eigengenes
  dat <- ModuleEigengenes(dat, group.by.vars = "acc_temp")
  
  # Identify differential module eigengenes between acclimation treatments
  DMEs <- FindDMEs(
    dat,
    barcodes1 = dat@meta.data |>  filter(acc_temp == "18°C") |> rownames(),
    barcodes2 = dat@meta.data |> filter(acc_temp == "25°C") |> rownames(),
    test.use = "MAST",
    wgcna_name = "hdWGCNA") |> 
    as_tibble()
  
  # Save outputs for this cell type
  saveRDS(
    DMEs,
    here::here(out_dir, paste0(ct_clean, "_DMEs.rds"))
  )
  saveRDS(
    GetModules(dat) |> as_tibble(), 
    here::here(out_dir, paste0(ct_clean, "_modules.rds"))
  )
  saveRDS(
    dat, 
    here::here(
      "data/processed/seurat_object/wgcna", 
      paste0("dat_wgcna_", ct_clean, ".rds"))
  )
  
  # Plot dendrogram
  png(here::here(out_dir, "dendrogram", paste0(ct_clean, "_dendrogram.png")),
      width = 1000, height = 600, units = "px", res = 100)
  PlotDendrogram(
    dat, 
    main = paste(ct, "hdWGCNA Dendrogram")
  )
  dev.off()
  
  return(dat)
  
} ## End function -----

# Get soft power threshold for all cell types ----
for (ct in cell_types) {
  setup_hdWGCNA(dat_all, ct = ct)
}

# Run hdWGCNA with defined soft powers for all cell types -----
cell_type_sp <- tibble(
  "cell_type" = cell_types,
  "soft_power" = c(7, 7, 7, 10, 8, 10, 5, 10, 8)
)

# Loop through all cell types
for (i in 9:nrow(cell_type_sp)) {
  print(paste("################################################"))
  print(paste("############",  cell_type_sp$cell_type[i], "############"))
  print(paste("################################################"))
  dat <- run_hdWGCNA(
    dat_all, 
    soft_power = cell_type_sp$sp[i], 
    ct = cell_type_sp$cell_type[i]
  )
}


################################################################################

# Processing and plotting some of these results

# Load the DMEs
DMEs <- list.files(here::here("output/wgcna/cell_type")) |> 
  str_subset("DMEs") %>%
  here::here("output/wgcna/cell_type", .) %>%
  set_names(basename(.)) |> 
  purrr::map_df(readRDS, .id = "cell_type") |> 
  mutate(cell_type = str_remove(cell_type, "_DMEs.rds")) |> 
  arrange(cell_type, desc(avg_log2FC)) |> 
  mutate(mod_number = row_number()) |> 
  group_by(cell_type) |> 
  mutate(mod_number_ct = row_number())

modules <- list.files(here::here("output/wgcna/cell_type")) |> 
  str_subset("modules") %>%
  here::here("output/wgcna/cell_type", .) %>%
  set_names(basename(.)) |> 
  purrr::map_df(readRDS, .id = "cell_type") |> 
  mutate(cell_type = str_remove(cell_type, "_modules.rds"))

DMEs <- modules |> 
  group_by(cell_type, module) |> 
  tally() |> 
  full_join(DMEs)

# Save DMEs list
saveRDS(DMEs, here::here("output/wgcna/celltype_DMEs.rds"))

DMEs <- readRDS(here::here("output/wgcna/celltype_DMEs.rds"))
DMEs |> 
  mutate(mod_number = factor(mod_number, levels = 66:1)) |> 
  mutate(cell_type = str_replace_all(cell_type, "_", " ")) |>
  mutate(cell_type = str_replace(cell_type, "prim", "prim.")) |> 
  filter(module != "grey") |> 
  arrange(avg_log2FC) |> 
  ggplot(aes(y = mod_number,
             x = -avg_log2FC)) +
  geom_vline(xintercept = 0, color = "grey50") +
  geom_segment(aes(xend = 0, 
                   yend = mod_number,
                   color = p_val_adj < 0.05),
               linewidth = 0.5) +
  geom_point(aes(fill = p_val_adj < 0.05, 
               size = n), 
           color = "grey50",
           shape = 21) +
  scale_size_continuous(range = c(1, 5),
                        name = "Number of\ngenes",
                        breaks = c(50, 100, 300, 600, 900)) +
  scale_fill_manual(values = c("grey80", "firebrick")) +
  scale_color_manual(values = c("grey80", "firebrick")) +
  scale_x_continuous(limits = c(-10, 10)) +
  guides(fill = "none", color = "none") +
  labs(y = "Module number",
       x = expression("log"[2]*"(fold-difference)")) +
  theme_bw() +
  theme(panel.grid.major.y = element_line(color = "grey95"),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        legend.position = c(0.85, 0.125)) +
  facet_wrap(~cell_type, scales = "free_y")


ggsave(here::here("output/figs/wgcna/cell_type_mods.pdf"),
       height = 15,
       width = 20,
       units = "cm")


modules_numbers <- DMEs |> 
  select(mod_number, cell_type, module) |> 
  arrange(mod_number, cell_type) |> 
  ungroup(cell_type) |> 
  mutate(mod_number = row_number()) |> 
  full_join(modules, by = c("cell_type", "module")) |> 
  select(-color)

DMEs_sig <- DMEs |> 
  filter(p_val_adj < 0.05)


for (i in DMEs_sig$mod_number) {
  modules_numbers |> 
    filter(mod_number == i) |> 
    select(gene_name) |>
    write_tsv(
      here::here("output/wgcna/cell_type/modules", 
                 paste0("module_", i, ".txt")),
      col_names = FALSE)
}

mod_go <- list.files(here::here("output/wgcna/cell_type/modules/go_ora")) |> 
  str_subset(".txt") %>%
  here::here("output/wgcna/cell_type/modules/go_ora", .) %>%
  set_names(basename(.)) |> 
  purrr::map_df(read_tsv, .id = "module") |> 
  mutate(module = str_remove(module, ".txt"))

p_mod_1 <- mod_go |> 
  filter(module %in% c("mod_18", "mod_19", "mod_29", "mod_30", "mod_38")) |>
  mutate(module = str_remove(module, "mod_")) |> 
  mutate(description = paste0(description, " (", module, ")")) |> 
  mutate(description = fct_reorder(description, enrichmentRatio, .fun = mean)) |> 
  ggplot(aes(y = description,
             x = enrichmentRatio)) +
  geom_point(aes(size = size,
                 fill = FDR), 
             color = "grey50",
             shape = 21) +
  scale_y_discrete(name = element_blank(),
                   labels = function(x) str_remove(x, " \\(.*\\)")) +
  scale_size_continuous(breaks = c(25, 50, 100, 150, 200)) +
  scale_x_continuous(limits = c(0, NA)) +
  scale_fill_gradient(low = "firebrick", high = "grey80", limits = c(0, 0.1)) +
  coord_cartesian(clip = "off") +
  theme_minimal() + 
  theme(strip.text.y = element_text(vjust = 1, angle = 0, face = "bold", size = 12,
                                    margin = margin(1, 1, 1, 1)),
        strip.background =  element_rect(fill = "grey99", color = "grey80"),
        panel.border = element_rect(color = "grey80", fill = NA, size = 0.25),
        legend.position = "none") +
  facet_grid(module~., scales = "free_y", space = "free_y")


p_mod_2 <- mod_go |> 
  filter(module %in% c("mod_39", "mod_43", "mod_52")) |>
  mutate(module = str_remove(module, "mod_")) |> 
  mutate(description = paste0(description, " (", module, ")")) |> 
  mutate(description = fct_reorder(description, enrichmentRatio, .fun = mean)) |> 
  ggplot(aes(y = description,
             x = enrichmentRatio)) +
  geom_point(aes(size = size,
                 fill = FDR), 
             color = "grey50",
             shape = 21) +
  scale_y_discrete(name = element_blank(),
                   labels = function(x) str_remove(x, " \\(.*\\)")) +
  scale_size_continuous(breaks = c(25, 50, 100, 150, 200)) +
  scale_x_continuous(limits = c(0, NA)) +
  scale_fill_gradient(low = "firebrick", high = "grey80", limits = c(0, 0.1)) +
  coord_cartesian(clip = "off") +
  theme_minimal() + 
  theme(strip.text.y = element_text(vjust = 1, angle = 0, face = "bold", size = 12,
                                    margin = margin(1, 1, 1, 1)),
        strip.background =  element_rect(fill = "grey99", color = "grey80"),
        panel.border = element_rect(color = "grey80", fill = NA, size = 0.25),
        legend.position = "none") +
  facet_grid(module~., scales = "free_y", space = "free_y")


p_mod_3 <- mod_go |> 
  filter(module %in% c("mod_53", "mod_60")) |>
  mutate(module = str_remove(module, "mod_")) |> 
  mutate(description = paste0(description, " (", module, ")")) |> 
  mutate(description = fct_reorder(description, enrichmentRatio, .fun = mean)) |> 
  ggplot(aes(y = description,
             x = enrichmentRatio)) +
  geom_point(aes(size = size,
                 fill = FDR), 
             color = "grey50",
             shape = 21) +
  scale_y_discrete(name = element_blank(),
                   labels = function(x) str_remove(x, " \\(.*\\)")) +
  scale_size_continuous(breaks = c(25, 50, 100, 150, 200),
                        name = "Number of\ngenes") +
  scale_x_continuous(limits = c(0, NA)) +
  scale_fill_gradient(low = "firebrick", high = "grey80", limits = c(0, 0.1),
                      breaks = c(0, 0.025, 0.05, 0.075, 0.1),
                      labels = c(0, 0.025, 0.05, 0.075, 0.1)) +
  coord_cartesian(clip = "off") +
  theme_minimal() + 
  theme(strip.text.y = element_text(vjust = 1, angle = 0, face = "bold", size = 12,
                                    margin = margin(1, 1, 1, 1)),
        strip.background =  element_rect(fill = "grey99", color = "grey80"),
        panel.border = element_rect(color = "grey80", fill = NA, size = 0.25),
        legend.position = "bottom",
        legend.direction = "vertical") +
  facet_grid(module~., scales = "free_y", space = "free_y")


cowplot::plot_grid(
  p_mod_1,
  p_mod_2,
  p_mod_3, 
  nrow = 1
)

ggsave(here::here("output/figs/wgcna/go_ora_mods.pdf"),
       height = 45,
       width = 50,
       units = "cm")


DMEs |> 
  group_by(avg_log2FC > 0) |> 
  tally()


DMEs |> 
  filter(p_val_adj < 0.05) |> 
  group_by(avg_log2FC > 0) |> 
  tally()
