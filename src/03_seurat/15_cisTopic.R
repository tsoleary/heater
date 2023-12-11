# ------------------------------------------------------------------------------
# cisTopic analysis for SCENIC+ input
# TS O"Leary
# ------------------------------------------------------------------------------

# Load libraries
library(tidyverse)
library(Seurat)
library(Signac)
library(cisTopic)
library(densityClust)
library(Rtsne)
library(AUCell)
library(org.Dm.eg.db)
library(TxDb.Dmelanogaster.UCSC.dm6.ensGene)

# Output figure directory
fig_dir <- "output/figs/cis_topic"

# Load in processed data and set up cisTopicObject -----------------------------
dat <- readRDS(here::here("data/processed/seurat_object/14_dat_tf_motif.rds"))

# Create peaks_matrix object
peaks_matrix <- dat@assays$peaks@counts

# Change region naming to chr:111-112 format
rownames(peaks_matrix) <- str_replace(
  rownames(peaks_matrix), 
  "mitochondrion-genome", 
  "mt_g") |> 
  str_replace("-", ":") |> 
  str_replace("mt_g", "mitochondrion-genome")

# Create cisTopic object from peaks_matrix
cisTopicObject <- createcisTopicObject(
  count.matrix = peaks_matrix, 
  project.name = "heater"
)

# Add in all meta.data from previous Seurat analysis
cisTopicObject <- cisTopicObject |> 
  addCellMetadata(cell.data = dat@meta.data)

# Build the models -------------------------------------------------------------

# Run CGS model – this takes a non-trivial amount of time to run.
cisTopicObject <- runCGSModels(
  cisTopicObject,
  topic = c(2, 5, 10, 15, 20, 25, 30, 35, 40), 
  seed = 987, 
  nCores = 9, 
  burnin = 120, 
  iterations = 150, 
  addModels = FALSE
)

# Visualize the likelihood stabilization
logLikelihoodByIter(
  cisTopicObject, 
  select = c(2, 5, 10, 15, 20, 25, 30, 35, 40)
)


# Select the model
cisTopicObject <- selectModel(
  cisTopicObject, 
  select = 35,
  type = "maximum"
)


# Interpreting the models ------------------------------------------------------
cisTopicObject <- runtSNE(
  cisTopicObject, 
  target = "cell", 
  seed = 123, 
  pca = FALSE, 
  method = "Probability"
)

# Get topic-cell assignments
cell_assign <- modelMatSelection(cisTopicObject, "cell", "Probability")

# Run t-SNE and assign clusters based on topic-cell assignments ----------------
set.seed(123)
DR <- Rtsne::Rtsne(t(cell_assign), pca = FALSE)
# Compute distance matrix
DRdist <- dist(DR$Y)
# Density cluster
dclust <- densityClust(DRdist, gaussian = TRUE)
# Find clusters
dclust <- findClusters(dclust, rho = 50, delta = 2.5)

# Plot to check threshold for rho and delta
png(here::here(fig_dir, "dClust_thresholds_rho_delta.png"),
    res = 300,
    height = 20,
    width = 25,
    units = "cm")
options(repr.plot.width = 6, repr.plot.height = 6)
plot(dclust$rho,dclust$delta,
     pch = 20,
     cex = 0.6,
     xlab = "rho", 
     ylab = "delta")
points(
  dclust$rho[dclust$peaks],
  dclust$delta[dclust$peaks],
  col = "red",
  pch = 20,
  cex = 0.8
)
text(
  dclust$rho[dclust$peaks] - 2,
  dclust$delta[dclust$peaks] + 1.5,
  labels = dclust$clusters[dclust$peaks]
)
abline(v = 50)
abline(h = 2.5)
dev.off()

# Add cluster information to meta data
densityClust <- as.data.frame(list(as.factor(dclust$clusters)),
                              row.names = cisTopicObject@cell.names, 
                              col.names = c("densityClust"))
cisTopicObject <- addCellMetadata(cisTopicObject, densityClust)

# Compare densityClust and cell-types
png(here::here(fig_dir, "densityClust_cell-type.png"),
    res = 300,
    height = 20,
    width = 25,
    units = "cm")
par(mfrow = c(1, 2))
plotFeatures(
  cisTopicObject, 
  method = "tSNE", 
  target = "cell", 
  topic_contr = NULL, 
  colorBy = c("densityClust", 
              "cell_type"), 
  cex.legend = 0.8, 
  factor.max = 0.75, 
  dim = 2, 
  legend = TRUE
)
dev.off()

# Plot cell-topic heatmap
png(here::here(fig_dir, "cell_topic_heatmap.png"),
    res = 300,
    height = 25,
    width = 25,
    units = "cm")
cellTopicHeatmap(
  cisTopicObject,
  method = "Probability", 
  colorBy = c("densityClust")
)
dev.off()

# Plot cell-topic heatmap by cell-type
png(here::here(fig_dir, "cell_topic_heatmap_cell-type.png"),
    res = 300,
    height = 25,
    width = 25,
    units = "cm")
cellTopicHeatmap(
  cisTopicObject,
  method = "Probability", 
  colorBy = c("cell_type")
)
dev.off()

# Color tSNE by topic score for each topic
png(here::here(fig_dir, "tSNE_topic-score.png"),
    res = 300,
    height = 25,
    width = 25,
    units = "cm")
par(mfrow = c(5, 7))
plotFeatures(
  cisTopicObject, 
  method = "tSNE", 
  target = "cell", 
  topic_contr = "Probability", 
  colorBy = NULL, 
  cex.legend = 0.8, 
  factor.max = 0.75, 
  dim = 2, 
  legend = TRUE
)
dev.off()


################################################################################
################################################################################
################################################################################
# # Create predictive distribution matrix to estimate dropouts and build 
# # cell-specific region rankings that can be used with AUCell
# pred_matrix <- predictiveDistribution(cisTopicObject)
# # Obtain signatures
# # Need to figure out exactly what this means but it seems like it has to do
# # having a predifined set of peaks from something like chip-seq for a specific
# # Transcription Factor....
# path_to_signatures <- here::here("data/processed/seq/all/raw_feature_bc_matrix")
# Bulk_ATAC_signatures <- paste(path_to_signatures, list.files(path_to_signatures), sep="")
# labels  <- gsub("._peaks.narrowPeak", "", list.files(path_to_signatures))
# cisTopicObject <- getSignaturesRegions(cisTopicObject, Bulk_ATAC_signatures, labels=labels, minOverlap = 0.4)
# 
# # To only keep unique peaks per signature
# cisTopicObject@signatures <- llply(1:length(cisTopicObject@signatures), function (i) cisTopicObject@signatures[[i]][-which(cisTopicObject@signatures[[i]] %in% unlist(as.vector(cisTopicObject@signatures[-i])))]) 
# names(cisTopicObject@signatures) <- labels
# 
# # Compute cell rankings (Reference time: 9 min)
# library(AUCell)
# aucellRankings <- AUCell_buildRankings(pred.matrix, plot=FALSE, verbose=FALSE)
# 
# # Check signature enrichment in cells (Reference time: 1 min)
# cisTopicObject <- signatureCellEnrichment(cisTopicObject, aucellRankings, selected.signatures="all", aucMaxRank = 0.3*nrow(aucellRankings), plot=FALSE)
# 
# # Plot
# par(mfrow=c(2,2))
# plotFeatures(cisTopicObject, method="tSNE", target="cell", topic_contr=NULL, colorBy=c("CD4Tcell", "Mono", "Bcell", "NKcell"), cex.legend = 0.4, factor.max=.75, dim=2, legend=TRUE, intervals=10)
# 
# # Analysis of the regulatory topics --------------------------------------------
cisTopicObject <- getRegionsScores(
  cisTopicObject,
  method = "NormTop",
  scale = TRUE
)
# 
# 
# # Bigwig files of the cisTopics
# getBigwigFiles(
#   cisTopicObject, 
#   path = "output/cisTopics_asBW", 
#   seqlengths = seqlengths(txdb)
# )
# 
# # Visualize the 
# par(mfrow = c(2, 5))
cisTopicObject <- binarizecisTopics(
  cisTopicObject,
  thrP = 0.975,
  plot = TRUE
)
# 
# # Get bed files of the cisTopics
# getBedFiles(
#   cisTopicObject, 
#   path = "output/cisTopics_asBed"
# )
# 
# 
# cisTopicObject <- runtSNE(
#   cisTopicObject, 
#   target = "region", 
#   perplexity = 200, 
#   check_duplicates = FALSE
# )
# 
# 
# plotFeatures(
#   cisTopicObject, 
#   method = "tSNE", 
#   target = "region", 
#   topic_contr = NULL, 
#   colorBy = c("nCells"), 
#   cex.legend = 0.8, 
#   factor.max = 0.75, 
#   dim = 2, 
#   legend = TRUE, 
#   col.low = "darkgreen", 
#   col.mid = "yellow", 
#   col.high = "brown1", 
#   intervals = 10
# )
# 
# par(mfrow = c(2, 5))
# plotFeatures(
#   cisTopicObject,
#   method = "tSNE", 
#   target = "region", 
#   topic_contr = "Z-score", 
#   colorBy = NULL, 
#   cex.legend = 0.8, 
#   factor.max = 0.75, 
#   dim = 2, 
#   legend=TRUE, 
#   col.low = "darkgreen", 
#   col.mid = "yellow", 
#   col.high = "brown1", 
#   intervals = 10
# )
# 
# # Obtain signatures (if it has not been run before)
# path_to_signatures <- paste0(pathTo10X, "Bulk_peaks/")
# Bulk_ATAC_signatures <- paste(path_to_signatures, list.files(path_to_signatures), sep="")
# labels  <- gsub("._peaks.narrowPeak", "", list.files(path_to_signatures))
# cisTopicObject <- getSignaturesRegions(cisTopicObject, Bulk_ATAC_signatures, labels=labels, minOverlap = 0.4)
# 
# # To only keep unique peaks per signature
# cisTopicObject@signatures <- llply(1:length(cisTopicObject@signatures), function (i) cisTopicObject@signatures[[i]][-which(cisTopicObject@signatures[[i]] %in% unlist(as.vector(cisTopicObject@signatures[-i])))]) 
# names(cisTopicObject@signatures) <- labels
# 
# signaturesHeatmap(cisTopicObject)
################################################################################
################################################################################
################################################################################
# # Read the cisTopic object to a file
# cisTopicObject <- readRDS(
#   file = here::here("data/processed/scenic/00_cisTopicObject.rds")
# )
################################################################################
################################################################################
################################################################################

# seqlevelsStyle(TxDb.Dmelanogaster.UCSC.dm6.ensGene) <- "UCSC"
# seqlevelsStyle(cisTopicObject@region.ranges) <- "UCSC"
# 
# cisTopicObject <- annotateRegions(
#   cisTopicObject, 
#   txdb = TxDb.Dmelanogaster.UCSC.dm6.ensGene, 
#   annoDb = org.Dm.eg.db
# )
# 
# # Plot heatmap with annotated regions
# par(mfrow = c(1, 1))
# signaturesHeatmap(
#   cisTopicObject, 
#   selected.signatures = "annotation"
# )
# plotFeatures(
#   cisTopicObject, 
#   method = "tSNE", 
#   target = "region", 
#   topic_contr = NULL, 
#   colorBy = c("annotation"), 
#   cex.legend = 0.8, 
#   factor.max = 0.75,
#   dim = 2, 
#   legend = TRUE, 
#   intervals = 20
# )
# 
# 
# 
# cisTopicObject <- GREAT(
#   cisTopicObject, 
#   genome = "hg19", 
#   fold_enrichment = 2, 
#   geneHits = 1, 
#   sign = 0.05, 
#   request_interval = 10
# )
# 
# 
# ontologyDotPlot(
#   cisTopicObject, 
#   top = 5,
#   topics = c(1, 3, 13, 26), 
#   var.y = "name", 
#   order.by = "Binom_Adjp_BH"
# )
# 
# 
# cisTopicObject <- binarizedcisTopicsToCtx(
#   cisTopicObject, 
#   genome = "hg19"
# )
# 
# cisTopicObject <- scoredRegionsToCtx(
#   cisTopicObject, 
#   genome = "hg19"
# )
# 
# 
# 
# pathToFeather <- "hg19-regions-9species.all_regions.mc8nr.feather"
# cisTopicObject <- topicsRcisTarget(cisTopicObject, genome="hg19", pathToFeather, reduced_database=FALSE, nesThreshold=3, rocthr=0.005, maxRank=20000, nCores=5)
# 
# 
# Topic1_motif_enr <- cisTopicObject@binarized.RcisTarget[[1]]
# DT::datatable(Topic1_motif_enr[,-c("enrichedRegions", "TF_lowConf"), 
#                                with=FALSE], escape = FALSE, filter="top", 
#               options=list(pageLength=5))
# 
# 
# 
# Topic3_motif_enr <- cisTopicObject@binarized.RcisTarget[[3]]
# DT::datatable(Topic3_motif_enr[,-c("enrichedRegions", "TF_lowConf"), with=FALSE], 
#               escape = FALSE, filter="top", options=list(pageLength=5))
# 
# cisTopicObject <- getCistromes(cisTopicObject, annotation = "Both", nCores=5)
# 
# 
# # Compute AUC rankings based on the predictive distribution
# pred_matrix <- predictiveDistribution(cisTopicObject)
# 
# aucellRankings <- AUCell_buildRankings(
#   pred_matrix, 
#   plot = FALSE,
#   verbose = FALSE
# )
# 
# 
# # Formation of cistromes -------------------------------------------------------
# # A cistrome is set of sequences enriched for motifs linked to a certain 
# # transcription factor
# cisTopicObject <- getCistromeEnrichment(
#   cisTopicObject, 
#   topic = 1, 
#   TFname = "SPI1", 
#   aucellRankings = aucellRankings, 
#   aucMaxRank = 0.05*nrow(aucellRankings), 
#   plot = FALSE
# )
# cisTopicObject <- getCistromeEnrichment(
#   cisTopicObject, 
#   topic = 3, 
#   TFname = "SPI1",
#   aucellRankings = aucellRankings, 
#   aucMaxRank = 0.05*nrow(aucellRankings),
#   plot = FALSE
# )
# 
# par(mfrow = c(1, 2))
# plotFeatures(cisTopicObject, method="tSNE", target="cell", topic_contr=NULL, 
#              colorBy=c("Topic1_SPI1 (675p)","Topic3_SPI1 (380p)"), cex.legend = 0.8, 
#              factor.max=.75, dim=2, legend=TRUE, intervals=10)
# 
# 
# # # Gene accesibility scores -----------------------------------------------------
# # region2gene <- cisTopicObject@region.data[,"SYMBOL", drop=FALSE]
# # region2gene <- split(region2gene, region2gene[,"SYMBOL"]) 
# # region2gene <- lapply(region2gene, rownames) 
# #
# # # From pretrained Garnett"s model markers on pBMC dataset from 10X (Pliner et al, 2019)
# # selectedGenes <- c("CD34", "THY1", "ENG", "KIT", "PROM1", #CD34+
# #                    "NCAM1", "FCGR3A", #NKcells
# #                    "CD14", "FCGR1A", "CD68", "S100A12", #Monocytes
# #                    "CD19", "MS4A1", "CD79A", #Bcells
# #                    "CD3D", "CD3E", "CD3G", #Tcells
# #                    "CD4", "FOXP3", "IL2RA", "IL7R", #CD4 Tcell
# #                    "CD8A", "CD8B", #CD8 Tcell
# #                    "IL3RA", "CD1C", "BATF3", "THBD", "CD209" #Dendritic cells
# # )
# # region2gene_subset <- region2gene[which(names(region2gene) %in% selectedGenes)]
# # predMatSumByGene <- sapply(region2gene_subset, 
# #                            function(x) apply(pred.matrix[x,, drop=FALSE], 2, sum))
# # rownames(predMatSumByGene) <- cisTopicObject@cell.names
# # 
# # 
# # # Add to cell data
# # cisTopicObject <- addCellMetadata(cisTopicObject, predMatSumByGene)
# #
# #
# # # Make plot with acessibility scores for specific genes
# # par(mfrow = c(1, 2))
# # plotFeatures(
# #   cisTopicObject, 
# #   method = "tSNE", 
# #   target = "cell", 
# #   topic_contr = NULL, 
# #   colorBy = c("KIT", "PROM1", "NCAM1", "FCGR3A","CD14","S100A12",
# #               "MS4A1", "CD79A", "CD3D", "CD3E", "CD4", "CD8A"), 
# #              
# #   cex.legend = 0.8, 
# #   factor.max = 0.75,
# #   dim = 2, 
# #   legend = TRUE, 
# #   intervals = 10
# # )
# 
# 



# Save the cisTopic object to a file
saveRDS(
  cisTopicObject,
  file = here::here("data/processed/scenic/00_cisTopicObject.rds")
)

# Convert cisTopic object for SCENIC+ in python --------------------------------

# Load library
library(arrow)

# Load data
cisTopicObject <- readRDS(
  here::here("data/processed/scenic/00_cisTopicObject.rds")
)

# Create cell_topic.feather file
modelMatSelection(cisTopicObject, "cell", "Probability") |> 
  as.data.frame() |> 
  write_feather(here::here("data/processed/scenic/cell_topic.feather"))

# Create cell_topic.feather file
modelMatSelection(cisTopicObject, "region", "Probability") |> 
  as.data.frame() |> 
  write_feather(here::here("data/processed/scenic/region_topic.feather"))

# Save the count matrix
cisTopicObject@count.matrix |> 
  as.data.frame() |> 
  write_feather(here::here("data/processed/scenic/count_matrix.feather"))


