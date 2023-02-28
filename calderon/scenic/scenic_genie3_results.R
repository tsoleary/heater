# ------------------------------------------------------------------------------
# SCENIC GENIE3 results
# November 07, 2022
# TS O'Leary
# ------------------------------------------------------------------------------

# Load libraries
require(tidyverse)
require(SCENIC)

# Load data
scenicOptions <- readRDS(here::here("calderon/scenic/int/scenicOptions.rds"))
exprMat_filtered_log <- readRDS(here::here("calderon/scenic/int/exprMat_filtered_log.rds"))
# Sample down to 3k cells
set.seed(1)
mat <- exprMat_filtered_log[, sample(1:ncol(exprMat_filtered_log), size = 3000)]

# For a quick run set this 
scenicOptions@settings$dbs <- scenicOptions@settings$dbs["1kb"]


scenicOptions <- runSCENIC_1_coexNetwork2modules(scenicOptions)
scenicOptions <- runSCENIC_2_createRegulons(scenicOptions, 
                                            coexMethod = c("top5perTarget"))
scenicOptions <- runSCENIC_3_scoreCells(scenicOptions, 
                                        mat)


# Looking at the results a bit ----
regulons <- loadInt(scenicOptions, "regulons")
regulons[c("Dlx5", "Irf1")]


regulons <- loadInt(scenicOptions, "aucell_regulons")
head(cbind(onlyNonDuplicatedExtended(names(regulons))))

regulonTargetsInfo <- loadInt(scenicOptions, "regulonTargetsInfo")
tableSubset <- regulonTargetsInfo[TF=="Stat6" & highConfAnnot==TRUE]
viewMotifs(tableSubset, options=list(pageLength=5)) 



# Load Transcription Factor modules
tf_mod <- readRDS(here::here("calderon/scenic/int/1.6_tfModules_asDF.Rds"))

# Total number of TFs included
length(unique(tf_mod$TF))

# Number of genes regulated by each TF in descending order
tf_mod %>%
  group_by(TF) %>%
  tally() %>%
  arrange(desc(n))

# Histogram of number of genes regulated by each TF
tf_mod %>%
  group_by(TF) %>%
  tally() %>%
  ggplot() +
  geom_histogram(aes(x = n),
                 color = "grey20",
                 fill = "grey50",
                 bins = 40) +
  labs(x = "Number of genes regulated by a transcription factor",
       y = "Count") +
  theme_classic()



