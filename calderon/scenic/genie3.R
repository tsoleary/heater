# ------------------------------------------------------------------------------
# GENIE3 Run of Calderon dat_comb data
# November 01, 2022
# TS O'Leary
# ------------------------------------------------------------------------------

# Load libraries
require(SCENIC)

# Load data
scenicOptions <- readRDS(here::here("calderon/scenic/int/scenicOptions.rds"))
exprMat_filtered_log <- readRDS(here::here("calderon/scenic/int/exprMat_filtered_log.rds"))

# Set the number of cores to match the cores allocated in slurm
scenicOptions@settings$nCores <- 10

# Sample down to 3k cells
set.seed(1)
mat <- exprMat_filtered_log[, sample(1:ncol(exprMat_filtered_log), size = 3000)]

# Run GENIE3
runGenie3(mat, scenicOptions)