# ------------------------------------------------------------------------------
# Driver script for all Seurat
# TS O'Leary
# ------------------------------------------------------------------------------

# Load libraries
library(tidyverse)
library(Seurat)
library(Signac)

# Source all files

dir <- here::here("src/03_seurat")

files <- list.files(dir, pattern = "[0-9]{2}")

for (i in files){
  source(here::here(dir, files[i]))
}
