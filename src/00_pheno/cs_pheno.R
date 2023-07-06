# ------------------------------------------------------------------------------
# Canton S acclimation phenotype figure
# March 31, 2023
# TS O'Leary
# ------------------------------------------------------------------------------

# Description -----
# Canton S phenotype figure

# Load libraries
library(tidyverse)

# Load data
dat <- readr::read_csv(here::here("data/raw/pheno/embryo_acc_pheno.csv"))

# Filter only Canton S data
dat <- dat |>
  dplyr::filter(Genotype == "CantonS") |>
  dplyr::filter(Stage == "Early")

# Count the number of eggs per acclimation treatment
dat |>
  dplyr::group_by(Acclimation) |>
  dplyr::summarise(total_eggs = sum(Number_Eggs))

# # Acclimation effect
# dat |>
#   glm(Survival ~ Acclimation, 
#       data = .,
#       weights = Number_Eggs,
#       family = binomial(link = "logit"))
