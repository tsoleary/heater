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
<<<<<<< HEAD
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
=======
dat <- read_csv(here::here("data/raw/pheno/embryo_acc_pheno.csv"))

# Filter only Canton S data
dat <- dat |>
  filter(Genotype == "CantonS") |>
  filter(Stage == "Early")

# Count the number of eggs per acclimation treatment
dat |>
  group_by(Acclimation) |>
  summarise(total_eggs = sum(Number_Eggs))

# Acclimation effect
dat |>
  glm(Survival ~ Acclimation, 
      data = .,
      weights = Number_Eggs,
      family = binomial(link = "logit"))
>>>>>>> a8de26aef8a500cdbd82fa277045d3e778c75e21
