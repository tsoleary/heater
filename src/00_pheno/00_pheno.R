# ------------------------------------------------------------------------------
# Canton S acclimation and acute heat shock phenotype statistics
# TS O'Leary
# ------------------------------------------------------------------------------

# Load libraries
library(tidyverse)

# Load data
dat <- read_csv(here::here("data/raw/pheno/acc_hs_survival.csv"))

# Count the number of eggs per acclimation treatment
dat |>
  group_by(acc_temp) |>
  summarise(total_eggs = sum(n_eggs))

# Acclimation effect
mod <- dat |> 
  with(glm(n_hatched/n_eggs ~ acc_temp,
      weights = n_eggs,
      family = quasibinomial())) |> 
  broom::tidy()

# Print out model results
mod

# Save results of model
saveRDS(mod, here::here("output/pheno/mod.rds"))
