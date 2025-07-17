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

# Proportion hatched
dat |>
  group_by(acc_temp) |>
  summarise(total_hatched = sum(n_hatched),
            total_eggs = sum(n_eggs)) |> 
  mutate(prop_hatched = total_hatched/total_eggs)

# Acclimation effect
mod <- dat |> 
  #mutate(acc_temp = factor(acc_temp)) |> 
  with(glm(n_hatched/n_eggs ~ acc_temp,
      weights = n_eggs,
      family = quasibinomial())) |> 
  broom::tidy()

# Print out model results
mod

# Save results of model
saveRDS(mod, here::here("output/pheno/mod.rds"))

# Test pairwise groups against eachother ----
emmeans::emmeans(
  glm(n_hatched/n_eggs ~ acc_temp,
      weights = n_eggs,
      family = quasibinomial(),
      data = dat |> mutate(acc_temp = factor(acc_temp))), ~ acc_temp) |>
  pairs()


