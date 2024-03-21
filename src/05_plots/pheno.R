# ------------------------------------------------------------------------------
# Acclimation temperature vs egg hatching success phenotype figure
# TS O'Leary
# ------------------------------------------------------------------------------

# Load libraries
library(tidyverse)

# Output dir
fig_dir <- "output/figs/pheno"

# Load data
dat <- read_csv(here::here("data/raw/pheno/acc_hs_survival.csv"))

# Plot Canton S hatching success after heat shock
dat |>
  mutate(acc_temp = paste0(acc_temp, "°C")) |>
  ggplot(aes(x = acc_temp,
             y = n_hatched/n_eggs*100)) +
  geom_boxplot(width = 0.3,
               color = "grey20",
               outlier.shape = NA) +
  ggbeeswarm::geom_beeswarm(size = 3,
                            cex = 2,
                            shape = 21,
                            fill = "grey50",
                            color = "grey20",
                            alpha = 0.9) +
  xlab("Acclimation temperature") +
  scale_y_continuous(name = "Egg hatching success after acute heat shock",
                     limits = c(0, 100),
                     expand = expansion(mult = c(0, 0.05)),
                     labels = function(x) paste0(x, "%")) + 
  cowplot::theme_minimal_hgrid()

# Save two sizes
ggsave(here::here(fig_dir, "cantonS_survival.png"), 
       height = 16, 
       width = 24,
       units = "cm")
ggsave(here::here(fig_dir, "cantonS_survival.pdf"), 
       height = 16, 
       width = 24,
       units = "cm")
ggsave(here::here(fig_dir, "cantonS_survival_small.png"), 
       height = 12, 
       width = 18,
       units = "cm")
ggsave(here::here(fig_dir, "cantonS_survival_small.pdf"), 
       height = 12, 
       width = 18,
       units = "cm")



# Plot Canton S hatching success after heat shock with spaced axis -------------

# Load data
dat <- read_csv(here::here("data/raw/pheno/acc_hs_survival.csv")) |> 
  mutate(acc_temp_factor = factor(paste0(acc_temp, "°C")))

dat |>
  ggplot(aes(x = acc_temp,
             y = n_hatched/n_eggs*100,
             group = acc_temp_factor)) +
  geom_boxplot(width = 1.5,
               color = "grey20",
               outlier.shape = NA) +
  ggbeeswarm::geom_beeswarm(size = 3,
                            cex = 2,
                            shape = 21,
                            fill = "grey50",
                            color = "grey20",
                            alpha = 0.9) +
  xlab("Acclimation temperature") +
  scale_y_continuous(name = "Egg hatching success after acute heat shock",
                     limits = c(0, 100),
                     expand = expansion(mult = c(0, 0.05)),
                     labels = function(x) paste0(x, "%")) + 
  scale_x_continuous(breaks = c(18, 25, 30),
                     labels = levels(dat$acc_temp_factor)) +
  cowplot::theme_minimal_hgrid()


# Save two sizes
ggsave(here::here(fig_dir, "cantonS_survival_space.png"), 
       height = 16, 
       width = 24,
       units = "cm")
ggsave(here::here(fig_dir, "cantonS_survival_space.pdf"), 
       height = 16, 
       width = 24,
       units = "cm")
ggsave(here::here(fig_dir, "cantonS_survival_small_space.png"), 
       height = 12, 
       width = 18,
       units = "cm")
ggsave(here::here(fig_dir, "cantonS_survival_small_space.pdf"), 
       height = 12, 
       width = 18,
       units = "cm")
