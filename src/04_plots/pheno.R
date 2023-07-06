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
dat <- read_csv(here::here("data/raw/pheno/embryo_acc_pheno.csv"))

# Filter only Canton S data
dat <- dat |>
  dplyr::filter(Genotype == "CantonS") |>
  dplyr::filter(Stage == "Early")

# Count the number of eggs per acclimation treatment
dat |>
  group_by(Acclimation) |>
  summarise(total_eggs = sum(Number_Eggs))

# Plot Canton S Survival
dat |>
  mutate(Acclimation = paste0(Acclimation, "°C")) |>
  ggplot(aes(x = Acclimation,
             y = Survival*100)) +
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
  scale_y_continuous(name = "Survival after acute heat shock",
                     limits = c(0, 100),
                     expand = expansion(mult = c(0, 0.05)),
                     labels = function(x) paste0(x, "%")) + 
  cowplot::theme_minimal_hgrid(font_family = "Myriad Pro")


# Save two sizes
ggsave(here::here("output/figs/pheno/cantonS_survival.png"), 
       height = 16, 
       width = 24,
       units = "cm")
ggsave(here::here("output/figs/pheno/cantonS_survival_small.png"), 
       height = 12, 
       width = 18,
       units = "cm")
