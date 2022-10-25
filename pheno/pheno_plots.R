# ------------------------------------------------------------------------------
# Acclimation Phenotypes for Guest Lecture
# October 20, 2022
# TS O'Leary
# ------------------------------------------------------------------------------

# Load libraries
require(tidyverse)

# Load data
dat <- read_csv(here::here("pheno/embryo_acc_pheno.csv"))

# Plot data
dat %>%
  filter(Genotype == "CantonS") %>%
  filter(Stage == "Early") %>%
  mutate(Acclimation = paste0(Acclimation, "°C")) %>%
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
  scale_y_continuous(name = "Survival after acute heat shock",
                     limits = c(0, 100),
                     labels = function(x) paste0(x, "%")) + 
  theme_classic(base_size = 14) +
  labs(x = "Acclimation temperature") +
  theme(axis.line = element_line(color = "grey60"),
        axis.ticks = element_line(color = "grey60"),
        panel.grid.major.y = element_line(color = "grey70", 
                                          size = 0.5),
        panel.grid.minor.y = element_line(color = "grey90", 
                                          size = 0.5))

# 
dat %>%
  #filter(Genotype == "CantonS") %>%
  filter(Stage == "Early") %>%
  mutate(Acclimation = paste0(Acclimation, "°C")) %>%
  ggplot(aes(x = Acclimation,
             y = Survival*100,
             fill = Region)) +
  geom_boxplot(width = 0.3,
               position = position_dodge(width = 0.5),
               color = "grey20",
               outlier.shape = NA) +
  scale_y_continuous(name = "Survival after acute heat shock",
                     limits = c(0, 100),
                     labels = function(x) paste0(x, "%")) + 
  theme_classic(base_size = 16) +
  labs(x = "Acclimation temperature") +
  theme(axis.line = element_line(color = "grey60"),
        axis.ticks = element_line(color = "grey60"),
        panel.grid.major.y = element_line(color = "grey70", 
                                          size = 0.5),
        panel.grid.minor.y = element_line(color = "grey90", 
                                          size = 0.5))

# Analyze results
mod <- dat %>%
  filter(Genotype == "CantonS") %>%
  filter(Temp == 38.75) %>%
  glm(Survival ~ Acclimation, 
      family = binomial,
      data = .)

summary(mod)
str(summary(mod))

TukeyHSD(mod)
