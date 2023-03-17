# ------------------------------------------------------------------------------
# Acclimation Phenotypes for Guest Lecture
# October 20, 2022
# TS O'Leary
# ------------------------------------------------------------------------------

# Load libraries
require(tidyverse)

# Load data
dat <- read_csv(here::here("pheno/data/embryo_acc_pheno.csv"))
dat_new <- read_csv(here::here("pheno/data/acclimation_survival_new_geno_summer_2022.csv"))
dat_ph <- read_csv(here::here("pheno/data/acclimation_pupation_new_geno_summer_2022.csv"))
dat_ws <- read_csv(here::here("pheno/data/acclimation_walking_speed_new_geno_summer_2022.csv"))


# Plot Canton S survival data
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
  theme(axis.line.y = element_line(color = "grey60"),
        axis.line.x = element_line(color = "grey100"),
        axis.ticks = element_line(color = "grey60"),
        panel.grid.major.y = element_line(color = "grey70", 
                                          size = 0.5),
        panel.grid.minor.y = element_line(color = "grey90", 
                                          size = 0.5))

# Plot data for individual genotypes
dat %>% 
  mutate(Genotype = factor(Genotype, 
                           levels = c("CantonS", "RFM 48", "RFM 6", 
                                      "VT 8", "VT 10", "Chiapas", 
                                       "Ghana", "Guam", "St. Kitts"))) %>%
  filter(Stage == "Early") %>%
  mutate(Acclimation = paste0(Acclimation, "°C")) %>%
  ggplot(aes(x = Acclimation,
             y = Survival*100,
             fill = Region)) +
  geom_boxplot(width = 0.3,
               color = "grey20",
               outlier.shape = NA) +
  ggbeeswarm::geom_beeswarm(size = 1.5,
                            cex = 2,
                            shape = 21,
                            color = "grey20",
                            alpha = 0.9) +
  scale_y_continuous(name = "Survival after acute heat shock",
                     limits = c(0, 100),
                     labels = function(x) paste0(x, "%")) + 
  theme_classic(base_size = 14) +
  scale_fill_manual(values = c("skyblue", "pink")) +
  labs(x = "Acclimation temperature") +
  theme(axis.line.y = element_line(color = "grey60"),
        axis.line.x = element_line(color = "grey100"),
        axis.ticks = element_line(color = "grey60"),
        panel.grid.major.y = element_line(color = "grey70", 
                                          size = 0.5),
        panel.grid.minor.y = element_line(color = "grey90", 
                                          size = 0.5)) +
  facet_wrap(~ Genotype,
             nrow = 2) 

# Diurnal variation in all these places
bio <- read_csv("~/R/projects/ox_stress/exp_design/lockwood_lab_stocks_bioclim.csv")
bio %>%
  distinct(Locale, .keep_all = TRUE) %>%
  select(Locale, bio2) %>%
  mutate(Mean_Diurnal_Range = bio2/10)



# Analyze results
mod <- dat %>%
  filter(Genotype == "CantonS") %>%
  filter(Temp == 38.75) %>%
  glm(Survival ~ Acclimation, 
      weights = Number_Eggs,
      family = quasibinomial,
      data = .) 

summary(mod)
broom::tidy(mod)

# New genotypes ----------------------------------------------------------------
dat_new %>% 
  mutate(Genotype = factor(Genotype,
                           levels = c("FRMO01", "FRMO02", "JPOK",
                                      "JPSC", "GH169", "GH174"))) %>%
  dplyr::filter(Stage == "Early") %>%
  mutate(Acclimation = paste0(Acclimation, "°C")) %>%
  ggplot(aes(x = Acclimation,
             y = N_Hatched/N_Eggs*100,
             fill = Region)) +
  geom_boxplot(width = 0.3,
               color = "grey20",
               outlier.shape = NA) +
  ggbeeswarm::geom_beeswarm(size = 1.5,
                            cex = 2,
                            shape = 21,
                            color = "grey20",
                            alpha = 0.9) +
  scale_y_continuous(name = "Survival after acute heat shock",
                     limits = c(0, 100),
                     labels = function(x) paste0(x, "%")) + 
  theme_classic(base_size = 14) +
  scale_fill_manual(values = c("skyblue", "pink")) +
  labs(x = "Acclimation temperature") +
  theme(axis.line.y = element_line(color = "grey60"),
        axis.line.x = element_line(color = "grey100"),
        axis.ticks = element_line(color = "grey60"),
        panel.grid.major.y = element_line(color = "grey70", 
                                          size = 0.5),
        panel.grid.minor.y = element_line(color = "grey90", 
                                          size = 0.5)) +
  facet_wrap(~ Genotype, nrow = 3) 

# Egg-to-pupa survival
dat_new %>%
  mutate(Genotype = factor(Genotype,
                           levels = c("FRMO01", "FRMO02", "JPOK",
                                      "JPSC", "GH169", "GH174"))) %>%
  dplyr::filter(Stage == "Early") %>%
  mutate(Acclimation = paste0(Acclimation, "°C")) %>%
  ggplot(aes(x = Acclimation,
             y = N_Pupated/N_Eggs_for_Pupa_Adult*100,
             fill = Region)) +
  geom_boxplot(width = 0.3,
               color = "grey20",
               outlier.shape = NA) +
  ggbeeswarm::geom_beeswarm(size = 1.5,
                            cex = 2,
                            shape = 21,
                            color = "grey20",
                            alpha = 0.9) +
  scale_y_continuous(name = "Egg-to-pupa survival after acute heat shock",
                     limits = c(0, 100),
                     labels = function(x) paste0(x, "%")) + 
  theme_classic(base_size = 14) +
  scale_fill_manual(values = c("skyblue", "pink")) +
  labs(x = "Acclimation temperature") +
  theme(axis.line.y = element_line(color = "grey60"),
        axis.line.x = element_line(color = "grey100"),
        axis.ticks = element_line(color = "grey60"),
        panel.grid.major.y = element_line(color = "grey70", 
                                          size = 0.5),
        panel.grid.minor.y = element_line(color = "grey90", 
                                          size = 0.5)) +
  facet_wrap(~ Genotype, nrow = 3) 

# Egg-to-adult survival
dat_new %>% 
  mutate(Genotype = factor(Genotype,
                           levels = c("FRMO01", "FRMO02", "JPOK",
                                      "JPSC", "GH169", "GH174"))) %>%
  dplyr::filter(Stage == "Early") %>%
  mutate(Acclimation = paste0(Acclimation, "°C")) %>%
  ggplot(aes(x = Acclimation,
             y = N_Eclosed/N_Eggs_for_Pupa_Adult*100,
             fill = Region)) +
  geom_boxplot(width = 0.3,
               color = "grey20",
               outlier.shape = NA) +
  ggbeeswarm::geom_beeswarm(size = 1.5,
                            cex = 2,
                            shape = 21,
                            color = "grey20",
                            alpha = 0.9) +
  scale_y_continuous(name = "Egg-to-adult survival after acute heat shock",
                     limits = c(0, 100),
                     labels = function(x) paste0(x, "%")) + 
  theme_classic(base_size = 14) +
  scale_fill_manual(values = c("skyblue", "pink")) +
  labs(x = "Acclimation temperature") +
  theme(axis.line.y = element_line(color = "grey60"),
        axis.line.x = element_line(color = "grey100"),
        axis.ticks = element_line(color = "grey60"),
        panel.grid.major.y = element_line(color = "grey70", 
                                          size = 0.5),
        panel.grid.minor.y = element_line(color = "grey90", 
                                          size = 0.5)) +
  facet_wrap(~ Genotype, nrow = 3) 


# Pupation height
dat_ph %>%
  mutate(Genotype = factor(Genotype,
                                levels = c("FRMO01", "FRMO02", "JPOK",
                                           "JPSC", "GH169", "GH174"))) %>%
  dplyr::filter(Stage == "Early") %>%
  mutate(Acclimation = paste0(Acclimation, "°C")) %>%
  ggplot(aes(x = Acclimation,
             y = Height,
             fill = Region)) +
  geom_boxplot(width = 0.3,
               color = "grey20",
               outlier.shape = NA) +
  ggbeeswarm::geom_beeswarm(size = 1.5,
                            cex = 2,
                            shape = 21,
                            color = "grey20",
                            alpha = 0.9) +
  scale_y_continuous(name = "Pupation height") +
  theme_classic(base_size = 14) +
  scale_fill_manual(values = c("skyblue", "pink")) +
  labs(x = "Acclimation temperature") +
  theme(axis.line.y = element_line(color = "grey60"),
        axis.line.x = element_line(color = "grey100"),
        axis.ticks = element_line(color = "grey60"),
        panel.grid.major.y = element_line(color = "grey70", 
                                          size = 0.5),
        panel.grid.minor.y = element_line(color = "grey90", 
                                          size = 0.5)) +
  facet_wrap(~ Genotype, nrow = 3) 

# Walking speed
dat_ws %>%
  dplyr::filter(Stage == "Early") %>%
  mutate(Acclimation = paste0(Acclimation, "°C")) %>%
  ggplot(aes(x = Acclimation,
             y = Speed,
             fill = Region)) +
  geom_boxplot(width = 0.3,
               color = "grey20",
               outlier.shape = NA) +
  ggbeeswarm::geom_beeswarm(size = 1.5,
                            cex = 2,
                            shape = 21,
                            color = "grey20",
                            alpha = 0.9) +
  scale_y_continuous(name = "Pupation height") +
  theme_classic(base_size = 14) +
  scale_fill_manual(values = c("skyblue", "pink")) +
  labs(x = "Acclimation temperature") +
  theme(axis.line.y = element_line(color = "grey60"),
        axis.line.x = element_line(color = "grey100"),
        axis.ticks = element_line(color = "grey60"),
        panel.grid.major.y = element_line(color = "grey70", 
                                          size = 0.5),
        panel.grid.minor.y = element_line(color = "grey90", 
                                          size = 0.5)) +
  facet_wrap(~ Genotype, nrow = 3) 
