####################################################################################################
####################################################################################################
## Script name: Diffusion Coefficient Model
##
## Purpose of script: A Bayesian mixed model to determine associations between variables obtained 
## from our phylogenetic analysis and difusion coefficient for each reassortant. Variables of interest
## include: 1. evolutionary rates, 2. persistence time, 3. region of origin, 4. persistence time in 
## wild birds, and 6. total number of species jumps
##
## This script implements a workflow for a zero-inflated lognormal model. 
##
## Date created: 2024-xx-xx
##
##
########################################## SYSTEM OPTIONS ##########################################
#options(scipen = 6, digits = 7) 
memory.limit(30000000) 


########################################## DEPENDENCIES ############################################
# Packages
library(tidyverse)
library(magrittr)
library(brms)
library(rstan)
library(broom)
library(broom.mixed)
library(tidybayes)
library(bayesplot)
library(emmeans)
library(marginaleffects)
library(ggmcmc)
library(performance)
library(DHARMa)
# User functions
#source('./scripts/figure_scripts/plot_settings.R')

############################################## DATA ################################################
combined_data <- read_csv('./2024Aug18/treedata_extractions/2024-09-20_combined_data.csv')
summary_data <- read_csv('./2024Aug18/treedata_extractions/summary_reassortant_metadata_20240904.csv') %>%
  select(-c(cluster_label,
            clade)) 
reassortant_ancestral_changes <- read_csv('./reassortant_ancestral_changes.csv')


############################################## MAIN ################################################

# Data preprocessing
diffusion_data <- combined_data %>%
  
  # select variables of interes
  dplyr::select(
    segment,
    cluster_profile,
    TMRCA,
    group2,
    weighted_diff_coeff,
    original_diff_coeff,
    evoRate,
    persist.time,
    collection_regionname,
    host_simplifiedhost,
    count_cross_species,
    starts_with('median'),
    starts_with('max')) %>%
  
  # Substitute NA values in diffusion coefficient with 0
  mutate(across(where(is.double), .fns = ~ replace_na(.x, 0))) %>%
  drop_na(collection_regionname) %>%
  
  rename_with(~gsub('-', '_', .x)) %>%
  mutate(collection_regionname = case_when(grepl('europe', collection_regionname) ~ 'europe',
                                           grepl('africa', collection_regionname) ~ 'africa',
                                           grepl('asia', collection_regionname) ~ 'asia',
                                           grepl('(central|northern) america', collection_regionname) ~ 'central & northern america',
                                           .default = collection_regionname
  )) %>%
  
  mutate(collection_regionname = factor(collection_regionname, levels = c('asia', 'africa', 'europe', 'central & northern america'))) %>%
  filter(!grepl('\\+', host_simplifiedhost)) %>%
  
  # join host richness
  left_join(summary_data %>% dplyr::select(c(cluster_profile, 
                                      host_richness)),
            by = join_by(cluster_profile)) %>%
  
  # season (breeding, migrating_spring, migrating_autumn, overwintering)
  mutate(collection_month = date_decimal(TMRCA) %>% format(., "%m") %>% as.integer(),
         collection_year = date_decimal(TMRCA) %>% format(., "%Y") %>% as.integer() %>% factor(levels = c(2017,2018,2019,2020,2021,2022,2023)),
         collection_season = case_when(collection_month %in% c(12,1,2) ~ 'overwintering', 
                            collection_month %in% c(3,4,5)  ~ 'migrating_spring', # Rename to spring migration
                            collection_month %in% c(6,7,8)  ~ 'breeding', 
                            collection_month %in% c(9,10,11)  ~ 'migrating_autumn' # Rename to autumn migration
         )) %>%
  
  # persistence time proportions
  mutate(across(c('median_anseriformes_wild', 
                  'median_charadriiformes_wild'), .fns= ~ ifelse(.x/persist.time > 1, 1, .x/persist.time ), .names = '{.col}_prop')) %>%
  
  
  # binary
  mutate(across(c('persist.time', 
                  'count_cross_species'), .fns= ~ log1p(.x), .names = '{.col}_log1p'))  %>%
  
  # Ancestral (WRT HA) Changes
  left_join(reassortant_ancestral_changes %>% select(cluster_profile, segments_changed,cluster_class)) %>%
  replace_na(list(segments_changed = 0)) %>%
  
  dplyr::select(cluster_profile,
                weighted_diff_coeff, 
                persist.time_log1p,
                collection_year,
                count_cross_species_log1p,
                median_anseriformes_wild_prop,
                median_charadriiformes_wild_prop,
                host_simplifiedhost,
                segments_changed,
                segment,
                cluster_class,
                collection_regionname) %>%

  filter(weighted_diff_coeff > 0) %>%
  filter(segment == 'ha') %>%
  drop_na()
  



# Model Formula
# We assume a hurdle lognormal model, in which the hurdle process is determined the season in which 
# the TMRCA is estimated to be, and the lognormal component is determined by persistence
# in anseriformes and charadriiformes. Both model components are conditional on segment from which 
# the measurement is taken and region of origin

#diffusion_formula <- bf(weighted_diff_coeff ~ 1 + median_anseriformes_wild +  median_charadriiformes_wild + collection_regionname + 
                       #   (1|segment),
                       # hu ~ 1 + collection_regionname +  season +(1|segment ))

#diffusion_formula <- bf(weighted_diff_coeff|trunc(lb=0) ~ 1 + median_anseriformes_wild +  median_charadriiformes_wild + collection_regionname + (1|segment) + (collection_regionname|collection_year/collection_season))
#diffusion_formula <- bf(weighted_diff_coeff ~ 0 +collection_regionname + median_anseriformes_wild + median_charadriiformes_wild +  collection_season + (1|segment),
                     #   shape ~ 0 +  collection_regionname + (1 | segment))
#diffusion_formula <- bf(weighted_diff_coeff ~ 0 + median_anseriformes_wild +  median_charadriiformes_wild + collection_regionname + collection_season +  (1|segment))



diffusion_formula <- bf(weighted_diff_coeff ~ 0 + collection_regionname + 
                          count_cross_species_log1p + 
                          median_anseriformes_wild_prop +
                          median_charadriiformes_wild_prop +
                          
                          collection_regionname:median_anseriformes_wild_prop + 
                          collection_regionname:median_charadriiformes_wild_prop +
                          collection_regionname:persist.time_log1p +
                          (1|collection_year),
                        shape ~  0 + collection_regionname )



# Define Priors
 get_prior(diffusion_formula,
           data = diffusion_data,
           family = Gamma(link = "log")) 


diffusionmodel1_priors <- c(set_prior("normal(0, 5)", class = 'b'),
                            set_prior('normal(13,2)', class = 'b', coef = 'collection_regionnameafrica'),
                            set_prior('normal(13,1.5)', class = 'b', coef = 'collection_regionnameasia'),
                            set_prior('normal(13,1.5)', class = 'b', coef = 'collection_regionnamecentral&northernamerica'),
                            set_prior('normal(12,1.5)', class = 'b', coef = 'collection_regionnameeurope'),
                            set_prior('normal(0,2)', class = 'b', coef = 'collection_regionnameafrica:median_anseriformes_wild_prop'),
                            set_prior('normal(0,2)', class = 'b', coef = 'collection_regionnameafrica:median_charadriiformes_wild_prop'),
                            set_prior('normal(0,2)', class = 'b', coef = 'collection_regionnameafrica:persist.time_log1p'),
                            set_prior('normal(0,5)',  dpar = 'shape'),
                            set_prior('exponential(0.5)', class = 'sd'))


# Set MCMC Options
CHAINS <- 4
CORES <- 4
ITER <- 4000
BURNIN <- ITER/10 # Discard 10% burn in from each chain
SEED <- 4472



# Prior Predictive Checks 
# Note that due to weakly informative priors, the prior predictive checks reveal little about the 
# ability of our priors to describe the data. 

diffusionmodel1_prior <- brm(
  diffusion_formula,
  data = diffusion_data,
  family = Gamma(link = "log"),
  prior = diffusionmodel1_priors,
  sample_prior = "only",
  chains = 2,
  threads = 2, 
  cores = 2, 
  iter = ITER,
  warmup = BURNIN,
  seed = SEED,
  control = list(adapt_delta = 0.95)
)


# Fit model to data
diffusionmodel1_fit_gamma_19<- brm(
  diffusion_formula,
  data = diffusion_data,
  prior = diffusionmodel1_priors,
  family = Gamma(link = "log"),
  chains = CHAINS,
  cores = CORES, 
  threads = 2, 
  backend = "cmdstanr",
  iter = ITER,
  warmup = BURNIN,
  seed = SEED,
  control = list(adapt_delta = 0.95))


# Post-fitting checks (including inspection of ESS, Rhat and posterior predictive)
performance(diffusionmodel1_fit_gamma_19)


saveRDS(diffusionmodel1_fit, './saved_models/diffusion_model_2.rds')
write_csv(diffusion_data, './saved_models/diffusion_model.csv')
############################################## END #################################################
####################################################################################################
####################################################################################################