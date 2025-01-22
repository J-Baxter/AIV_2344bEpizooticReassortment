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
library(broom)
library(broom.mixed)
library(tidybayes)
library(bayesplot)
library(emmeans)
library(marginaleffects)
library(magrittr)
library(ggmcmc)

# User functions
source('./scripts/figure_scripts/plot_settings.R')

############################################## DATA ################################################
combined_data <- read_csv('./2024Aug18/treedata_extractions/2024-09-20_combined_data.csv')
summary_data <- read_csv('./2024Aug18/treedata_extractions/summary_reassortant_metadata_20240904.csv') %>%
  select(-c(cluster_label,
            clade)) 


############################################## MAIN ################################################

# Data preprocessing
diffusion_data <- combined_data %>%
  
  # select variables of interes
  select(
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
  left_join(summary_data %>% select(c(cluster_profile, 
                                      host_richness)),
            by = join_by(cluster_profile)) %>%
  
  # season (breeding, migrating_spring, migrating_autumn, overwintering)
  mutate(collection_month = date_decimal(TMRCA) %>% format(., "%m") %>% as.integer(),
         season = case_when(collection_month %in% c(12,1,2) ~ 'overwintering', 
                            collection_month %in% c(3,4,5)  ~ 'migrating_spring', # Rename to spring migration
                            collection_month %in% c(6,7,8)  ~ 'breeding', 
                            collection_month %in% c(9,10,11)  ~ 'migrating_autumn' # Rename to autumn migration
         )) %>%
  
  select(weighted_diff_coeff, 
         median_anseriformes_wild,
         median_charadriiformes_wild,
         segment, 
         collection_regionname, 
         season)



# Model Formula
# We assume a hurdle lognormal model, in which the hurdle process is determined the season in which 
# the TMRCA is estimated to be, and the lognormal component is determined by persistence
# in anseriformes and charadriiformes. Both model components are conditional on segment from which 
# the measurement is taken and region of origin

diffusion_formula <- bf(weighted_diff_coeff ~ 1 + median_anseriformes_wild +  median_charadriiformes_wild + collection_regionname + 
                          (1|segment),
                        hu ~ 1 + collection_regionname +  season +(1|segment ))


# Define Priors
diffusionmodel1_priors <- get_prior(diffusion_formula,
                                    data = diffusion_data,
                                    family = hurdle_lognormal()) 

diffusionmodel1_priors$prior[c(1,12)] <- "normal(0,5)"

#diffusionmodel1_priors$prior[9] <- "student_t(3, 0, 3)"
#diffusionmodel1_priors


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
  family = hurdle_lognormal(),
  prior = diffusionmodel1_priors,
  sample_prior = "only",
  chains = CHAINS,
  threads = 2, 
  cores = CORES, 
  iter = ITER,
  warmup = BURNIN,
  seed = SEED,
  control = list(adapt_delta = 0.95)
)


# Fit model to data
diffusionmodel1_fit <- brm(
  diffusion_formula,
  data = diffusion_data,
  prior = diffusionmodel1_priors,
  family = hurdle_lognormal(),
  chains = CHAINS,
  cores = CORES, 
  threads = 2, 
  iter = ITER,
  warmup = BURNIN,
  seed = SEED,
  control = list(adapt_delta = 0.95))


# Post-fitting checks (including inspection of ESS, Rhat and posterior predictive)
tidy_diffusionmodel1 <- tidy(diffusionmodel1_fit)
posteriorpredictive <-pp_check(diffusionmodel1_fit, ndraws = 500)


# Misc evaluations
performance(diffusionmodel1_fit) # tibble output of model metrics including R2, ELPD, LOOIC, RMSE
plot(diffusionmodel1_fit) # default output plot of brms showing posterior distributions of
prior_summary(diffusionmodel1_fit) #obtain dataframe of priors used in model.

mcmc_diffusion <- ggs(diffusionmodel1_fit) # Warning message In custom.sort(D$Parameter) : NAs introduced by coercion


saveRDS(diffusionmodel1_fit, './saved_models/diffusion_model.rds')
write_csv(diffusion_data, './saved_models/diffusion_model.csv')
############################################## END #################################################
####################################################################################################
####################################################################################################