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

# User functions
#source('./scripts/figure_scripts/plot_settings.R')

############################################## DATA ################################################
combined_data <- read_csv('./2024Aug18/treedata_extractions/2024-09-20_combined_data.csv')
summary_data <- read_csv('./2024Aug18/treedata_extractions/summary_reassortant_metadata_20240904.csv') %>%
  select(-c(cluster_label,
            clade)) 


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
         collection_year = date_decimal(TMRCA) %>% format(., "%Y") %>% as.integer(),
         collection_season = case_when(collection_month %in% c(12,1,2) ~ 'overwintering', 
                            collection_month %in% c(3,4,5)  ~ 'migrating_spring', # Rename to spring migration
                            collection_month %in% c(6,7,8)  ~ 'breeding', 
                            collection_month %in% c(9,10,11)  ~ 'migrating_autumn' # Rename to autumn migration
         )) %>%
  
  # format persistence times
  
  # binary
  mutate(in_charadriiformes = ifelse(median_charadriiformes_wild >0, '1', '0'),
         in_anseriformes = ifelse(median_anseriformes_wild >0, '1', '0'),
         median_anseriformes_wild_log = log1p(median_anseriformes_wild),
         median_charadriiformes_wild_log = log1p(median_charadriiformes_wild)) %>%
  
  group_by(cluster_profile) %>%
  mutate(wdc_sd = sd(weighted_diff_coeff) ) %>%
  ungroup() %>%
  
  dplyr::select(weighted_diff_coeff, 
                wdc_sd,
                in_charadriiformes,in_anseriformes,
         median_anseriformes_wild,
         median_charadriiformes_wild,
         cluster_profile,
         count_cross_species,
         host_simplifiedhost,
         segment, 
         collection_regionname, 
         collection_year,
         collection_season,
         group2) %>%
  
  filter(weighted_diff_coeff > 0)



# Model Formula
# We assume a hurdle lognormal model, in which the hurdle process is determined the season in which 
# the TMRCA is estimated to be, and the lognormal component is determined by persistence
# in anseriformes and charadriiformes. Both model components are conditional on segment from which 
# the measurement is taken and region of origin

#diffusion_formula <- bf(weighted_diff_coeff ~ 1 + median_anseriformes_wild +  median_charadriiformes_wild + collection_regionname + 
                       #   (1|segment),
                       # hu ~ 1 + collection_regionname +  season +(1|segment ))

#diffusion_formula <- bf(weighted_diff_coeff|trunc(lb=0) ~ 1 + median_anseriformes_wild +  median_charadriiformes_wild + collection_regionname + (1|segment) + (collection_regionname|collection_year/collection_season))
diffusion_formula <- bf(weighted_diff_coeff ~ 0 +collection_regionname + median_anseriformes_wild + median_charadriiformes_wild +  collection_season + (1|segment),
                        shape ~ 0 +  collection_regionname + (1 | segment))
#diffusion_formula <- bf(weighted_diff_coeff ~ 0 + median_anseriformes_wild +  median_charadriiformes_wild + collection_regionname + collection_season +  (1|segment))

int_step <- function(x){
  ifelse(x >0, 1,0)
}
diffusion_formula <- bf(weighted_diff_coeff ~ 0 + collection_regionname + int_step(median_anseriformes_wild) + median_anseriformes_wild + int_step(median_charadriiformes_wild) + median_charadriiformes_wild + int_step(count_cross_species) + log1p(count_cross_species) + 
                          collection_regionname:median_anseriformes_wild + 
                          collection_regionname:median_charadriiformes_wild + (1|segment) ,
                        shape ~  0 + collection_regionname + (1 | collection_regionname))



# Define Priors
diffusionmodel1_priors <- get_prior(diffusion_formula,
                                    data = diffusion_data,
                                    family = Gamma(link = "log")) 

#diffusionmodel1_priors$prior[2:5] <- "normal(0,5)"
#diffusionmodel1_priors$prior[6:11] <- "normal(0,1)"
#diffusionmodel1_priors$prior[16:19] <- "normal(0,5)"

diffusionmodel1_priors$prior[2:5] <- "normal(0,5)"
diffusionmodel1_priors$prior[6:17] <- "normal(0,1)"
diffusionmodel1_priors$prior[21:25] <- "normal(0,5)"


#diffusionmodel1_priors$prior[1:14] <- "normal(0,5)"
#diffusionmodel1_priors$prior[18:22] <- "normal(0,5)"
diffusionmodel1_priors

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
diffusionmodel1_fit_gamma_18<- brm(
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
tidy_diffusionmodel1 <- tidy(diffusionmodel1_fit)
posteriorpredictive <-pp_check(diffusionmodel1_fit, ndraws = 500)
performance(diffusionmodel1_fit_gamma_14)

loo_compare(diffusionmodel1_fit, diffusionmodel2_fit)
# Misc evaluations
# performance(diffusionmodel1_fit) # tibble output of model metrics including R2, ELPD, LOOIC, RMSE
plot(diffusionmodel1_fit) # default output plot of brms showing posterior distributions of
prior_summary(diffusionmodel1_fit) #obtain dataframe of priors used in model.

mcmc_diffusion <- ggs(diffusionmodel1_fit) # Warning message In custom.sort(D$Parameter) : NAs introduced by coercion


saveRDS(diffusionmodel1_fit, './saved_models/diffusion_model_2.rds')
write_csv(diffusion_data, './saved_models/diffusion_model.csv')
############################################## END #################################################
####################################################################################################
####################################################################################################