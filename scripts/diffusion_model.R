####################################################################################################
####################################################################################################
# The aim of this model is to determine associations between variables obtained from our phylogenetic
# analysis and difusion coefficient for each reassortant

# variables of interest include: 
# 1. evolutionary rates
# 2. persistence time
# 3. region of origin
# 4. persistence time in wild birds
# 5. persistence time in domestic birds
# 6. total number of species jumps

# This script implements a workflow for a zero-inflated lognormal model. Initial exploratory analysis 
# is conducted using OLS regression, thereafter progressing to Bayesian analysis using BRMS.

########################################## DEPENDENCIES ############################################
library(brms)
library(broom)
library(broom.mixed)
library(tidybayes)
library(bayesplot)


########################################### IMPORT DATA ############################################

combined_data <- read_csv()


########################################### FORMAT DATA ############################################

model_data <- combined_data %>%
  
  # select variables of interes
  select(
    cluster_profile,
    group2,
    weighted_diff_coeff,
    original_diff_coeff,
    evoRate,
    persist.time,
    collection_regionname,
    host_simplifiedhost,
    count_cross_species) %>%
  
  # Substitute NA values in diffusion coefficient with 0
  replace_na(list(original_diff_coeff = 0,
                  weighted_diff_coeff  = 0)) %>%

  
################################### INITIAL EXPLORATORY MODELS #####################################
# Plot logged data to show zero/non-zero segregation

model_data %>%
  ggplot() +
  geom_histogram(aes(x = log1p(weighted_diff_coeff), fill = weighted_diff_coeff>0)) +
  scale_fill_brewer(palette = 'Dark2', 'Is Zero') +
  scale_x_continuous('Weighted Diffusion Coefficient') +
  theme_minimal()


# OLS model of 




####################################### START BRMS PIPELINE ########################################

# Set Priors
diffusionmodel1_priors <- c()


# Set MCMC Options
CHAINS <- 4
CORES <- 4
ITER <- 4000
BURNIN <- ITER/10 # Discard 10% burn in from each chain
SEED <- 4472


# Prior Predictive Checks 
diffusionmodel1_priorpredictive <- brm(
  bf(weighted_diff_coeff ~ collection_regionname,
     hu ~ 1),
  data = model_data,
  family = hurdle_lognormal(),
  sample_prior = "yes",
  chains = CHAINS,
  cores = CORES, 
  iter = ITER,
  warmup = BURNIN,
  seed = SEED#,
  #opencl = opencl(c(0, 0)) # Enables GPU accelerated computation, remove if not applicable
  )


diffusionmodel1_priorpreds <- posterior_predict(diffusionmodel1_priorpredictive)
n <- sample(1:nrow(diffusionmodel1_priorpreds), 50)

color_scheme_set("green")
ppc_dens_overlay(y = log1p(model_data$weighted_diff_coeff),
                 yrep = log1p(diffusionmodel1_priorpreds[1:10,]))


# Fit Model
diffusionmodel1_fit <- brm(
  bf(weighted_diff_coeff ~ collection_regionname,
     hu ~ 1),
  data = model_data,
  family = hurdle_lognormal(),
  chains = CHAINS,
  cores = CORES, 
  iter = ITER,
  warmup = BURNIN,
  seed = SEED#,
  #opencl = opencl(c(0, 0)) # Enables GPU accelerated computation, remove if not applicable
)


# Posterior Predictive Checks
diffusionmodel1_posteriorpreds <- posterior_predict(diffusionmodel1_fit)
n <- sample(1:nrow(diffusionmodel1_posteriorpreds), 50)

color_scheme_set("green")
ppc_dens_overlay(y = log1p(model_data$weighted_diff_coeff),
                 yrep = log1p(diffusionmodel1_posteriorpreds[1:10,]))


# Extract Model Terms and Marginal/Conditional Effects
diffusionmodel1_fit_tidy <- tidy(model_hurdle)


####################################################################################################
####################################################################################################