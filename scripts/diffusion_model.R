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


########################################## DEPENDENCIES ############################################
library(brms)
library(broom)
library(broom.mixed)
library(tidybayes)
library(bayesplot)


########################################### IMPORT DATA ############################################

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

  
################################### INITIAL EXPLANANTORY MODELS ####################################





####################################### START BRMS PIPELINE ########################################

# Set Priors


# Prior Predictive Checks 


# Fit Model


# Posterior Predictive Checks


# Posterior Predictive Checks


# Extract Marginal/Conditional Effects




model_hurdle <- brms::brm(
  bf(weighted_diff_coeff ~ group2,
     hu ~ 1),
  data = model_data,
  family = hurdle_lognormal(),
  chains = 2,
  iter = 5000,
  warmup= 500,
  seed =123
)
  mutate(log_weighted_diff_coeff = log1p(weighted_diff_coeff))



# logged Data
model_data %>%
  ggplot() +
  geom_histogram(aes(x = log1p(weighted_diff_coeff), fill = weighted_diff_coeff>0)) +
  scale_fill_brewer(palette = 'Dark2', 'Is Zero') +
  scale_x_continuous('Weighted Diffusion Coefficient') +
  theme_minimal()


# prior predictive checks hurdle model
model_hurdle_priorpredictive <- brms::brm(
  bf(weighted_diff_coeff ~ host_simplifiedhost,
     hu ~ 1),
  data = model_data,
  family = hurdle_lognormal(),
  sample_prior = "yes",
  chains = 2,
  cores = 2, 
  iter = 2000,
  warmup= 200,
  seed =123
)

pred <- posterior_predict(model_hurdle_priorpredictive)

color_scheme_set("green")
ppc_dens_overlay(y = log1p(model_data$weighted_diff_coeff),
                 yrep = log1p(pred[1:50,]))

ppc_dens_overlay()

tidy(model_hurdle)

pred <- posterior_predict(model_hurdle)

color_scheme_set("green")
ppc_dens_overlay(y = log1p(model_data$weighted_diff_coeff),
                 yrep = log1p(pred[1:10,]))

####################################################################################################
####################################################################################################