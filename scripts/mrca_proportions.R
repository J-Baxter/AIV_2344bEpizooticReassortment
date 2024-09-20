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

combined_data <- read_csv('./2024Aug18/treedata_extractions/2024-09-20_combined_data.csv')


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
    count_cross_species,
    starts_with('median'),
    starts_with('max')) %>%
  
  # Substitute NA values in diffusion coefficient with 0
  replace_na(list(original_diff_coeff = 0,
                  weighted_diff_coeff = 0,
                  `median_galliformes-domestic` = 0,
                  `median_other-bird` = 0,
                  `median_anseriformes-wild` = 0,
                  `median_galliformes-wild` = 0,
                  median_environment = 0,
                  median_mammal = 0,
                  `median_charadriiformes-wild` = 0,
                  median_human = 0,
                  `median_anseriformes-domestic` = 0)) %>%
  rename_with(~gsub('-', '_', .x))


################################### PLOT HOST ANCESTOR #####################################
combined_data %>%
  mutate(collection_regionname = case_when(grepl('europe', collection_regionname) ~ 'europe',
                                           grepl('africa', collection_regionname) ~ 'africa',
                                           grepl('asia', collection_regionname) ~ 'asia',
                                           grepl('(central|northern) america', collection_regionname) ~ 'central & northern america',
                                           .default = collection_regionname
  )) %>%
  filter(!is.na(collection_regionname)) %>% 
  filter(!grepl('\\+', host_simplifiedhost)) %>%
  gsub()
  ggplot(aes(x = collection_regionname,
             fill = host_simplifiedhost)) +
  
  # Geom objects
  geom_bar(position = "fill") +
  scale_y_continuous('Number of Reassortants', labels = scales::percent) +
  
  # Scales
  scale_x_discrete('Region of Origin', labels = function(x) str_wrap(x, width = 20) %>% str_to_title())+
  scale_fill_brewer('Reassortant Class', direction = -1) +
  
  
  
  # Graphical
  facet_grid(rows = vars(segment)) + 
  theme_minimal() 



#######################################  Plot 2: Region MRCA  ########################################
combined_data %>%
  mutate(collection_regionname = case_when(grepl('europe', collection_regionname) ~ 'europe',
                                           grepl('africa', collection_regionname) ~ 'africa',
                                           grepl('asia', collection_regionname) ~ 'asia',
                                           grepl('(central|northern) america', collection_regionname) ~ 'central & northern america',
                                           .default = collection_regionname
  )) %>%
  filter(!is.na(collection_regionname)) %>%
  ggplot(aes(x = collection_regionname,
             fill = group2)) +
  
  # Geom objects
  geom_bar(position = "fill") +
  scale_y_continuous('Number of Reassortants', labels = scales::percent) +
  
  # Scales
  scale_x_discrete('Region of Origin', labels = function(x) str_wrap(x, width = 20) %>% str_to_title())+
  scale_fill_brewer('Reassortant Class', direction = -1) +
  
  
  
  # Graphical
  facet_grid(rows = vars(segment)) + 
  theme_minimal() 




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
  bf(weighted_diff_coeff ~ collection_regionname + median_anseriformes_wild,
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