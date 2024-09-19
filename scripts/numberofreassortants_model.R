####################################################################################################
####################################################################################################
# The aim of this model is to determine whether reassortants emerging from some regions are more
# likely and more severe than those from others

# variables of interest include: 
# 1. region of origin
# 3. persistence time (persistence time increases uncertainty)
# 4. number of sequences per region
# 5. diveresity measure?

# This script implements a workflow for a negative binomial model. 

########################################## DEPENDENCIES ############################################
library(brms)
library(broom)
library(broom.mixed)
library(tidybayes)
library(bayesplot)


########################################### IMPORT DATA ############################################

combined_data <- read_csv()


########################################### FORMAT DATA ############################################


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
  n_reassortants ~ n_sequences + collection.region.name,
  data = grouped_lines,
  family = negbinomial(), sample_prior = "yes",
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
  n_reassortants ~ n_sequences + collection.region.name,
  data = grouped_lines,
  family = negbinomial(),
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

#deprecated 
bayes_mod <- 
cbind.data.frame(collection.region.name = c('europe', 'asia', 'america', 'africa'), n_sequences = 100) %>%
  as_tibble()%>%
  add_epred_draws(bayes_mod, ndraws = 100) %>%
  ggplot(.,  aes(x = .epred, y = collection.region.name)) +
  stat_halfeye() +
  theme_minimal()


cbind.data.frame(collection.region.name = c('europe', 'asia', 'america', 'africa'), n_sequences = 100) %>%
  as_tibble()%>%
  add_epred_draws(bayes_mod, ndraws = 100) %>%
  mean_hdi()


bayes_mod %>%
  emmeans(~ collection.region.name,
          at = list(collection.region.name =  c('europe', 'asia', 'america', 'africa')),
          re_formula = NA) %>%
  contrast(method = 'revpairwise', ref= 'asia', type = 'response') %>%
  gather_emmeans_draws() %>%
  mutate(`.value` = exp(`.value`)) %>% mean_hdi()

ggplot(., aes(x = `.value`, y = contrast, fill = contrast)) + 
  stat_halfeye(point_interval = mean_hdi,.width = c(.90, .95))

reassortant_metadata_formatted %>%
  mutate(collection.date = ymd(collection.date)) %>%
  ggplot(aes(x = collection.date, y = collection.region.name, fill = as.factor(cluster.genome)))+
  geom_density()+
  facet_wrap(.~collection.region.name)