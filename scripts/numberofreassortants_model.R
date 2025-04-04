####################################################################################################
####################################################################################################
## Script name: Number and Class of Reassortants Model
##
## Purpose of script: to determine whether reassortants emerging from some regions are more
## likely and more severe than those from others. Variables of interest interest include: region of
## origin, persistence time (persistence time increases uncertainty), number of sequences per region
##
## This script implements a workflow for a multivariate negative binomial model. 
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

FormatContinent <- function(dataframe){
  dataframe %<>%
    mutate(collection_regionname = case_when(grepl('europe', 
                                                   collection_regionname, 
                                                   ignore.case = TRUE) ~ 'europe',
                                             
                                             grepl('africa',
                                                   collection_regionname, 
                                                   ignore.case = TRUE) ~ 'africa',
                                             
                                             grepl('asia', 
                                                   collection_regionname, 
                                                   ignore.case = TRUE) ~ 'asia',
                                             
                                             grepl('(central|northern) america', 
                                                   collection_regionname, 
                                                   ignore.case = TRUE) ~ 'central & northern america',
                                             
                                             .default = NA_character_ ))
  
  return(dataframe)
}

############################################## DATA ################################################
combined_data <- read_csv('./2024Aug18/treedata_extractions/2024-09-20_combined_data.csv')


############################################## MAIN ################################################
# Data pre-processing
# Data pre-processing
data_processed <- data %>%
  
  # remove reassortant classes
  select(-c(ends_with('class'), time_since_last_dominant)) %>%
  
  # scaling
  mutate(across(c(collection_dec, collection_year), .fns = ~ subtract(.x, 2015))) %>%
  
  # set default NA values
  replace_na(list(woah_cases = 0, 
                  woah_susceptibles = 0,
                  woah_deaths = 0,
                  n_sequences = 0,
                  n_reassortants = 0,
                  time_since_last_dominant = 0, # this may be quite a strong assumption (accidentally)
                  minor = 0,
                  major = 0,
                  dominant = 0)) %>%
  
  mutate(across(starts_with('woah'), ~log1p(.x), .names = "{.col}_log1p")) %>%
  mutate(n_sequences = log1p(n_sequences)) %>%
  #rename_with(~gsub('_', '-' ,.x)) %>%
  
  filter(collection_regionname != 'south america')



# Model: hierachical count/detection model. 
# 1) count  -> the 'true' reassortant process, predicted by ecological variables
# 2) detection -> the probability that of true reassortants are observed

# First try with poisson
model_poisson <- brm(
  bf(n_reassortants ~  collection_regionname + woah_susceptibles_log1p +  collection_season ),
  data = data_processed,
  family = poisson(),
  chains = 4, iter = 4000, 
  backend = "cmdstanr", 
  refresh = 0
)


model_negbinom <- brm(
  bf(n_reassortants ~ collection_regionname + woah_susceptibles_log1p +  collection_season ),
  data = data_processed,
  family = negbinomial(),
  chains = 4, 
  iter = 4000, 
  backend = "cmdstanr",
  refresh = 0
)

model_zi_pois <- brm(
  bf(n_reassortants ~  collection_regionname + woah_susceptibles_log1p +  collection_season ),
  data = data_processed,
  family = zero_inflated_poisson(),
  chains = 4, 
  cores = 4,
  iter = 4000, 
  backend = "cmdstanr",
  refresh = 0
)

model_zi_negbinom <- brm(
  bf(n_reassortants ~  collection_regionname + woah_susceptibles_log1p +  collection_season  ),
  data = data_processed,
  family = zero_inflated_negbinomial(),
  chains = 4, 
  cores = 4,
  iter = 4000, 
  backend = "cmdstanr",
  refresh = 0
)

pp_check(model_poisson, ndraws = 100, type = 'bars')
model_performance(model_poisson)

pp_check(model_negbinom, ndraws = 100, type = 'bars')
model_performance(model_negbinom)

pp_check(model_zi_pois, ndraws = 100, type = 'bars')
model_performance(model_zi_pois)

pp_check(model_zi_negbinom, ndraws = 100, type = 'bars')
model_performance(model_zi_pois)


model_zi_pois_full <- brm(
  bf(n_reassortants ~ 0 + collection_regionname + (collection_regionname|collection_year/collection_season)) ,
  data = data_processed,
  family = zero_inflated_poisson(),
  chains = 4, 
  cores = 4,
  iter = 4000, 
  backend = "cmdstanr",
  refresh = 0
)


pp_check(model_zi_pois_full, ndraws = 100, type = 'bars')
model_performance(model_zi_pois_full)

preds <- posterior_predict(model_zi_pois_full, nsamples = 250, summary = FALSE)
preds <- t(preds)

res <- createDHARMa(
  simulatedResponse = t(posterior_predict(model_zi_pois_full)),
  observedResponse = data_processed$n_reassortants,
  fittedPredictedResponse = apply(t(posterior_epred(model_zi_pois_full)), 1, mean),
  integerResponse = FALSE)

plot(res)
###################################################

# Define Priors
mv <- mvbf(bf1, bf2)

countmodel_mv_priors <- get_prior(mv, count_data)

#countmodel_mv_priors[c(1, 4:15,25:29, 38:44, 53:59, 69:74),1] <- "normal(0,5)"
countmodel_mv_priors
countmodel_mv_priors[c(1, 4:15, 25:29, 38:44, 53:59, 69:74),1] <- "normal(0,5)"


# Set MCMC Options
CHAINS <- 4
CORES <- 4
ITER <- 16000
BURNIN <- ITER/10 # Discard 10% burn in from each chain
SEED <- 4472
THREADS <- 2


# Prior Predictive Checks 
# Note that due to weakly informative priors, the prior predictive checks reveal little about the 
# ability of our priors to describe the data. 
countmodel_priorpredictive <- brm(
  mv,
  data = count_data,           
  prior = countmodel_mv_priors,
  chains = CHAINS,
  threads = THREADS, 
  cores = CORES, 
  iter = 16000,
  warmup = 1600,
  seed = SEED,
  backend = "cmdstanr",
  sample_prior = "only",
  control = list(adapt_delta = 0.95)
)


# Fit model to data
countmodel_mv_fit <- brm(
  mv,
  data = count_data,           
  prior = countmodel_mv_priors,
  chains = CHAINS,
  threads = THREADS, 
  cores = CORES, 
  iter = 16000,
  warmup = 1600,
  seed = SEED,
  backend = "cmdstanr",
  control = list(adapt_delta = 0.97))


# Post-fitting diagnostics (including inspection of ESS, Rhat and posterior predictive)
tidy_countmodel_fit <- tidy(countmodel_mv_fit)
#posteriorpredictive <-pp_check(countmodel_fit, ndraws = 500)


# Misc evaluations
#print(prior_summary(countmodel_fit), show_df = FALSE)
performance(countmodel_mv_fit) # tibble output of model metrics including R2, ELPD, LOOIC, RMSE
plot(countmodel_mv_fit) # default output plot of brms showing posterior distributions of
prior_summary(countmodel_mv_fit) #obtain dataframe of priors used in model.

mcmc_countmodel_fit <- ggs(countmodel_mv_fit) # Warning message In custom.sort(D$Parameter) : NAs introduced by coercion

posteriorpredictive_group <-pp_check(countmodel_mv_fit, ndraws = 500, resp = 'reassortantclass')
posteriorpredictive_nreassortants <-pp_check(countmodel_mv_fit, ndraws = 500, resp = 'nreassortants')





############################################## WRITE ###############################################
saveRDS(countmodel_mv_fit, file = './saved_models/count_model.rds')
write_csv(count_data, './saved_models/count_data.csv')


############################################## END #################################################
####################################################################################################
####################################################################################################
# Deprecated

# Marginal probability density of the number of unique reassortants/region
# ie, irrespective of class
count_prob <- jointmodel_temp %>% 
  epred_draws(newdata = count_data %>%
                select(collection_regionname) %>%
                distinct(),resp = "nreassortants",
              re_formula = NA)

count_prob %>% 
  median_hdi(.epred)

ggplot(count_prob, aes(x = .epred, fill = collection_regionname , y = collection_regionname)) +
  stat_halfeye() +
  theme_minimal() + 
  scale_x_continuous('Number of Unique Reassortants/Year') +
  
  # Scales
  scale_y_discrete('Region of Origin', labels = function(x) str_wrap(x, width = 20) %>% str_to_title())+
  scale_colour_manual(
    'Reassortant Class',
    values = region_colour %>% pull(Trait, name = Name))  +
  theme(legend.position = 'none') + 
  coord_cartesian(xlim = c(0,10))

# Marginal probability of reassortant class/region
# ie, irrespective of number
class_prob <- jointmodel_temp %>% 
  epred_draws(newdata = count_data %>%
                select(collection_regionname) %>%
                distinct(),
              resp = "group2",
              re_formula = NA)
class_prob %>% 
  median_hdi(.epred)

ggplot(class_prob, aes(x = collection_regionname, colour = .category, y = .epred)) +
  geom_boxplot() +
  theme_minimal() + 
  scale_y_continuous('Posterior Probability of Reassortant Class', labels = scales::percent) +
  
  # Scales
  scale_x_discrete('Region of Origin', labels = function(x) str_wrap(x, width = 20) %>% str_to_title())+
  scale_colour_manual(
    'Reassortant Class',
    values = riskgroup_colour %>% pull(Trait, name = group2)) +
  theme(legend.position = 'bottom')


#### Contrasts ####
fit_year %>% 
  emmeans(., ~ collection_regionname ,
          epred = TRUE,
          resp = 'group2') %>% 
  filter(rep.meas == 'nreassortants') %>%
  contrast(method = 'revpairwise')



joint_prob <- count_prob %>%
  left_join(class_prob, 
            by = c(".draw",  '.row', "collection_regionname"), 
            suffix = c("_count", "_class")) %>%
  mutate(joint_epred = .epred_count * .epred_class)