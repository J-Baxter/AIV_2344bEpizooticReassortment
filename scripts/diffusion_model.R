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
  mutate(across(c('median_anseriformes_wild', 
                  'median_charadriiformes_wild', 
                  'count_cross_species'), .fns= ~ log1p(.x), .names = '{.col}_log1p'))  %>%
  
  dplyr::select(weighted_diff_coeff, 
                median_anseriformes_wild_log1p,
                median_charadriiformes_wild_log1p,
                count_cross_species_log1p,
                host_simplifiedhost,
                segment, 
                collection_regionname) %>%
  
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
#diffusion_formula <- bf(weighted_diff_coeff ~ 0 +collection_regionname + median_anseriformes_wild + median_charadriiformes_wild +  collection_season + (1|segment),
                     #   shape ~ 0 +  collection_regionname + (1 | segment))
#diffusion_formula <- bf(weighted_diff_coeff ~ 0 + median_anseriformes_wild +  median_charadriiformes_wild + collection_regionname + collection_season +  (1|segment))

int_step <- function(x){
  ifelse(x >0, 1,0)
}
diffusion_formula <- bf(weighted_diff_coeff ~ 0 + collection_regionname + int_step(median_anseriformes_wild_log1p) + median_anseriformes_wild_log1p + int_step(median_charadriiformes_wild_log1p) + median_charadriiformes_wild_log1p + int_step(count_cross_species_log1p) + count_cross_species_log1p + 
                          collection_regionname:median_anseriformes_wild_log1p + 
                          collection_regionname:median_charadriiformes_wild_log1p + (1|segment) ,
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
performance(diffusionmodel1_fit_gamma_14)


# Misc evaluations
# performance(diffusionmodel1_fit) # tibble output of model metrics including R2, ELPD, LOOIC, RMSE
plot(diffusionmodel1_fit) # default output plot of brms showing posterior distributions of


#### Check Residuals using DHARMA ####
# sample from the Posterior Predictive Distribution
preds <- posterior_predict(diffusionmodel1_fit_gamma_19, nsamples = 250, summary = FALSE)
preds <- t(preds)

res <- createDHARMa(
  simulatedResponse = t(posterior_predict(diffusionmodel1_fit_gamma_18)),
  observedResponse = diffusion_data$weighted_diff_coeff,
  fittedPredictedResponse = apply(t(posterior_epred(diffusionmodel1_fit_gamma_18)), 1, mean),
  integerResponse = FALSE)

plot(res, quantreg = FALSE)


#### Check Ratio of Effective Population Size to Total Sample Size #### 
# values <0.1 should raise concerns about autocorrelation

neff_ratio(diffusionmodel1_fit_gamma_19) %>% as_tibble(rownames = 'param') %>%
  mutate(param = fct_reorder(param, desc(value))) %>%
  ggplot() + 
  geom_segment(aes(yend = value,
                   xend=param, 
                   y=0,
                   x = param,
                   colour = value > 0.1)) +
  geom_point(aes(y = value,
                 x = param,
                 colour = value > 0.1)) + 
  geom_hline(aes(yintercept = 0.1), linetype = 'dashed') + 
  geom_hline(aes(yintercept = 0.5), linetype = 'dashed') + 
  scale_colour_manual(values = c( '#0047AB',  'red')) + 
  scale_y_continuous(limits = c(0, 1), expand = c(0,0),
                     expression(N["eff"]/N)) + 
  scale_x_discrete(expand= c(0.1,0), 'Fitted Parameter') + 
  theme_classic() + 
  coord_flip()  + 
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = 'none') 
 
  
trank <- as_draws_df(diffusionmodel1_fit_gamma_19)  %>%  
  mcmc_rank_overlay(regex_pars = '^b_') 

trank$data %>%
  rename(.variable = parameter) %>%
  mutate(label = case_when(.variable == 'b_collection_regionnameasia'~ "beta['asia']",
                           .variable == 'b_collection_regionnameafrica'~ "beta['africa']",
                           .variable == 'b_collection_regionnameeurope'~ "beta['europe']",
                           .variable == "b_collection_regionnamecentral&northernamerica"~ "beta['americas']",
                           
                           .variable == 'b_int_stepmedian_anseriformes_wild_log1p'~ "beta['step_anseriformes']",
                           .variable == 'b_median_anseriformes_wild_log1p'~ "beta['persist_anseriformes']",
                           .variable == 'b_int_stepmedian_charadriiformes_wild_log1p'~ "beta['step_charadriiformes']",
                           .variable == 'b_median_charadriiformes_wild_log1p'~ "beta['persist_charadriiformes']",
                           .variable == 'b_int_stepcount_cross_species_log1p'~ "beta['step_hostjump']",
                           .variable == 'b_count_cross_species_log1p'~ "beta['num_hostjump']",
                           
                           .variable == 'b_collection_regionnameafrica:median_anseriformes_wild_log1p'~ "beta['africa-anser']",
                           .variable == 'b_collection_regionnameeurope:median_anseriformes_wild_log1p'~ "beta['europe-anser']",
                           .variable == 'b_collection_regionnamecentral&northernamerica:median_anseriformes_wild_log1p'~ "beta['americas-anser']",
                           
                           .variable == 'b_collection_regionnameafrica:median_charadriiformes_wild_log1p'~ "beta['africa-charad']",
                           .variable == 'b_collection_regionnameeurope:median_charadriiformes_wild_log1p'~ "beta['europe-charad']",
                           .variable == 'b_collection_regionnamecentral&northernamerica:median_charadriiformes_wild_log1p'~ "beta['americas-charad']",
                           
                           .variable == 'b_shape_collection_regionnameasia'~ "alpha['asia']",
                           .variable == 'b_shape_collection_regionnameafrica'~ "alpha['africa']",
                           .variable == 'b_shape_collection_regionnameeurope'~ "alpha['europe']",
                           .variable == 'b_shape_collection_regionnamecentral&northernamerica'~ "alpha['americas']")) %>%
  ggplot() + 
  geom_step(aes(x = bin_start, 
                y = n,
                colour = chain),
            linewidth = 0.8) + 
  scale_y_continuous(expand = c(0,0), limits = c(150, 225)) +
  scale_x_continuous(expand = c(0,0)) + 
  facet_wrap(~label,  ncol = 3, labeller = label_parsed) +
  scale_colour_brewer(palette = 'GnBu', 'Chains') +
  theme_minimal() + 
  theme(legend.position = 'bottom',
        axis.title = element_blank())



saveRDS(diffusionmodel1_fit, './saved_models/diffusion_model_2.rds')
write_csv(diffusion_data, './saved_models/diffusion_model.csv')
############################################## END #################################################
####################################################################################################
####################################################################################################