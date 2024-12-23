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
  
  rename_with(~gsub('-', '_', .x)) %>%
  mutate(collection_regionname = case_when(grepl('europe', collection_regionname) ~ 'europe',
                                           grepl('africa', collection_regionname) ~ 'africa',
                                           grepl('asia', collection_regionname) ~ 'asia',
                                           grepl('(central|northern) america', collection_regionname) ~ 'central & northern america',
                                           .default = collection_regionname
  )) %>%
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
         ))



# Model Formula
# We assume a hurdle lognormal model, in which the hurdle process is determined the season in which 
# the TMRCA is estimated to be, and the lognormal component is determined by persistence
# in anseriformes and charadriiformes. Both model components are conditional on segment from which 
# the measurement is taken and region of origin

diffusion_formula <- bf(weighted_diff_coeff ~ 1 + median_anseriformes_wild +  median_charadriiformes_wild + 
                          (1|segment+ collection_regionname),
                        hu ~ 1 +  season +(1|segment + collection_regionname))


# Define Priors
diffusionmodel1_priors <- get_prior(diffusion_formula,
                                    data = diffusion_data,
                                    family = hurdle_lognormal()) 

diffusionmodel1_priors$prior[c(1,11)] <- "normal(0,5)"

diffusionmodel1_priors$prior[9] <- "student_t(3, 0, 3)"
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


# Posterior Predictions and Marginal Effects 
# HU part of model

## 1. Global Average
global_average <- diffusionmodel1_fit %>%
  epred_draws(re_formula = NA,
              newdata =  .$data)


## 2. predicted probability of diffusion coefficient == 0 | season
# need to change to draws for density plot
predict_hu_season <- diffusionmodel1_fit %>%
  emmeans(~ season,
          var = "weighted_diff_coeff",
          at = list(continent = unique(diffusion_data$season)),
          #epred = TRUE, 
          dpar = "hu",
          re_formula = NA, regrid = "response",
          allow_new_levels = TRUE)

# Print estimates as 1-x, to infer P(X>=1)
predict_hu_season %>%
  as_tibble() %>% 
  mutate(across(where(is.numeric), .fns = ~ 1 -.x))

# Pairwise Contrast
diffusionmodel1_fit %>%
  emmeans(~ 1 + season,
          var = "weighted_diff_coeff",
          at = list(continent = unique(diffusion_data$season)),
          # epred = TRUE, 
          dpar = "hu",
          re_formula = NA ,
          allow_new_levels = TRUE) %>%
  contrast(method = "pairwise",
           type = 'response')


## 3. Prediction for region groups 
regional_average <- diffusionmodel1_fit %>%
  emmeans(~ 1 + collection_regionname,
          var = "weighted_diff_coeff",
          at = list(continent = unique(diffusion_data$collection_regionname)),
          # epred = TRUE, 
          dpar = "hu",
          re_formula = ~ 1|collection_regionname , regrid = "response",
          allow_new_levels = TRUE)

# Print estimates as 1-x, to infer P(X>=1)
regional_average %>% as_tibble() %>% mutate(across(where(is.numeric), .fns = ~ 1 -.x))

# Pairwise Contrast
diffusionmodel1_fit %>%
  emmeans(~ 1 + collection_regionname,
          var = "weighted_diff_coeff",
          at = list(continent = unique(diffusion_data$collection_regionname)),
          # epred = TRUE, 
          dpar = "hu",
          re_formula = ~ 1|collection_regionname ,
          allow_new_levels = TRUE) %>%
  contrast(method = "pairwise",
           type = 'response')



## MU part of model
# MU ~  anseriformes wild + (1|region)
anser_byregion_mu <- diffusionmodel1_fit %>%
  emmeans(~ 1 + median_anseriformes_wild + collection_regionname,
          var = "weighted_diff_coeff",
          at = list(continent = unique(diffusion_data$collection_regionname),
                    median_anseriformes_wild = seq(0.5, 10, by = 0.5)),
          #epred = TRUE,
          dpar = "mu",
          re_formula = ~ 1|collection_regionname ,  
          regrid = "response",
          tran = "log", 
          type = "response",
          allow_new_levels = TRUE)


# MU ~  charadriiformes wild + (1|region)
charad_byregion_mu <- diffusionmodel1_fit %>%
  emmeans(~ 1 + median_charadriiformes_wild + collection_regionname,
          var = "weighted_diff_coeff",
          at = list(continent = unique(diffusion_data$collection_regionname),
                    median_charadriiformes_wild = seq(0.5, 10, by = 0.5)),
          #epred = TRUE,
          dpar = "mu",
          re_formula = ~ 1|collection_regionname ,  
          regrid = "response",
          tran = "log", 
          type = "response",
          allow_new_levels = TRUE)


# MU ~ 1|Region Predictions
diffusionmodel1_fit %>%
  emmeans(~ 1 + collection_regionname,
          var = "weighted_diff_coeff",
          at = list(continent = unique(diffusion_data$collection_regionname)),
          # epred = TRUE, 
          dpar = "mu",
          re_formula = ~ 1|collection_regionname ,  regrid = "response",
          tran = "log", 
          type = "response",
          allow_new_levels = TRUE)


# Marginal effect of charadriiformes on MU
diffusionmodel1_fit %>%
  emtrends(~ 1 + median_charadriiformes_wild,
           var = "median_charadriiformes_wild",
           at = list(median_charadriiformes_wild = seq(0, 6, by = 0.5)),
           # epred = TRUE, 
           dpar = "mu", regrid= "response") 


# Combined mixture model (ie HU + MU)



############################################## WRITE ###############################################
# print(prior_summary(diffusionmodel1_fit), show_df = FALSE)

# Plot Diffusion Coefficient Distribution
plt_diffusiondata <- diffusion_data %>%
  ggplot() +
  geom_histogram(aes(x = log1p(weighted_diff_coeff), 
                     fill = weighted_diff_coeff>0)) +
  scale_fill_brewer(palette = 'Dark2', 'Is Zero') +
  scale_x_continuous('Log Weighted Diffusion Coefficient') +
  scale_y_continuous('Count') + 
  theme_minimal(base_size = 8) + 
  theme(legend.position = 'inside',
        legend.position.inside = c(0.8,0.8))


# Plot Posterior Predictive Check
plt_posteriorpredictive <- plot_posteriorpredictive$data  %>%
  mutate(value = log1p(value)) %>%
  ggplot() + 
  geom_density(aes(x = value, 
                   group= rep_id, 
                   alpha = is_y,  
                   linewidth = is_y, 
                   colour = is_y), 
               key_glyph = draw_key_path) + 
  scale_alpha_manual(values = c('FALSE' = 0.00001, 
                                'TRUE' = 1)) + 
  scale_linewidth_manual(values = c('FALSE' = 0.1, 
                                    'TRUE' = 1),
                         labels  = c(expression(italic('y')['rep']),
                                     expression(italic('y'))))+ 
  scale_colour_manual(values = c('FALSE' = '#cbc9e2', 
                                 'TRUE' = '#54278f'),
                      labels  = c(expression(italic('y')['rep']),
                                  expression(italic('y')))) + 
  
  l
guides(alpha= 'none', 
       colour=guide_legend()) + 
  scale_x_continuous('Weighted Diffusion Coefficient') + 
  scale_y_continuous('Density') + 
  
  theme_minimal(base_size = 8) + 
  theme(legend.position = 'inside',
        legend.title = element_blank(),
        legend.position.inside = c(0.8, 0.6))



# Plot MCMC chains for Key parameters
plt_mcmc <-mcmc_diffusion %>% 
  filter(Iteration > 400) %>%
  mutate(Parameter = factor(Parameter,
                            levels = c("b_Intercept", "b_hu_Intercept", "b_median_anseriformes_wild",
                                       "b_median_charadriiformes_wild", "b_hu_seasonmigrating_spring", "b_hu_seasonmigrating_autumn",
                                       "b_hu_seasonoverwintering", "sd_collection_regionname__Intercept",
                                       "sd_segment__Intercept", "sd_collection_regionname__hu_Intercept", 
                                       "sd_segment__hu_Intercept", "sigma"),
                            labels = c(expression(paste(beta[0])),
                                       expression(paste(beta['hu'])),
                                       expression(paste(beta['anseriformes'])),
                                       expression(paste(beta['charadriiformes'])),
                                       expression(paste(beta['hu_seasonmigrating_spring'])),
                                       expression(paste(beta['hu_seasonmigrating_autumn'])),
                                       expression(paste(beta['hu_seasonoverwinterin'])),
                                       expression(paste(sigma['collection_regionname'])),
                                       expression(paste(sigma[mu])),
                                       expression(paste(sigma['hu_collection_regionname'])),
                                       expression(paste(sigma['hu_segment'])),
                                       expression(paste(sigma[epsilon]))
                            )
  ))  %>% drop_na(Parameter) %>%
  ggplot(aes(x = Iteration,
             y = value, 
             col = as.factor(Chain)))+
  geom_line(alpha = 0.8)+
  facet_wrap(~ Parameter,
             ncol = 2,
             scale  = 'free_y',
             switch = 'y',
             labeller = label_parsed)+
  scale_colour_brewer(palette = 'GnBu', 'Chains') +
  theme_minimal(base_size = 8) + 
  theme(legend.position = 'bottom')


# Plot Prior & posterior distributions of key parameters
posterior_beta_draws <- as.data.frame(diffusionmodel1_fit) %>%
  as_tibble() %>%
  pivot_longer(., cols = everything(),
               names_to = 'parameter',
               values_to = 'estimate') %>% filter(!grepl('lp__|lprior', parameter)) %>%
  #filter(grepl('b_', parameter)) %>%
  mutate(draw = 'posterior')


priors <- expand_grid(parameter = unique(posterior_beta_draws$parameter),
                      mu = NA_real_,
                      sigma = NA_real_,
                      df = NA_real_) %>%
  mutate(mean = case_when(grepl('^b', parameter) ~ 0,
                          grepl('^sd_(segment|collection)|^sigma$', parameter) ~ 0,
                          grepl('Intercept_hu', parameter) ~ 0),
         sd = case_when(grepl('^b', parameter) ~ 10,
                        grepl('^sd_(segment|collection)|^sigma$', parameter) ~ 11.7,
                        grepl('Intercept_hu', parameter) ~ 1),
         df = case_when(grepl('^sd_(segment|collection)|^sigma$', parameter) ~ 3),
         dist = case_when(grepl('^b', parameter) ~ 'norm',
                          grepl('^sd_(segment|collection)|^sigma$', parameter) ~ 'student_t',
                          grepl('Intercept_hu', parameter) ~ 'logistic')) %>%
  drop_na(mean) 


plt_params <- ggplot() + 
  geom_histogram(data = posterior_beta_draws %>% 
                   filter(parameter %in% priors$parameter), 
                 aes(x = estimate,
                     y = after_stat(density)),
                 inherit.aes = F, 
                 bins = 50, 
                 fill = '#1b9e77') + 
  
  stat_function(fun = dnorm,
                data = tibble(parameter = "b_Intercept" ),
                args = list(mean = 0, sd = 3),
                fill = '#d95f02',
                geom = 'area',
                alpha = 0.5) +
  
  stat_function(fun = dnorm,
                data = tibble(parameter = "b_hu_Intercept"),
                args = list(mean = 0, sd = 3),
                fill = '#d95f02',
                geom = 'area',
                alpha = 0.5) +
  
  stat_function(fun = dnorm,
                data = tibble(parameter = "b_host_richness"),
                args = list(mean = 0, sd = 3),
                fill = '#d95f02',
                geom = 'area',
                alpha = 0.5) +
  
  stat_function(fun = dnorm,
                data = tibble(parameter = "b_median_anseriformes_wild"),
                args = list(mean = 0, sd = 3),
                fill = '#d95f02',
                geom = 'area',
                alpha = 0.5) +
  
  stat_function(fun = dnorm,
                data = tibble(parameter = "b_median_charadriiformes_wild"),
                args = list(mean = 0, sd = 3),
                fill = '#d95f02',
                geom = 'area',
                alpha = 0.5) +
  
  stat_function(fun = dnorm,
                data = tibble(parameter = "b_seasonmigrating_spring"),
                args = list(mean = 0, sd = 3),
                fill = '#d95f02',
                geom = 'area',
                alpha = 0.5) +
  
  stat_function(fun = dnorm,
                data = tibble(parameter = "b_seasonmigrating_autumn"),
                args = list(mean = 0, sd = 3),
                fill = '#d95f02',
                geom = 'area',
                alpha = 0.5) +
  
  stat_function(fun = dnorm,
                data = tibble(parameter = "b_seasonoverwintering"),
                args = list(mean = 0, sd = 3),
                fill = '#d95f02',
                geom = 'area',
                alpha = 0.5) +
  
  stat_function(fun = dnorm,
                data = tibble(parameter = "b_hu_seasonmigrating_spring"),
                args = list(mean = 0, sd = 3),
                fill = '#d95f02',
                geom = 'area',
                alpha = 0.5) +
  
  stat_function(fun = dnorm,
                data = tibble(parameter = "b_hu_seasonmigrating_autumn"),
                args = list(mean = 0, sd = 3),
                fill = '#d95f02',
                geom = 'area',
                alpha = 0.5) +
  
  stat_function(fun = dnorm,
                data = tibble(parameter = "b_hu_seasonoverwintering"),
                args = list(mean = 0, sd = 3),
                fill = '#d95f02',
                geom = 'area',
                alpha = 0.5) +
  
  stat_function(fun = dstudent_t,
                data = tibble(parameter = "sd_collection_regionname__Intercept"),
                args = list(mu = 0, sigma = 5, df = 3),
                fill = '#d95f02',
                geom = 'area',
                alpha = 0.5) +
  
  stat_function(fun = dstudent_t,
                data = tibble(parameter = "sd_collection_regionname__hu_Intercept"),
                args = list(mu = 0, sigma = 5, df = 3),
                fill = '#d95f02',
                geom = 'area',
                alpha = 0.5) +
  
  stat_function(fun = dstudent_t,
                data = tibble(parameter = "sd_segment__Intercept"),
                args = list(mu = 0, sigma = 5, df = 3),
                fill = '#d95f02',
                geom = 'area',
                alpha = 0.5) +
  
  stat_function(fun = dstudent_t,
                data = tibble(parameter = "sd_segment__hu_Intercept"),
                args = list(mu = 0, sigma = 5, df = 3),
                fill = '#d95f02',
                geom = 'area',
                alpha = 0.5) +
  
  stat_function(fun = dstudent_t,
                data = tibble(parameter = "sigma"),
                args = list(mu = 0, sigma = 5, df = 3),
                fill = '#d95f02',
                geom = 'area',
                alpha = 0.5) +
  
  xlim(c(-15,15)) + 
  facet_wrap(~parameter, scales = 'free_y') +
  theme_minimal()


# Posterior Prediction HU ~ Season
plt_hu_season <- predict_hu_season %>%
  gather_emmeans_draws() %>%
  ggplot(aes(x  = 1-.value,
             fill = season)) +
  stat_halfeye(alpha = 0.5) +
  scale_fill_brewer(palette = 'Set1') +
  scale_x_continuous('p(Weighted Diffusion Coefficient > 0)') + 
  scale_y_continuous('Probability Density') + 
  theme_minimal(base_size = 8) + 
  theme(legend.position = 'inside',
        legend.position.inside = c(0.8,0.8))


# Posterior Prediction HU ~ 1|Region
plt_hu_region <- regional_average %>%
  gather_emmeans_draws() %>%
  ggplot(aes(x  = 1-.value,
             y = collection_regionname, 
             fill = collection_regionname)) +
  stat_halfeye() +
  scale_fill_brewer(palette = 'Dark2') +
  scale_x_continuous('p(Weighted Diffusion Coefficient > 0)') + 
  scale_y_discrete('Region', 
                   labels = function(x) str_wrap(x, width = 10)) + 
  theme_minimal(base_size = 8) + 
  theme(legend.position = 'none')


# Posterior Prediction MU ~ anseriformes + (1 | Region)
plt_mu_anseriformes <- anser_byregion_mu %>%
  gather_emmeans_draws() %>%
  median_hdi(., .width =  .95) %>%
  mutate(across(c(.value, .lower, .upper),
                .fns = ~ exp(.x))) %>%
  ggplot() +
  geom_line(aes(x  = median_anseriformes_wild, 
                y = .value,
                colour = collection_regionname)) + 
  geom_ribbon(aes(x  = median_anseriformes_wild, 
                  ymin = .lower, 
                  ymax = .upper, 
                  fill = collection_regionname), 
              alpha = 0.05) + 
  #ggplot(aes(x  = median_anseriformes_wild, y = .value, colour = collection_regionname, fill = collection_regionname)) +
  #stat_lineribbon(.width =  0.95) +
  scale_colour_brewer(palette = 'Dark2',) +
  scale_fill_brewer(palette = 'Dark2') +
  scale_x_continuous('Persistence Time in Anseriformes (Wild)') + 
  scale_y_continuous('Weighted Diffusion Coefficient') + 
  theme_minimal(base_size = 8) + 
  theme(legend.position = 'none')


# Posterior Prediction MU ~ charadriiformes + (1 | Region)
plt_mu_charadriiformes <- charad_byregion_mu %>%
  gather_emmeans_draws() %>%
  group_by(median_charadriiformes_wild, collection_regionname) %>%
  median_hdi(., .width =  .95) %>%
  mutate(across(c(.value, .lower, .upper), .fns = ~ exp(.x))) %>%
  ggplot() +
  geom_line(aes(x  = median_charadriiformes_wild,
                y = .value, 
                colour = collection_regionname)) + 
  geom_ribbon(aes(x  = median_charadriiformes_wild, 
                  ymin = .lower, ymax = .upper,
                  fill = collection_regionname), alpha = 0.05) + 
  #stat_lineribbon(.width =  0.95) +
  scale_colour_brewer(palette = 'Dark2') +
  scale_fill_brewer(palette = 'Dark2') +
  scale_x_continuous('Persistence Time in Charadriiformes (Wild)') + 
  scale_y_continuous('Weighted Diffusion Coefficient') + 
  theme_minimal(base_size = 8) + 
  theme(legend.position = 'none')


############################################## END #################################################
####################################################################################################
####################################################################################################