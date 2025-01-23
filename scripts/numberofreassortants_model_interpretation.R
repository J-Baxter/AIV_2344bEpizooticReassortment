####################################################################################################
####################################################################################################
## Script name: Count and Class Model Interpretation and Figures
##
## Purpose of script: To extract, interpret and display the parameter values and predictions of the
## hurdle model fitted in ./scripts/numberofreassortants_model.R
##
## Date created: 2025-01-22
##
##
########################################## SYSTEM OPTIONS ##########################################
options(scipen = 6, digits = 7) 
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
scientific_10 <- function(x) {
  parse(text=gsub("e\\+*", " %*% 10^", scales::scientific_format()(x)))
}



############################################## DATA ################################################
countmodel_mv_fit <- readRDS('./saved_models/count_model.rds')
count_data <- read_csv('./saved_models/count_data.csv')


# MCMC chains
mcmc_countmodel <- ggs(countmodel_mv_fit) # Warning message In custom.sort(D$Parameter) : NAs introduced by coercion


# Posterior predictive distribution
posteriorpredictive_reassortantclass <-pp_check(countmodel_fit, ndraws = 500, resp = 'reassortantclass')
posteriorpredictive_nreassortants <-pp_check(countmodel_fit, ndraws = 500, resp = 'nreassortants')

                                             

############################################## MAIN ################################################




# Conditional prediction draws  + Conditional marginal means stratified by continent

# Required outputs: Pre/post epizootic, marginal effect of region ,
# 


# Posterior prediction of the number and type of reassortant at 'average' parameters (global and stratified by region)


# marginal effect of of continent on number and type of reassortant (ie how many more reassortants in a vs b)


# reassortant 'transition' probabilities (ie based on lagged reassortant type), global and stratified by region


# reassortant  probabilities according to time since last dominant (stratified by region)
# Random Effects:

# Year (does the probability of number and type change over each year per continent)




# Marginal effect of region


# Marginal class probability  


# Marginal class probabilty ~ 1| Region


# Joint Probability



# Joint Probability ~ 1|Region


# Marginal effect of region



##### Plot Posterior Predictive Check #####
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
  
  guides(alpha= 'none', 
         colour=guide_legend()) + 
  scale_x_continuous('Weighted Diffusion Coefficient') + 
  scale_y_continuous('Density') + 
  
  theme_minimal(base_size = 8) + 
  theme(legend.position = 'inside',
        legend.title = element_blank(),
        legend.position.inside = c(0.8, 0.6))



#### Plot MCMC chains for Key parameters #####
plt_mcmc <- mcmc_diffusion %>% 
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


#### Plot Prior & posterior distributions of key parameters ####
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


############################################## WRITE ###############################################




############################################## END #################################################
####################################################################################################
####################################################################################################