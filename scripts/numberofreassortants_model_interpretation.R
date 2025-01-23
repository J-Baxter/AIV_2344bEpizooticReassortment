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
posteriorpredictive_reassortantclass <-pp_check(countmodel_mv_fit, ndraws = 500, resp = 'reassortantclass')
posteriorpredictive_nreassortants <-pp_check(countmodel_mv_fit, ndraws = 500, resp = 'nreassortants')

countmodel_posteriorpredictive <- bind_rows(posteriorpredictive_nreassortants$data %>%
  mutate(resp = 'nreassortants'),
  posteriorpredictive_reassortantclass$data %>%
    mutate(resp = 'reassortantclass'))
  

                                             

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
pp_a <- countmodel_posteriorpredictive  %>%
  filter(resp == 'nreassortants') %>%
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
  scale_x_continuous('Number of Reassortants',
                     expand = c(0,0),
  ) +
  scale_y_continuous('Density',expand = c(0,0)) + 
  global_theme+ 
  theme(legend.position = 'inside',
        legend.title = element_blank(),
        legend.position.inside = c(0.8, 0.6))


pp_b <- countmodel_posteriorpredictive  %>%
  filter(resp == 'reassortantclass') %>%
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
  scale_x_continuous('Reassortant Class',
                     expand = c(0,0),
  ) +
  scale_y_continuous('Density',expand = c(0,0)) + 
  global_theme+ 
  theme(legend.position = 'inside',
        legend.title = element_blank(),
        legend.position.inside = c(0.8, 0.6))

cowplot::plot_grid(pp_a,  pp_b, nrow=1,align='h',axis='tb',labels='AUTO')


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
posterior_beta_draws <- as.data.frame(countmodel_mv_fit) %>%
  as_tibble() %>%
  pivot_longer(., cols = everything(),
               names_to = 'parameter',
               values_to = 'estimate') %>% 
  filter(!grepl('lp__|lprior|^r_|^cor|^sd', parameter)) %>%
  #filter(grepl('b_', parameter)) %>%
  mutate(draw = 'posterior') 


plt_params <- ggplot() + 
  geom_histogram(data = posterior_beta_draws , 
                 aes(x = estimate,
                     y = after_stat(density)),
                 inherit.aes = F, 
                 binwidth = 0.1, 
                 fill = '#1b9e77') + 
  
  # Intercepts
  stat_function(fun = dlogis,
                data = tibble(parameter = "Intercept_zi_nreassortants" ),
                args = list(location = 0, scale = 1),
                fill = '#d95f02',
                geom = 'area',
                alpha = 0.5) +
  
  stat_function(fun = dstudent_t,
                data = tibble(parameter = "Intercept_nreassortants"),
                args = list(mu = 3, sigma = -2.3, df = 2.5),
                fill = '#d95f02',
                geom = 'area',
                alpha = 0.5) +
  
  stat_function(fun = dstudent_t,
                data = tibble(parameter = "Intercept_mumajor_reassortantclass"),
                args = list(mu = 3, sigma = 0, df = 2.5),
                fill = '#d95f02',
                geom = 'area',
                alpha = 0.5) +
  
  stat_function(fun = dstudent_t,
                data = tibble(parameter = "Intercept_muminor_reassortantclass"),
                args = list(mu = 3, sigma = 0, df = 2.5),
                fill = '#d95f02',
                geom = 'area',
                alpha = 0.5) +
  
  stat_function(fun = dstudent_t,
                data = tibble(parameter = "Intercept_munone_reassortantclass"),
                args = list(mu = 3, sigma = 0, df = 2.5),
                fill = '#d95f02',
                geom = 'area',
                alpha = 0.5) +
  
  # Beta - number of reassortants
  stat_function(fun = dnorm,
                data = tibble(parameter = "b_nreassortants_Intercept" ),
                args = list(mean = 0, sd = 5),
                fill = '#d95f02',
                geom = 'area',
                alpha = 0.5) +
  
  stat_function(fun = dnorm,
                data = tibble(parameter = "b_nreassortants_collection_regionnameasia" ),
                args = list(mean = 0, sd = 5),
                fill = '#d95f02',
                geom = 'area',
                alpha = 0.5) +
  
  stat_function(fun = dnorm,
                data = tibble(parameter = "b_nreassortants_collection_regionnamecentral&northernamerica" ),
                args = list(mean = 0, sd = 5),
                fill = '#d95f02',
                geom = 'area',
                alpha = 0.5) +
  
  stat_function(fun = dnorm,
                data = tibble(parameter = "b_nreassortants_collection_regionnameeurope" ),
                args = list(mean = 0, sd = 5),
                fill = '#d95f02',
                geom = 'area',
                alpha = 0.5) +
  
  stat_function(fun = dnorm,
                data = tibble(parameter = "b_nreassortants_n_cases" ),
                args = list(mean = 0, sd = 5),
                fill = '#d95f02',
                geom = 'area',
                alpha = 0.5) +
  
  stat_function(fun = dnorm,
                data = tibble(parameter = "b_nreassortants_previous_reassortant_classmajor" ),
                args = list(mean = 0, sd = 5),
                fill = '#d95f02',
                geom = 'area',
                alpha = 0.5) +
  
  stat_function(fun = dnorm,
                data = tibble(parameter = "b_nreassortants_previous_reassortant_classminor" ),
                args = list(mean = 0, sd = 5),
                fill = '#d95f02',
                geom = 'area',
                alpha = 0.5) +
  
  stat_function(fun = dnorm,
                data = tibble(parameter = "b_nreassortants_previous_reassortant_classnone" ),
                args = list(mean = 0, sd = 5),
                fill = '#d95f02',
                geom = 'area',
                alpha = 0.5) +
  
  stat_function(fun = dnorm,
                data = tibble(parameter = "b_nreassortants_collection_seasonmigrating_autumn" ),
                args = list(mean = 0, sd = 5),
                fill = '#d95f02',
                geom = 'area',
                alpha = 0.5) +
  
  stat_function(fun = dnorm,
                data = tibble(parameter = "b_nreassortants_collection_seasonmigrating_spring" ),
                args = list(mean = 0, sd = 5),
                fill = '#d95f02',
                geom = 'area',
                alpha = 0.5) +
  
  stat_function(fun = dnorm,
                data = tibble(parameter = "b_nreassortants_collection_seasonoverwintering" ),
                args = list(mean = 0, sd = 5),
                fill = '#d95f02',
                geom = 'area',
                alpha = 0.5) +
  
  stat_function(fun = dnorm,
                data = tibble(parameter = "b_nreassortants_previous_reassortant_classnone" ),
                args = list(mean = 0, sd = 5),
                fill = '#d95f02',
                geom = 'area',
                alpha = 0.5) +
  
  
  # Beta - ZI nreassortants
  stat_function(fun = dnorm,
                data = tibble(parameter = "b_zi_nreassortants_Intercept" ),
                args = list(mean = 0, sd = 5),
                fill = '#d95f02',
                geom = 'area',
                alpha = 0.5) +
  
  stat_function(fun = dnorm,
                data = tibble(parameter = "b_zi_nreassortants_collection_regionnameasia" ),
                args = list(mean = 0, sd = 5),
                fill = '#d95f02',
                geom = 'area',
                alpha = 0.5) +
  
  stat_function(fun = dnorm,
                data = tibble(parameter = "b_zi_nreassortants_collection_regionnamecentral&northernamerica" ),
                args = list(mean = 0, sd = 5),
                fill = '#d95f02',
                geom = 'area',
                alpha = 0.5) +
  
  stat_function(fun = dnorm,
                data = tibble(parameter = "b_zi_nreassortants_collection_regionnameeurope" ),
                args = list(mean = 0, sd = 5),
                fill = '#d95f02',
                geom = 'area',
                alpha = 0.5) +
  
  stat_function(fun = dnorm,
                data = tibble(parameter = "b_zi_nreassortants_n_cases" ),
                args = list(mean = 0, sd = 5),
                fill = '#d95f02',
                geom = 'area',
                alpha = 0.5) +
  
  # Beta - reassortant class
  
  stat_function(fun = dnorm,
                data = tibble(parameter = "b_mumajor_reassortantclass_Intercept" ),
                args = list(mean = 0, sd = 5),
                fill = '#d95f02',
                geom = 'area',
                alpha = 0.5) +
  
  stat_function(fun = dnorm,
                data = tibble(parameter = "b_muminor_reassortantclass_Intercept" ),
                args = list(mean = 0, sd = 5),
                fill = '#d95f02',
                geom = 'area',
                alpha = 0.5) +
  
  stat_function(fun = dnorm,
                data = tibble(parameter = "b_munone_reassortantclass_Intercept" ),
                args = list(mean = 0, sd = 5),
                fill = '#d95f02',
                geom = 'area',
                alpha = 0.5) +
  
  stat_function(fun = dnorm,
                data = tibble(parameter = "b_mumajor_reassortantclass_collection_regionnameasia" ),
                args = list(mean = 0, sd = 5),
                fill = '#d95f02',
                geom = 'area',
                alpha = 0.5) +
  
  stat_function(fun = dnorm,
                data = tibble(parameter = "b_mumajor_reassortantclass_collection_regionnamecentral&northernamerica" ),
                args = list(mean = 0, sd = 5),
                fill = '#d95f02',
                geom = 'area',
                alpha = 0.5) +
  
  stat_function(fun = dnorm,
                data = tibble(parameter = "b_mumajor_reassortantclass_collection_regionnameeurope" ),
                args = list(mean = 0, sd = 5),
                fill = '#d95f02',
                geom = 'area',
                alpha = 0.5) +
  
  stat_function(fun = dnorm,
                data = tibble(parameter = "b_mumajor_reassortantclass_previous_reassortant_classmajor"),
                args = list(mean = 0, sd = 5),
                fill = '#d95f02',
                geom = 'area',
                alpha = 0.5) +
  
  stat_function(fun = dnorm,
                data = tibble(parameter = "b_mumajor_reassortantclass_previous_reassortant_classminor"),
                args = list(mean = 0, sd = 5),
                fill = '#d95f02',
                geom = 'area',
                alpha = 0.5) +
  
  stat_function(fun = dnorm,
                data = tibble(parameter = "b_mumajor_reassortantclass_previous_reassortant_classnone"),
                args = list(mean = 0, sd = 5),
                fill = '#d95f02',
                geom = 'area',
                alpha = 0.5) +
  
  stat_function(fun = dnorm,
                data = tibble(parameter = "b_muminor_reassortantclass_collection_regionnameasia" ),
                args = list(mean = 0, sd = 5),
                fill = '#d95f02',
                geom = 'area',
                alpha = 0.5) +
  
  stat_function(fun = dnorm,
                data = tibble(parameter = "b_muminor_reassortantclass_collection_regionnamecentral&northernamerica" ),
                args = list(mean = 0, sd = 5),
                fill = '#d95f02',
                geom = 'area',
                alpha = 0.5) +
  
  stat_function(fun = dnorm,
                data = tibble(parameter = "b_muminor_reassortantclass_collection_regionnameeurope" ),
                args = list(mean = 0, sd = 5),
                fill = '#d95f02',
                geom = 'area',
                alpha = 0.5) +
  
  stat_function(fun = dnorm,
                data = tibble(parameter = "b_muminor_reassortantclass_previous_reassortant_classmajor"),
                args = list(mean = 0, sd = 5),
                fill = '#d95f02',
                geom = 'area',
                alpha = 0.5) +
  
  stat_function(fun = dnorm,
                data = tibble(parameter = "b_muminor_reassortantclass_previous_reassortant_classminor"),
                args = list(mean = 0, sd = 5),
                fill = '#d95f02',
                geom = 'area',
                alpha = 0.5) +
  
  stat_function(fun = dnorm,
                data = tibble(parameter = "b_muminor_reassortantclass_previous_reassortant_classnone"),
                args = list(mean = 0, sd = 5),
                fill = '#d95f02',
                geom = 'area',
                alpha = 0.5) +
  
  
  stat_function(fun = dnorm,
                data = tibble(parameter = "b_munone_reassortantclass_collection_regionnameasia" ),
                args = list(mean = 0, sd = 5),
                fill = '#d95f02',
                geom = 'area',
                alpha = 0.5) +
  
  stat_function(fun = dnorm,
                data = tibble(parameter = "b_munone_reassortantclass_collection_regionnamecentral&northernamerica" ),
                args = list(mean = 0, sd = 5),
                fill = '#d95f02',
                geom = 'area',
                alpha = 0.5) +
  
  stat_function(fun = dnorm,
                data = tibble(parameter = "b_munone_reassortantclass_collection_regionnameeurope" ),
                args = list(mean = 0, sd = 5),
                fill = '#d95f02',
                geom = 'area',
                alpha = 0.5) +
  
  stat_function(fun = dnorm,
                data = tibble(parameter = "b_munone_reassortantclass_previous_reassortant_classmajor"),
                args = list(mean = 0, sd = 5),
                fill = '#d95f02',
                geom = 'area',
                alpha = 0.5) +
  
  stat_function(fun = dnorm,
                data = tibble(parameter = "b_munone_reassortantclass_previous_reassortant_classminor"),
                args = list(mean = 0, sd = 5),
                fill = '#d95f02',
                geom = 'area',
                alpha = 0.5) +
  
  stat_function(fun = dnorm,
                data = tibble(parameter = "b_munone_reassortantclass_previous_reassortant_classnone"),
                args = list(mean = 0, sd = 5),
                fill = '#d95f02',
                geom = 'area',
                alpha = 0.5) +
  
  
  stat_function(fun = dgamma,
                data = tibble(parameter = "shape_nreassortants"),
                args = list(shape = 0.01, rate = 0.01),
                fill = '#d95f02',
                geom = 'area',
                alpha = 0.5) +
  
  
  
  
  xlim(c(-15,20)) + 
  facet_wrap(~parameter, scales = 'free_y',  ncol = 4) +
  theme_minimal(base_size = 8)

############################################## WRITE ###############################################




############################################## END #################################################
####################################################################################################
####################################################################################################