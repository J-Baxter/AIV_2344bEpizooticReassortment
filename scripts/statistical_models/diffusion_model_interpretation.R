####################################################################################################
####################################################################################################
## Script name: Diffusion Model Interpretation and Figures
##
## Purpose of script: To extract, interpret and display the parameter values and predictions of the
## hurdle model fitted in ./scripts/diffusion_model.R
##
## Date created: 2025-01-22
##
##
########################################## SYSTEM OPTIONS ##########################################
options(scipen = 6, digits = 7) 
memory.limit(30000000) 
options('marginaleffects_posterior_interval' = 'hdi')
  
########################################## DEPENDENCIES ############################################
# Packages
library(tidyverse)
library(magrittr)
library(brms)
library(broom)
library(broom.mixed)
library(tidybayes) # missing
library(bayesplot)
library(emmeans) # missing
library(marginaleffects) # missing
library(magrittr)
library(ggmcmc) # missing
library(bayestestR)
library(modelbased)

# User functions
scientific_10 <- function(x) {
  parse(text=gsub("e\\+*", " %*% 10^", scales::scientific_format()(x)))
  }

############################################## DATA ################################################
#diffusionmodel_fit <- readRDS('./saved_models/diffusion_model_2.rds')sw 
#diffusion_data <- read_csv('./saved_models/diffusion_model.csv')

host_colours <- c(
  'anseriformes-domestic' = '#a6cee3',
  'anseriformes-wild' = '#1f78b4',
  'galliformes-domestic' = '#b2df8a',
  'galliformes-wild' = '#33a02c',
  'mammal' = '#fb9a99',
  'human' = '#e31a1c',
  'charadriiformes-wild' = '#fdbf6f',
  'other-bird' = '#ff7f00',
  'unknown' = '#cab2d6',
  'environment' = '#6a3d9a')


region_colours <- c('europe' = '#1b9e77',
                    'asia' ='#d95f02',
                    'africa' ='#7570b3',
                    'australasia' = '#e7298a',
                    'central & northern america' ='#66a61e',
                    'south america' ='#e6ab02')


global_theme <- theme_classic()+
  theme(
    #text = element_text(size=10),
    axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0), size = 10),
    axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0), size = 10),
    axis.text = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0), size = 8),
    legend.text = element_text(size = 8),
    legend.position = 'none', 
    panel.spacing = unit(2, "lines"), 
    strip.background = element_blank()
  )


############################################## MAIN ################################################
# 1. PP Check - deprecated to evaluation script
#pp_check(diffusionmodel1_fit_gamma_12, ndraws = 100,type = 'stat_grouped', group = 'collection_regionname', stat= 'mean')
#posteriorpredictive <-pp_check(diffusionmodel1_fit_gamma_19, ndraws = 100, type = 'dens_overlay_grouped', group = 'collection_regionname' )

# P2
avg_predictions(diffusionmodel1_fit_gamma_19)
avg_predictions(diffusionmodel1_fit_gamma_19, by = 'collection_regionname')

#P3
avg_predictions(diffusionmodel1_fit_gamma_19, 
                variables = list('count_cross_species_log1p' = log1p(0:10)))


diffusionmodel1_fit_gamma_19 %>%
  avg_comparisons(variables = list("count_cross_species_log1p" = log1p(1)), # a 1 month increase
                  newdata = datagrid(count_cross_species_log1p = log1p(c(0,1,2,3,4,5,6,7,8,9,10)), persist.time_log1p = log1p(1),
                                     grid_type = 'counterfactual'), 
                  type = 'response')

avg_predictions(diffusionmodel1_fit_gamma_19, by = 'collection_year')


# P4
avg_predictions(diffusionmodel1_fit_gamma_19, 
                variables = list('path_prop_charadriiformes_wild' = seq(0,1,by = 0.1),
                                 'persistence' = 1),
                by = c('path_prop_charadriiformes_wild')) %>%
  as_tibble()

avg_predictions(diffusionmodel1_fit_gamma_19, 
                variables = list('path_prop_anseriformes_wild' = seq(0, 1, by = 0.1), 'persistence' = 1),
                by = c('path_prop_anseriformes_wild'))%>%
  as_tibble()

avg_predictions(diffusionmodel1_fit_gamma_19, 
                variables = list('path_prop_galliformes_domestic' = seq(0, 1, by = 0.1), 'persistence' = 1),
                by = c('path_prop_galliformes_domestic'))%>%
  as_tibble()


#P5
diffusionmodel1_fit_gamma_19 %>%
  avg_comparisons(variables = list("path_prop_anseriformes_wild" = 0.05), 
                  by = "collection_regionname",
                  newdata =datagrid(collection_regionname = unique(diffusion_data$collection_regionname,
                                                                   persist.time_log1p = log1p(1)),
                                    grid_type = 'counterfactual'), 
                  type = 'response')

diffusionmodel1_fit_gamma_19 |> 
  avg_comparisons(variables = list("path_prop_anseriformes_wild" = 0.05), 
                  by = "collection_regionname",
                  newdata =datagrid(collection_regionname = unique(diffusion_data$collection_regionname), grid_type = 'counterfactual'), 
                  type = 'response')

############################################## WRITE ###############################################
diffusionmodel1_fit_gamma_19 %>%
  avg_comparisons(variables = list('path_prop_anseriformes_wild' = 0.1), 
                #by = "segments_changed",
                type = 'response')
############################################## END #################################################
####################################################################################################
####################################################################################################