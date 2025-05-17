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

# average posterior predictions 
avg_predictions(diffusionmodel1_fit_gamma_19)
avg_predictions(diffusionmodel1_fit_gamma_19, by = 'collection_regionname')
avg_predictions(diffusionmodel1_fit_gamma_19, by = 'collection_regionname', newdata = 'balanced')

# contrasts for continent, all else equal
avg_comparisons(diffusionmodel1_fit_gamma_19, variables = list("collection_regionname" = 'pairwise'))


# 2 the number of species 
# 2a. average posterior predictions marginalised across empirical distribution of predictors
avg_predictions(diffusionmodel1_fit_gamma_19, variables = list('count_cross_species' = c(0,2,4,6,8,10)))

# 2b. average marginal effect at the mean
avg_slopes(diffusionmodel1_fit_gamma_19, variables = 'count_cross_species', newdata = 'balanced')

# 2c. by continent
avg_slopes(diffusionmodel1_fit_gamma_19, variables = 'count_cross_species', by = 'collection_regionname', newdata = 'balanced')


# 3. binary charadriiformes
# 3a.

# 4 median_charadriiformes_wild
# 4a. average posterior predictions marginalised across empirical distribution of predictors
avg_predictions(diffusionmodel1_fit_gamma_19, variables = list(list('median_charadriiformes_wild_log1p' = log1p(c(0.08, 0.25, 0.5, 1, 1.5)))))

# 4b. average marginal effect at the mean
avg_slopes(diffusionmodel1_fit_gamma_19, variables = 'median_charadriiformes_wild_log1p', newdata = 'balanced')

# 4c. by continent
avg_slopes(diffusionmodel1_fit_gamma_19, variables = 'median_charadriiformes_wild_log1p', by = 'collection_regionname', newdata = 'balanced')


# 5 binary anseriformes
# 5a.

# median_anseriformes_wild 
# 6a. average posterior predictions
avg_predictions(diffusionmodel1_fit_gamma_19, variables = list('median_anseriformes_wild' = c(0.25, 0.5, 1, 1.5,2)))

# 6b. average marginal effect
avg_slopes(diffusionmodel1_fit_gamma_19, variables = 'median_anseriformes_wild', newdata = 'balanced')

# 6c. by continent
avg_slopes(diffusionmodel1_fit_gamma_19, variables = 'median_anseriformes_wild_log1p', by = 'collection_regionname', newdata = 'balanced')

slopes(diffusionmodel1_fit_gamma_19, newdata = datagrid(median_anseriformes_wild_log1p = log1p(c(0.08, 0.25, 0.5, 1, 1.5)), grid_type = 'counterfactual'), variables = 'median_anseriformes_wild_log1p', by =c('median_anseriformes_wild_log1p', 'collection_regionname'))
slopes(diffusionmodel1_fit_gamma_19, newdata = datagrid(median_charadriiformes_wild_log1p = log1p(c(0.08, 0.25, 0.5, 1, 1.5)), grid_type = 'counterfactual'), variables = 'median_charadriiformes_wild_log1p', by =c('median_charadriiformes_wild_log1p', 'collection_regionname'))
slopes(diffusionmodel1_fit_gamma_19, newdata = datagrid(count_cross_species_log1p = log1p(c(1, 10, 25, 50)), grid_type = 'counterfactual'), variables = 'count_cross_species_log1p', by =c('count_cross_species_log1p'), eps = 0.001)

diffusionmodel1_fit_gamma_19 |> 
  avg_comparisons(variables = list("median_charadriiformes_wild_prop" = 0.05), 
                  by = "collection_regionname",
                  newdata =datagrid(collection_regionname = unique(diffusion_data$collection_regionname), grid_type = 'counterfactual'), 
                  type = 'response')

diffusionmodel1_fit_gamma_19 |> 
  avg_comparisons(variables = list("median_anseriformes_wild_prop" = 0.05), 
                  by = "collection_regionname",
                  newdata =datagrid(collection_regionname = unique(diffusion_data$collection_regionname), grid_type = 'counterfactual'), 
                  type = 'response')

diffusionmodel1_fit_gamma_19 |> 
  avg_comparisons(variables = list("persist.time_log1p" = log1p(0.1)), 
                  by = "persist.time_log1p",
                  newdata =datagrid(persist.time_log1p =  log1p(c(0.25,0.5,1,2)),
                                    grid_type = 'counterfactual'), 
                  type = 'response')



class_model |> 
  avg_comparisons(variables = list("time_since_previous" = 0.08), # a 1 month increase
                  by = "time_since_previous", 
                  newdata = datagrid(time_since_previous = c(0.08, 0.25), 
                                     grid_type = 'counterfactual'), 
                  type = 'response')

diffusionmodel1_fit_gamma_19 %>% plot_predictions(condition = c("median_anseriformes_wild_prop"))
diffusionmodel1_fit_gamma_19 %>% plot_predictions(condition = c("median_charadriiformes_wild_prop" ))
diffusionmodel1_fit_gamma_19 %>% plot_predictions(condition = c("persist.time_log1p", 'collection_regionname')) + geom_point(data = diffusion_data, aes(x = persist.time_log1p, y = weighted_diff_coeff), inherit.aes = F) + scale_y_continuous(trans = 'log10')
diffusionmodel1_fit_gamma_19 %>% plot_predictions(condition = c("count_cross_species_log1p"))


diffusionmodel1_fit_gamma_19 |> 
  avg_comparisons(variables = list("count_cross_species_log1p" = log1p(1)), # a 1 month increase
                  newdata = datagrid(count_cross_species_log1p = log1p(c(0,1,2,3,4,5,6,7,8,9,10)), persist.time_log1p = log1p(1),
                                     grid_type = 'counterfactual'), 
                  type = 'response')


avg_predictions(diffusionmodel1_fit_gamma_19)

plot_predictions(class_model, condition = c("time_since_previous", "collection_regionname")) + facet_wrap(.~group)

# a 1% increase in the time in X leads the following change in diffusion coefficient at a fixed persistence time
diffusionmodel1_fit_gamma_19 |> 
  avg_comparisons(variables = list("median_anseriformes_wild_prop" = 0.05), # a 1 month increase
                  newdata = datagrid(median_anseriformes_wild_prop = c(0,0.25,0.5,1), persist.time_log1p = log1p(1), collection_regionname = diffusion_data %>% pull(collection_regionname) %>% unique(), 
                                     grid_type = 'counterfactual'), by = c('collection_regionname', 'median_anseriformes_wild_prop'),
                  type = 'response')

avg_predictions(diffusionmodel1_fit_gamma_19, by = "count_cross_species_log1p")


diffusionmodel1_fit_gamma_20 %>%
  emtrends(~ median_charadriiformes_wild_prop + persist.time_log1p + collection_regionname,
           var = "median_charadriiformes_wild_prop",
         at = list(median_charadriiformes_wild_prop = c(0,0.25,0.5,1), 
                   persist.time_log1p = log1p(0.5),
                   collection_regionname = diffusion_data %>% pull(collection_regionname) %>% unique() %>% droplevels()),
         regrid = "response", delta.var = 0.01) 


plot_slopes(diffusionmodel1_fit_gamma_20, # model
       newdata = datagrid(median_anseriformes_wild_prop = c(0,0.25,0.5,1), 
                          persist.time_log1p = log1p( 0.5),
                          grid_type = 'counterfactual'), 
       variables = 'median_anseriformes_wild_prop',
       by =c('median_anseriformes_wild_prop', 'collection_regionname'))
############################################## WRITE ###############################################

############################################## END #################################################
####################################################################################################
####################################################################################################