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
diffusionmodel_fit <- readRDS('./saved_models/diffusion_model_2.rds')
diffusion_data <- read_csv('./saved_models/diffusion_model.csv')


############################################## MAIN ################################################
# 1. PP Check
#pp_check(diffusionmodel1_fit_gamma_12, ndraws = 100,type = 'stat_grouped', group = 'collection_regionname', stat= 'mean')
posteriorpredictive <-pp_check(diffusionmodel1_fit_gamma_16, ndraws = 100, type = 'dens_overlay_grouped', group = 'collection_regionname' )

# average posterior predictions 
avg_predictions(diffusionmodel1_fit_gamma_16)
avg_predictions(diffusionmodel1_fit_gamma_17, by = 'collection_regionname')
avg_predictions(diffusionmodel1_fit_gamma_17, by = 'collection_regionname', newdata = 'balanced')

# contrasts for continent, all else equal
avg_comparisons(diffusionmodel1_fit_gamma_17, variables = list("collection_regionname" = 'pairwise'))


# the number of species 
# 2a. average posterior predictions
avg_predictions(diffusionmodel1_fit_gamma_17, variables = list('count_cross_species' = c(0,2,4,6,8,10)))

# 2b. average marginal effect
avg_slopes(diffusionmodel1_fit_gamma_14, variables = 'count_cross_species')

# 2c. by continent
avg_slopes(diffusionmodel1_fit_gamma_14, variables = 'count_cross_species', by = 'collection_regionname')


# median_charadriiformes_wild
# 3a. average posterior predictions
avg_predictions(diffusionmodel1_fit_gamma_17, variables = list('median_charadriiformes_wild' = c(0.08, 0.25, 0.5, 1, 1.5)))

# 3b. average marginal effect
avg_slopes(diffusionmodel1_fit_gamma_14, variables = 'median_charadriiformes_wild')

# 3c. by continent
avg_slopes(diffusionmodel1_fit_gamma_14, variables = 'median_charadriiformes_wild', by = 'collection_regionname')



# median_anseriformes_wild
# 3a. average posterior predictions
avg_predictions(diffusionmodel1_fit_gamma_17, variables = list('median_anseriformes_wild' = c(0.25, 0.5, 1, 1.5,2)))

# 3b. average marginal effect
avg_slopes(diffusionmodel1_fit_gamma_14, variables = 'median_anseriformes_wild')

# 3c. by continent
avg_slopes(diffusionmodel1_fit_gamma_14, variables = 'median_anseriformes_wild', by = 'collection_regionname')


############################################## WRITE ###############################################
posteriorpredictive$data  %>%
  mutate(value = log1p(value)) %>%
  #filter(!is_y) %>%
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
         colour=guide_legend()) + facet_wrap(~group) + 
  scale_x_continuous(expression(paste('Predicted Weighted Diffusion Coefficient (',Km**2~year**-1, ')' )),
                     breaks = log1p(c(0, 10^(seq(from = 2, to = 10, by = 4)))),
                     labels = expression(0,  1%*%10^2,  1%*%10^6, 1%*%10^10),
                     limits = c(-0.01, log1p(10^9)),
                     expand = c(0.02,0.02))+
  scale_y_continuous('Density',expand = c(0,0)) 




############################################## END #################################################
####################################################################################################
####################################################################################################