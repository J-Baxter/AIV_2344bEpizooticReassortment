####################################################################################################
####################################################################################################
## Script name: Diffusion Model Interpretation and Figures
##
## Purpose of script: To extract, interpret and display the parameter values and predictions of the
## hurdle model fitted in ./scripts/diffusion_model.R
##
## Date created: 2025-03-20
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

# User functions
scientific_10 <- function(x) {
  parse(text=gsub("e\\+*", " %*% 10^", scales::scientific_format()(x)))
}


############################################## DATA ################################################
diffusionmodel1_fit <- readRDS('./saved_models/diffusion_model_2.rds')
diffusion_data <- read_csv('./saved_models/diffusion_model.csv')



############################################## MAIN ################################################
mcmc_diffusion <- ggs(diffusionmodel1_fit) # Warning message In custom.sort(D$Parameter) : NAs introduced by coercion
posteriorpredictive <-pp_check(diffusionmodel1_fit, ndraws = 500)


# 1. Data plots
means = diffusion_data %>% 
  summarise(mean = median(log(weighted_diff_coeff)), .by = collection_regionname)
diffusion_data %>%
  ggplot(aes(x = log(weighted_diff_coeff))) +
  geom_histogram(binwidth = 0.5) +
  geom_vline(aes(xintercept = mean), data = means) + 
  facet_grid(rows = vars(collection_regionname)) +
  scale_x_continuous(expression(paste('Weighted Diffusion Coefficient (',Km**2~year**-1, ')' )),
                     breaks = log1p(c(0, 10^(seq(from = 1, to = 10, by = 1)))),
                     labels = expression(0, 1%*%10^1, 1%*%10^2, 1%*%10^3, 1%*%10^4,1%*%10^5,1%*%10^6,1%*%10^7,1%*%10^8, 1%*%10^9,1%*%10^10),
                     limits = c(-0.01, log1p(10^10.5)),
                     expand = c(0.02,0.02))
                                                                                                                                                                                                                                   breaks = log1p(c(0, 10^(seq(from = 1, to = 10, by = 1)))),
                                                                                                                                                                                                                                   labels = expression(0, 1%*%10^1, 1%*%10^2, 1%*%10^3, 1%*%10^4,1%*%10^5,1%*%10^6,1%*%10^7,1%*%10^8, 1%*%10^9,1%*%10^10),
                                                                                                                                                                                                                                   limits = c(-0.01, log1p(10^10.5)),
                                                                                                                                                                                                                                   expand = c(0.02,0.02))

# 2. PP Check
pp_check(diffusionmodel1_fit_gamma_5, ndraws = 500,type = 'stat_grouped', group = 'collection_regionname', stat= 'mean')
posteriorpredictive <-pp_check(diffusionmodel1_fit_gamma_5, ndraws = 500 )
posteriorpredictive$data  %>%
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
  scale_x_continuous(expression(paste('Predicted Weighted Diffusion Coefficient (',Km**2~year**-1, ')' )),
                     breaks = log1p(c(0, 10^(seq(from = 2, to = 10, by = 4)))),
                     labels = expression(0,  1%*%10^2,  1%*%10^6, 1%*%10^10),
                     #limits = c(-0.01, log1p(10^9)),
                     expand = c(0.02,0.02))+
  scale_y_continuous('Density',expand = c(0,0)) 


# 2. Plot (average marginal means)



# 3. Fixed effect parameter 'significance'
bayesfactor_parameters(diffusionmodel1_fit_gamma)

avg_slopes(diffusionmodel1_fit_gamma_5, by = "collection_regionname")
# 4. Fixed effect parameter contrasts
diffusionmodel1_fit_gamma_6%>%
  emmeans(~ collection_regionname,
          epred = TRUE, 
          #dpar = "mu",
          re_formula = NA ) %>%
  pairs()

group_diff_prior <- pairs(emmeans(unupdate(diffusionmodel1_fit_gamma_5), ~collection_regionname,  re_formula = NA))

diffusionmodel1_fit_gamma_5%>%
  emmeans(~ collection_regionname,
          #dpar = "mu",
          re_formula = NA ) %>%
  pairs() %>%
  bayesfactor_rope(prior = group_diff_prior)
############################################## WRITE ###############################################

diffusionmodel1_fit_gamma_5 %>%
  emtrends(~ median_charadriiformes_wild,
           var = 'median_charadriiformes_wild',
           at = list(median_charadriiformes_wild = c(0.08, 0.25, 0.5, 1)),
           delta.var = 0.001) %>% 
  gather_emmeans_draws() %>%
  mutate(var = 'median_charadriiformes_wild') %>%
  rename(var_change = median_charadriiformes_wild)


############################################## END #################################################
####################################################################################################
####################################################################################################