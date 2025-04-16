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
# 1. PP Check
#pp_check(diffusionmodel1_fit_gamma_12, ndraws = 100,type = 'stat_grouped', group = 'collection_regionname', stat= 'mean')
posteriorpredictive <-pp_check(diffusionmodel1_fit_gamma_19, ndraws = 100, type = 'dens_overlay_grouped', group = 'collection_regionname' )

# average posterior predictions 
avg_predictions(diffusionmodel1_fit_gamma_17)
avg_predictions(diffusionmodel1_fit_gamma_17, by = 'collection_regionname')
avg_predictions(diffusionmodel1_fit_gamma_17, by = 'collection_regionname', newdata = 'balanced')

# contrasts for continent, all else equal
avg_comparisons(diffusionmodel1_fit_gamma_17, variables = list("collection_regionname" = 'pairwise'))


# 2 the number of species 
# 2a. average posterior predictions marginalised across empirical distribution of predictors
avg_predictions(diffusionmodel1_fit_gamma_17, variables = list('count_cross_species' = c(0,2,4,6,8,10)))

# 2b. average marginal effect at the mean
avg_slopes(diffusionmodel1_fit_gamma_17, variables = 'count_cross_species', newdata = 'balanced')

# 2c. by continent
avg_slopes(diffusionmodel1_fit_gamma_17, variables = 'count_cross_species', by = 'collection_regionname', newdata = 'balanced')


# 3. binary charadriiformes
# 3a.

# 4 median_charadriiformes_wild
# 4a. average posterior predictions marginalised across empirical distribution of predictors
avg_predictions(diffusionmodel1_fit_gamma_17, variables = list('median_charadriiformes_wild' = c(0.08, 0.25, 0.5, 1, 1.5)))

# 4b. average marginal effect at the mean
avg_slopes(diffusionmodel1_fit_gamma_17, variables = 'median_charadriiformes_wild', newdata = 'balanced')

# 4c. by continent
avg_slopes(diffusionmodel1_fit_gamma_17, variables = 'median_charadriiformes_wild', by = 'collection_regionname', newdata = 'balanced')


# 5 binary anseriformes
# 5a.

# median_anseriformes_wild 
# 6a. average posterior predictions
avg_predictions(diffusionmodel1_fit_gamma_17, variables = list('median_anseriformes_wild' = c(0.25, 0.5, 1, 1.5,2)))

# 6b. average marginal effect
avg_slopes(diffusionmodel1_fit_gamma_17, variables = 'median_anseriformes_wild', newdata = 'balanced')

# 6c. by continent
avg_slopes(diffusionmodel1_fit_gamma_17, variables = 'median_anseriformes_wild', by = 'collection_regionname', newdata = 'balanced')


############################################## WRITE ###############################################

# A
plt_5a <- diffusion_data %>%
  ggplot() +
  geom_histogram(aes(x = log1p(weighted_diff_coeff), 
                     fill = weighted_diff_coeff>0,
                     colour = weighted_diff_coeff>0,
                     y = after_stat(density)), 
                 binwidth = 0.5,
                 alpha = 0.7) +
  scale_fill_brewer(palette = 'Accent', 'Is Zero') +
  scale_colour_brewer(palette = 'Accent', 'Is Zero') +
  #scale_x_continuous('Log1p Weighted Diffusion Coefficient' , expand = c(0.01,0.01)) +
  scale_x_continuous(expression(paste('Weighted Diffusion Coefficient (',Km**2~year**-1, ')' )),
                     breaks = log1p(c(0, 10^(seq(from = 1, to = 7, by = 2)))),
                     labels = expression(0,  1%*%10^1,  1%*%10^3, 1%*%10^5,  1%*%10^7),
                     expand = c(0.02,0.02))+
  
  scale_y_continuous('Probability Density',expand = c(0,0)) + 
  global_theme + 
  theme(legend.position = 'none')

# B - anseriformes persistence
plt_5b <- diffusion_data %>%
  ggplot() +
  geom_histogram(aes(x = median_anseriformes_wild,
                     y = after_stat(density)), 
                 binwidth = 0.2,
                 fill = host_colours['anseriformes-wild'],
                 colour = host_colours['anseriformes-wild'],
                 alpha = 0.7) +
  scale_x_continuous('Persistence in wild Anseriformes', 
                     expand = c(0.02,0.02), 
                     breaks = seq(from = 0, to = 10, by = 2)) +
  scale_y_continuous('Probability Density',expand = c(0,0)) + 
  coord_cartesian(xlim = c(0, 10))  +
  global_theme + 
  theme(legend.position = 'none')


# C - charadriiformes persistence
plt_5c  <- diffusion_data %>%
  ggplot() +
  geom_histogram(aes(x = median_charadriiformes_wild,
                     y = after_stat(density)),
                 binwidth = 0.2,
                 fill = host_colours['charadriiformes-wild'],
                 colour = host_colours['charadriiformes-wild'],
                 alpha = 0.7) +
  scale_x_continuous('Persistence in wild Charadriiformes', 
                     expand = c(0.01,0.01), 
                     breaks = seq(from = 0, to = 4, by = 1),
                     limits = c(-0.1,4.5)) +
  scale_y_continuous('Probability Density',expand = c(0,0)) + 
  #coord_cartesian(xlim = c(0, 6))  +
  global_theme + 
  theme(legend.position = 'none')


# D - Host switches
plt_5d  <- diffusion_data %>%
  ggplot() +
  geom_histogram(aes(x = log1p(count_cross_species),
                     y = after_stat(density)),
                 binwidth = 0.5,
                 fill = '#beaed4',
                 colour = '#beaed4',

                 alpha = 0.7) +
  scale_x_continuous('Independent Host Jumps (Ln)', 
                     expand = c(0.01,0.01), 
                     breaks = seq(from = 0, to = 6, by = 1)) +
  scale_y_continuous('Probability Density',expand = c(0,0)) + 
  #coord_cartesian(xlim = c(0, 6))  +
  global_theme + 
  theme(legend.position = 'none')


# E - Posterior predictive draws
plt_5e  <- diffusionmodel1_fit_gamma_17 %>%
  #  back-transformed linear predictive draws, equivalent to add_linpred (and in this case equviv to epred)
  # this uses the empirical data distribution by default
  avg_predictions(by = 'collection_regionname',  type = 'response') %>% 
  get_draws() %>%
  as_tibble() %>%
  ggplot() + 
  geom_histogram(aes(x = log1p(draw), y = after_stat(density), colour = collection_regionname, fill = collection_regionname), 
                 binwidth = 0.1,
                 alpha = 0.7) + 
  scale_colour_manual(values = region_colours)+
  scale_fill_manual(values = region_colours) + 
  scale_x_continuous(expression(paste('Predicted Weighted Diffusion Coefficient (',Km**2~year**-1, ')' )),
                     breaks =log(10^seq(5, 6.5, b = 0.5)),
                     labels = expression(1%*%10^5, 1%*%10^5.5,  1%*%10^6, 1%*%10^6.5),
                     limits = c(11.5, 15),
                     expand = c(0.02,0.02)
                     )+
  scale_y_continuous('Probability Density' ,
                     expand = c(0,0))+
  facet_grid(
    cols = vars(collection_regionname),
    labeller =  labeller(collection_regionname=str_to_title)) +
  #geom_vline(aes(xintercept = emmean, colour = collection_regionname), data = averages, linetype = 'dashed') +
  #geom_text(aes(label =  paste0("E*'('*X*'|'*X*'>'*0*') = '*", label, "~km^2"), 
          #      colour = collection_regionname),
         #   parse = T,
           # x = 17.5, 
          #  y = 0.6,
           # size = 2.5,
           # data = averages) + 
  global_theme + 
  theme(strip.placement  = 'inside',
        strip.text = element_text(face = 'bold', size = 10),
        strip.background = element_blank(),
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 10),
        legend.text = element_text(size = 8))


# F - Persistence Effect

diffusionmodel1_fit_gamma_17 %>%
            #  back-transformed linear predictive draws, equivalent to add_linpred (and in this case equviv to epred)
            # this uses the empirical data distribution by default
            avg_slopes(variables = c('median_charadriiformes_wild','median_anseriformes_wild'), by = 'collection_regionname', newdata = 'balanced') %>%
            get_draws() %>%
            as_tibble()

diffusionmodel1_fit_gamma_17 %>%
  #  back-transformed linear predictive draws, equivalent to add_linpred (and in this case equviv to epred)
  # this uses the empirical data distribution by default
  #avg_predictions(variables = list('median_charadriiformes_wild' = c(0.08, 0.25, 0.5, 1, 1.5)), by = c('collection_regionname', 'median_charadriiformes_wild'))  %>%
  avg_predictions(diffusionmodel1_fit_gamma_17, variables = list('count_cross_species' = c(0,2,4,6,8,10)), by = c('collection_regionname', 'count_cross_species')) %>%

  get_draws() %>%
  as_tibble() %>%
  ggplot(aes(y = draw,
             x = count_cross_species,
             colour = collection_regionname)) + 
  geom_point(alpha = 0.2) + 
  stat_smooth(method = 'glm') + 
  facet_grid(
    cols = vars(collection_regionname),
    labeller =  labeller(collection_regionname=str_to_title))


  
  ggplot(aes(x  = draw,
           y = as.factor(term),
           slab_colour = var,
           slab_fill = var)) +
  stat_halfeye(slab_alpha = 0.7 ,
               p_limits = c(0.001, 0.999),
               point_interval = "median_hdi",
               linewidth = 1.5,
               .width =  0.95)  +
  geom_vline(aes(xintercept = 0), linetype = 'dashed')  +
  scale_x_continuous(expression(paste('Change in Weighted Diffusion Coefficient (',Km**2~year**-1, ')' )),
                     breaks = seq(from = -1*10**6, to = 5*10**6, by = 2*10**6),
                     labels = scientific_10,
                     limits = c(-1*10**6, 5*10**6),
                     expand = c(0.02,0.02))+
  scale_fill_manual(values = host_colours, aesthetics = 'slab_fill') +
  scale_colour_manual(values = host_colours, aesthetics = 'slab_colour') + 
  scale_y_discrete('Persistence in Host (Years)', 
                   labels = function(x) str_to_title(x) %>% str_wrap(., width = 10)) +
  global_theme + 
  theme(legend.position = 'none')

  ggplot() + 
  geom_histogram(aes(x = log1p(draw), y = after_stat(density), colour = collection_regionname, fill = collection_regionname), 
                 binwidth = 0.1,
                 alpha = 0.7) + 
  scale_colour_manual(values = region_colours)+
  scale_fill_manual(values = region_colours) + 
  scale_x_continuous(expression(paste('Predicted Weighted Diffusion Coefficient (',Km**2~year**-1, ')' )),
                     breaks =log(10^seq(5, 6.5, b = 0.5)),
                     labels = expression(1%*%10^5, 1%*%10^5.5,  1%*%10^6, 1%*%10^6.5),
                     limits = c(11.5, 15),
                     expand = c(0.02,0.02)
  )+
  scale_y_continuous('Probability Density' ,
                     expand = c(0,0))+
  facet_grid(
    cols = vars(collection_regionname),
    labeller =  labeller(collection_regionname=str_to_title)) +
  #geom_vline(aes(xintercept = emmean, colour = collection_regionname), data = averages, linetype = 'dashed') +
  #geom_text(aes(label =  paste0("E*'('*X*'|'*X*'>'*0*') = '*", label, "~km^2"), 
  #      colour = collection_regionname),
  #   parse = T,
  # x = 17.5, 
  #  y = 0.6,
  # size = 2.5,
  # data = averages) + 
  global_theme + 
  theme(strip.placement  = 'inside',
        strip.text = element_text(face = 'bold', size = 10),
        strip.background = element_blank(),
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 10),
        legend.text = element_text(size = 8))


# G - Host Jump Effect


# H - Marginal 



# Supplementary plots


############################################## END #################################################
####################################################################################################
####################################################################################################