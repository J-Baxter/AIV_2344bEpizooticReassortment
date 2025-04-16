####################################################################################################
####################################################################################################
## Script name: Diffusion Model Evaluation
##
## Purpose of script:
##
## Date created: 2025-04-11
##
##
########################################## SYSTEM OPTIONS ##########################################
options(scipen = 6, digits = 7) 
memory.limit(30000000) 

  
########################################## DEPENDENCIES ############################################
# Packages
library(tidyverse)
library(magrittr)
library(broom.mixed)
library(tidybayes)
library(bayesplot)
library(scales)
library(DHARMa)
library(ggdist)
library(ggmcmc)


# User functions

######################################### DATA & MODEL #############################################
fitted_model <- diffusionmodel1_fit_gamma_19

diffusion_data <- diffusion_data
###################################### MCMC Diagnostics #############################################

#MCMC chain convergence (Visual Inspection)
# Plot Chains
ggs(diffusionmodel1_fit_gamma_19) %>% 
  from_ggmcmc_names() %>%
  filter(.iteration > 400) %>%
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
                               .variable == 'b_shape_collection_regionnamecentral&northernamerica'~ "alpha['americas']",
                               
                               
                               .variable == 'sd_segment__Intercept'~ "sigma[gamma]",
                               .variable == 'sd_collection_regionname__shape_Intercept'~ "sigma[zeta]",
                               
                               .variable == 'r_segment[ha,Intercept] '~ "gamma['ha']",
                               .variable == "r_segment[mp,Intercept]"~ "gamma['mp']",
                               .variable == "r_segment[np,Intercept]"~ "gamma['np']",
                               .variable == "r_segment[ns,Intercept]"~ "gamma['ns']",
                               .variable == "r_segment[nx,Intercept]"~ "gamma['nx']",
                               .variable == "r_segment[pa,Intercept]"~ "gamma['pa']",
                               .variable == "r_segment[pb1,Intercept]"~ "gamma['pb1']",
                               .variable == "r_segment[pb2,Intercept]"~ "gamma['pb2']",
                               .variable == 'r_collection_regionname__shape[asia,Intercept]'~ "zeta['asia']",
                               .variable == 'r_collection_regionname__shape[africa,Intercept]'~ "zeta['africa']",
                               .variable == 'r_collection_regionname__shape[central.&.northern.america,Intercept]'~ "zeta['europe']",
                               .variable == 'collection_regionname__shape_central.&.northern.america'~ "zeta['central&northernamerica']",
                               
                               .variable == 'lprior' ~ 'prior',
                               .variable == 'lp__' ~ 'log~probability')) %>%
  drop_na(label) %>%
  ggplot(aes(x = .iteration,
             y = .value, 
             col = as.factor(.chain)))+
  geom_line(alpha = 0.8)+
  facet_wrap(~ label,
             ncol = 2,
             scale  = 'free_y',
             strip.position = 'left',
             labeller = label_parsed)+
  scale_colour_brewer(palette = 'GnBu', 'Chains') +
  scale_x_continuous('Iteration') + 
  scale_y_continuous('Parameter Value') + 
  theme_minimal() + 
  theme(legend.position = 'bottom',
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 10),
        legend.text = element_text(size = 8))

ggsave('~/Downloads/flu_plots/diffusion_trace.jpeg',
       dpi = 360,
       height = 38,
       width = 20,
       units = 'cm')

# Plot ranked traces
# If chains are exploring the same space efficiently, the traces should be similar to one another 
# and largely overlapping.
as_draws_df(diffusionmodel1_fit_gamma_19) %>% 
  
  # Selecting only beta coefficients
  mcmc_rank_overlay(regex_pars = '^b_') %>%
  
  # Extract data
  .$data %>%
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
                           .variable == 'b_shape_collection_regionnamecentral&northernamerica'~ "alpha['americas']",
                           
                           
                           .variable == 'sd_segment__Intercept'~ "sigma[gamma]",
                           .variable == 'sd_collection_regionname__shape_Intercept'~ "sigma[zeta]",
                           
                           .variable == 'r_segment[ha,Intercept] '~ "gamma['ha']",
                           .variable == "r_segment[mp,Intercept]"~ "gamma['mp']",
                           .variable == "r_segment[np,Intercept]"~ "gamma['np']",
                           .variable == "r_segment[ns,Intercept]"~ "gamma['ns']",
                           .variable == "r_segment[nx,Intercept]"~ "gamma['nx']",
                           .variable == "r_segment[pa,Intercept]"~ "gamma['pa']",
                           .variable == "r_segment[pb1,Intercept]"~ "gamma['pb1']",
                           .variable == "r_segment[pb2,Intercept]"~ "gamma['pb2']",
                           .variable == 'r_collection_regionname__shape[asia,Intercept]'~ "zeta['asia']",
                           .variable == 'r_collection_regionname__shape[africa,Intercept]'~ "zeta['africa']",
                           .variable == 'r_collection_regionname__shape[central.&.northern.america,Intercept]'~ "zeta['europe']",
                           .variable == 'collection_regionname__shape_central.&.northern.america'~ "zeta['central&northernamerica']",
                           
                           .variable == 'lprior' ~ 'prior',
                           .variable == 'lp__' ~ 'log~probability')) %>%
  
  # Plot 
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


#MCMC chain resolution (ESS)
model_resolution <- list(tidy(diffusionmodel1_fit_gamma_19, ess = T, rhat = T, effects = 'ran_vals') %>%
                           dplyr::select( ess, rhat, group, level) %>%
                           unite(term, group, level, sep = '_'),
                         tidy(diffusionmodel1_fit_gamma_19, ess = T, rhat = T) %>%  dplyr::select(term, ess, rhat,term) ) %>%
  bind_rows() %>%
  
  mutate(term = case_when(term == 'collection_regionnameasia'~ "beta['asia']",
                          term == 'collection_regionnameafrica'~ "beta['africa']",
                          term == 'collection_regionnameeurope'~ "beta['europe']",
                          term == "collection_regionnamecentral&northernamerica"~ "beta['americas']",
                           
                          term == 'int_stepmedian_anseriformes_wild_log1p'~ "beta['step_anseriformes']",
                          term == 'median_anseriformes_wild_log1p'~ "beta['persist_anseriformes']",
                          term == 'int_stepmedian_charadriiformes_wild_log1p'~ "beta['step_charadriiformes']",
                          term == 'median_charadriiformes_wild_log1p'~ "beta['persist_charadriiformes']",
                          term == 'int_stepcount_cross_species_log1p'~ "beta['step_hostjump']",
                          term == 'count_cross_species_log1p'~ "beta['num_hostjump']",
                           
                          term == 'collection_regionnameafrica:median_anseriformes_wild_log1p'~ "beta['africa-anser']",
                          term == 'collection_regionnameeurope:median_anseriformes_wild_log1p'~ "beta['europe-anser']",
                          term == 'collection_regionnamecentral&northernamerica:median_anseriformes_wild_log1p'~ "beta['americas-anser']",
                           
                          term == 'collection_regionnameafrica:median_charadriiformes_wild_log1p'~ "beta['africa-charad']",
                          term == 'collection_regionnameeurope:median_charadriiformes_wild_log1p'~ "beta['europe-charad']",
                          term == 'collection_regionnamecentral&northernamerica:median_charadriiformes_wild_log1p'~ "beta['americas-charad']",
                           
                          term == 'shape_collection_regionnameasia'~ "alpha['asia']",
                          term == 'shape_collection_regionnameafrica'~ "alpha['africa']",
                          term == 'shape_collection_regionnameeurope'~ "alpha['europe']",
                          term == 'shape_collection_regionnamecentral&northernamerica'~ "alpha['americas']",
                          
                          term == 'sd__(Intercept)'~ "sigma[gamma]",
                          term == 'sd__shape_(Intercept)'~ "sigma[zeta]",
                          
                          term == 'segment_ha'~ "gamma['ha']",
                          term == 'segment_mp'~ "gamma['mp']",
                          term == 'segment_np'~ "gamma['np']",
                          term == 'segment_ns'~ "gamma['ns']",
                          term == 'segment_nx'~ "gamma['nx']",
                          term == 'segment_pa'~ "gamma['pa']",
                          term == 'segment_pb1'~ "gamma['pb1']",
                          term == 'segment_pb2'~ "gamma['pb2']",
                          term == 'collection_regionname__shape_asia'~ "zeta['asia']",
                          term == 'collection_regionname__shape_africa'~ "zeta['africa']",
                          term == 'collection_regionname__shape_europe'~ "zeta['europe']",
                          term == 'collection_regionname__shape_central.&.northern.america'~ "zeta['central&northernamerica']")) 
  
 

model_resolution %>% 
  ggplot(aes(x = term, y = ess)) + 
  geom_bar(stat = 'identity') + 
  scale_x_discrete(labels= label_parse(), 'Parameter', expand = c(0.05, 0)) + 
  scale_y_continuous(expand = c(0,0), 'Effective Sample Size')


# NB - technically a convergence check
model_resolution %>% 
  ggplot(aes(x = term, y = rhat)) + 
  geom_bar(stat = 'identity') + 
  scale_x_discrete(labels= label_parse(), 'Parameter', expand = c(0.05, 0)) + 
  scale_y_continuous(expand = c(0,0), expression(paste("Potential Reduction in Scale Factor (", hat(R), ')')))


#MCMC Autocorrelation
diffusionmodel1_fit_gamma_19 %>% 
  mcmc_acf() %>% 
  .$data %>%
  as_tibble() %>%
  mutate(Parameter = case_when(Parameter == 'b_collection_regionnameasia'~ "beta['asia']",
                               Parameter == 'b_collection_regionnameafrica'~ "beta['africa']",
                               Parameter == 'b_collection_regionnameeurope'~ "beta['europe']",
                               Parameter == "b_collection_regionnamecentral&northernamerica"~ "beta['americas']",
                               
                               Parameter == 'b_int_stepmedian_anseriformes_wild_log1p'~ "beta['step_anseriformes']",
                               Parameter == 'b_median_anseriformes_wild_log1p'~ "beta['persist_anseriformes']",
                               Parameter == 'b_int_stepmedian_charadriiformes_wild_log1p'~ "beta['step_charadriiformes']",
                               Parameter == 'b_median_charadriiformes_wild_log1p'~ "beta['persist_charadriiformes']",
                               Parameter == 'b_int_stepcount_cross_species_log1p'~ "beta['step_hostjump']",
                               Parameter == 'b_count_cross_species_log1p'~ "beta['num_hostjump']",
                               
                               Parameter == 'b_collection_regionnameafrica:median_anseriformes_wild_log1p'~ "beta['africa-anser']",
                               Parameter == 'b_collection_regionnameeurope:median_anseriformes_wild_log1p'~ "beta['europe-anser']",
                               Parameter == 'b_collection_regionnamecentral&northernamerica:median_anseriformes_wild_log1p'~ "beta['americas-anser']",
                               
                               Parameter == 'b_collection_regionnameafrica:median_charadriiformes_wild_log1p'~ "beta['africa-charad']",
                               Parameter == 'b_collection_regionnameeurope:median_charadriiformes_wild_log1p'~ "beta['europe-charad']",
                               Parameter == 'b_collection_regionnamecentral&northernamerica:median_charadriiformes_wild_log1p'~ "beta['americas-charad']",
                               
                               Parameter == 'b_shape_collection_regionnameasia'~ "alpha['asia']",
                               Parameter == 'b_shape_collection_regionnameafrica'~ "alpha['africa']",
                               Parameter == 'b_shape_collection_regionnameeurope'~ "alpha['europe']",
                               Parameter == 'b_shape_collection_regionnamecentral&northernamerica'~ "alpha['americas']",
                               
                               
                               Parameter == 'sd_segment__Intercept'~ "sigma[gamma]",
                               Parameter == 'sd_collection_regionname__shape_Intercept'~ "sigma[zeta]",
                               
                               Parameter == 'r_segment[ha,Intercept] '~ "gamma['ha']",
                               Parameter == "r_segment[mp,Intercept]"~ "gamma['mp']",
                               Parameter == "r_segment[np,Intercept]"~ "gamma['np']",
                               Parameter == "r_segment[ns,Intercept]"~ "gamma['ns']",
                               Parameter == "r_segment[nx,Intercept]"~ "gamma['nx']",
                               Parameter == "r_segment[pa,Intercept]"~ "gamma['pa']",
                               Parameter == "r_segment[pb1,Intercept]"~ "gamma['pb1']",
                               Parameter == "r_segment[pb2,Intercept]"~ "gamma['pb2']",
                               Parameter == 'r_collection_regionname__shape[asia,Intercept]'~ "zeta['asia']",
                               Parameter == 'r_collection_regionname__shape[africa,Intercept]'~ "zeta['africa']",
                               Parameter == 'r_collection_regionname__shape[central.&.northern.america,Intercept]'~ "zeta['europe']",
                               Parameter == 'collection_regionname__shape_central.&.northern.america'~ "zeta['central&northernamerica']",
                               
                               Parameter == 'lprior' ~ 'prior',
                               Parameter == 'lp__' ~ 'log~probability')) %>%
  drop_na(Parameter) %>%
  ggplot(aes(y = AC, x = Lag, colour = as.factor(Chain))) + 
  geom_path() + 
  facet_wrap(~Parameter, labeller = label_parsed) + 
  theme_minimal()  + 
  scale_colour_brewer('Chain') +
  theme(legend.position = 'bottom',
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 10),
        legend.text = element_text(size = 8))

ggsave('~/Downloads/flu_plots/diffusion_autocorrelation.jpeg',
       dpi = 360,
       height = 30,
       width = 30,
       units = 'cm')


# Check Ratio of Effective Population Size to Total Sample Size 
# values <0.1 should raise concerns about autocorrelation
neff_ratio(diffusionmodel1_fit_gamma_19) %>% 
  as_tibble(rownames = 'param') %>%
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



###################################### Posterior Checks ############################################

#Posterior predictive check
pp_check(diffusionmodel1_fit_gamma_19, 
         ndraws = 100, 
         type = 'dens_overlay_grouped', 
         group = 'collection_regionname' ) %>%
  .$data %>%
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
  scale_linewidth_manual('Data', values = c('FALSE' = 0.1, 
                                            'TRUE' = 1),
                         labels  = c(expression(italic('y')['rep']),
                                     expression(italic('y'))))+ 
  scale_colour_manual('Data',values = c('FALSE' = '#cbc9e2', 
                                        'TRUE' = '#54278f'),
                      labels  = c(expression(italic('y')['rep']),
                                  expression(italic('y')))) + 
  
  guides(alpha= 'none', 
         colour=guide_legend()) + facet_wrap(~group, labeller  = as_labeller(str_to_title)) + 
  scale_x_continuous(expression(paste('Predicted Weighted Diffusion Coefficient (',Km^{2}, ' ', year^{-1}, ')' )),
                     breaks = log1p(c(0, 10^(seq(from = 2, to = 8, by = 2)))),
                     labels = expression(0,  1%*%10^{2},  1%*%10^{4},1%*%10^{6},1%*%10^{8}),
                     limits = c(-0.01, log1p(10^8.5)),
                     expand = c(0.02,0.02))+
  scale_y_continuous('Density',expand = c(0,0)) +
  theme_classic() + 
  theme(strip.text = element_text(face = 'bold', size = 10),
        strip.background = element_blank(),
        legend.title = element_blank(),
        legend.position = 'bottom',
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 10),
        legend.text = element_text(size = 8))


ggsave('~/Downloads/flu_plots/diffusion_ppc.png',
       dpi = 360,
       device = 'png' ,
       height = 17,
       width = 20, 
       units = 'cm')

#Step 3B. Summarize posterior of variables
t <- get_variables(diffusionmodel1_fit_gamma_19)

beta_draws <- diffusionmodel1_fit_gamma_19 %>%
  gather_draws(., !!!syms(t)) %>%
  mutate(type = 'posterior') %>%
  filter(grepl('^b_', .variable)) %>%
  
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
                           .variable == 'b_shape_collection_regionnamecentral&northernamerica'~ "alpha['americas']",
                           
                           
                           .variable == 'sd_segment__Intercept'~ "sigma[gamma]",
                           .variable == 'sd_collection_regionname__shape_Intercept'~ "sigma[zeta]",
                           
                           .variable == 'r_segment[ha,Intercept] '~ "gamma['ha']",
                           .variable == "r_segment[mp,Intercept]"~ "gamma['mp']",
                           .variable == "r_segment[np,Intercept]"~ "gamma['np']",
                           .variable == "r_segment[ns,Intercept]"~ "gamma['ns']",
                           .variable == "r_segment[nx,Intercept]"~ "gamma['nx']",
                           .variable == "r_segment[pa,Intercept]"~ "gamma['pa']",
                           .variable == "r_segment[pb1,Intercept]"~ "gamma['pb1']",
                           .variable == "r_segment[pb2,Intercept]"~ "gamma['pb2']",
                           .variable == 'r_collection_regionname__shape[asia,Intercept]'~ "zeta['asia']",
                           .variable == 'r_collection_regionname__shape[africa,Intercept]'~ "zeta['africa']",
                           .variable == 'r_collection_regionname__shape[central.&.northern.america,Intercept]'~ "zeta['europe']",
                           .variable == 'collection_regionname__shape_central.&.northern.america'~ "zeta['central&northernamerica']",
                           
                           .variable == 'lprior' ~ 'prior',
                           .variable == 'lp__' ~ 'log~probability'))



 ggplot() + 
  geom_histogram(data = beta_draws, 
                 aes(x = .value,
                     y = after_stat(density)),
                 inherit.aes = F, 
                 binwidth = 0.1, 
                 fill = '#1b9e77') + 
  
  stat_function(fun = dnorm,
                data = tibble(label = "beta['asia']"),
                args = list(mean = 0, sd = 5),
                fill = '#d95f02',
                geom = 'area',
                alpha = 0.5) +
  
  stat_function(fun = dnorm,
                data = tibble(label = "beta['africa']"),
                args = list(mean = 0, sd = 5),
                fill = '#d95f02',
                geom = 'area',
                alpha = 0.5) +
  
  stat_function(fun = dnorm,
                data =tibble(label = "beta['europe']"),
                args = list(mean = 0, sd = 5),
                fill = '#d95f02',
                geom = 'area',
                alpha = 0.5) +
  
  stat_function(fun = dnorm,
                data = tibble(label = "beta['americas']"),
                args = list(mean = 0, sd = 5),
                fill = '#d95f02',
                geom = 'area',
                alpha = 0.5) +
  
  
  stat_function(fun = dnorm,
                data = tibble(label = "beta['step_anseriformes']"),
                args = list(mean = 0, sd = 1),
                fill = '#d95f02',
                geom = 'area',
                alpha = 0.5) +
  
  stat_function(fun = dnorm,
                data = tibble(label = "beta['persist_anseriformes']"),
                args = list(mean = 0, sd = 1),
                fill = '#d95f02',
                geom = 'area',
                alpha = 0.5) +
  
  
  stat_function(fun = dnorm,
                data = tibble(label = "beta['step_charadriiformes']"),
                args = list(mean = 0, sd = 1),
                fill = '#d95f02',
                geom = 'area',
                alpha = 0.5) +
  
  
  stat_function(fun = dnorm,
                data = tibble(label = "beta['persist_charadriiformes']"),
                args = list(mean = 0, sd = 1),
                fill = '#d95f02',
                geom = 'area',
                alpha = 0.5) +
  
  
  stat_function(fun = dnorm,
                data = tibble(label = "beta['step_hostjump']"),
                args = list(mean = 0, sd = 1),
                fill = '#d95f02',
                geom = 'area',
                alpha = 0.5) +
  
  
  stat_function(fun = dnorm,
                data = tibble(label = "beta['num_hostjump']"),
                args = list(mean = 0, sd = 1),
                fill = '#d95f02',
                geom = 'area',
                alpha = 0.5) +
  
  stat_function(fun = dnorm,
                data = tibble(label = "beta['africa-anser']"),
                args = list(mean = 0, sd = 1),
                fill = '#d95f02',
                geom = 'area',
                alpha = 0.5) +
  stat_function(fun = dnorm,
                data = tibble(label = "beta['europe-anser']"),
                args = list(mean = 0, sd = 1),
                fill = '#d95f02',
                geom = 'area',
                alpha = 0.5) +
  stat_function(fun = dnorm,
                data =tibble(label = "beta['americas-anser']"),
                args = list(mean = 0, sd = 1),
                fill = '#d95f02',
                geom = 'area',
                alpha = 0.5) +
  stat_function(fun = dnorm,
                data = tibble(label = "beta['africa-charad']"),
                args = list(mean = 0, sd = 1),
                fill = '#d95f02',
                geom = 'area',
                alpha = 0.5) +
  stat_function(fun = dnorm,
                data = tibble(label = "beta['europe-charad']"),
                args = list(mean = 0, sd = 1),
                fill = '#d95f02',
                geom = 'area',
                alpha = 0.5) +
  stat_function(fun = dnorm,
                data =tibble(label = "beta['americas-charad']"),
                args = list(mean = 0, sd = 1),
                fill = '#d95f02',
                geom = 'area',
                alpha = 0.5) +
  stat_function(fun = dnorm,
                data = tibble(label = "alpha['asia']",),
                args = list(mean = 0, sd = 5),
                fill = '#d95f02',
                geom = 'area',
                alpha = 0.5) +
  
  stat_function(fun = dnorm,
                data = tibble(label = "alpha['africa']",),
                args = list(mean = 0, sd = 5),
                fill = '#d95f02',
                geom = 'area',
                alpha = 0.5) +
  
  stat_function(fun = dnorm,
                data = tibble(label = "alpha['europe']",),
                args = list(mean = 0, sd = 5),
                fill = '#d95f02',
                geom = 'area',
                alpha = 0.5) +
  
  stat_function(fun = dnorm,
                data = tibble(label = "alpha['americas']",),
                args = list(mean = 0, sd = 5),
                fill = '#d95f02',
                geom = 'area',
                alpha = 0.5) +
  xlim(c(-15,15)) + 
  facet_wrap(~label, scales = 'free_y',  ncol = 3, labeller = label_parsed) +
  theme_minimal() +
   theme(axis.text = element_text(size = 8),
         axis.title = element_text(size = 10),
         legend.text = element_text(size = 8))
 

 


###################################### Residuals Checks ############################################
# Check Residuals using DHARMA 
# sample from the Posterior Predictive Distribution
preds <- posterior_predict(diffusionmodel1_fit_gamma_19, nsamples = 250, summary = FALSE)
preds <- t(preds)

res <- createDHARMa(
  simulatedResponse = t(posterior_predict(diffusionmodel1_fit_gamma_19)),
  observedResponse = diffusion_data$weighted_diff_coeff,
  fittedPredictedResponse = apply(t(posterior_epred(diffusionmodel1_fit_gamma_19)), 1, mean),
  integerResponse = FALSE)


# QQ plot
qq_data <- data.frame(
  sample = sort(res$scaledResiduals),
  theoretical = sort(ppoints(length(res$scaledResiduals)))
)

ggplot(qq_data, aes(sample = sample)) +
  stat_qq(distribution = stats::qunif) +
  stat_qq_line(distribution = stats::qunif) +
  scale_y_continuous('Observed')+
  scale_x_continuous('Expected')+
  theme_minimal() +
  theme(axis.text = element_text(size = 8),
        axis.title = element_text(size = 10),
        legend.text = element_text(size = 8))

ggsave('~/Downloads/flu_plots/diffusion_qq.jpeg',
       dpi = 360,
       height = 15,
       width = 15,
       units = 'cm')



############################################## WRITE ###############################################




############################################## END #################################################
####################################################################################################
####################################################################################################