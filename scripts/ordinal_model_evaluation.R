####################################################################################################
####################################################################################################
## Script name: Ordinal Model Evaluation
##
## Purpose of script:
##
## Date created: 2025-05-17
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
fitted_model <- ordinal_model

diffusion_data <- ordinal_data
###################################### MCMC Diagnostics #############################################

#MCMC chain convergence (Visual Inspection)
# Plot Chains
ggs(ordinal_model) %>% 
  from_ggmcmc_names() %>%
  filter(.iteration > 400) %>% 
  mutate(label = case_when(.variable == 'b_Intercept[1]'~ "tau['minor']",
                           .variable == 'b_Intercept[2]'~ "tau['moderate']", 
                           
                           .variable == 'b_parent_classminor'~ "beta['minor']",
                           .variable == "b_parent_classmoderate"~ "beta['moderate']",
                           
                           .variable == 'b_cluster_regionasia'~ "beta['asia']",
                           .variable == 'b_cluster_regioncentral&northernamerica'~ "beta['americas']",
                           .variable == 'b_cluster_regioneurope'~ "beta['europe']",
                           .variable == 'b_segments_changed'~ "beta['segments-changed']",
                           #.variable == 'bs_stime_since_last_major_1'~ "beta['smooth']",
                           .variable == 's_stime_since_last_major_1[1]'~"gamma['1']",
                          .variable == 's_stime_since_last_major_1[2]'~"gamma['2']",
                          .variable == 's_stime_since_last_major_1[3]'~"gamma['3']",
                          .variable == 's_stime_since_last_major_1[4]'~"gamma['4']",
                          .variable == 's_stime_since_last_major_1[5]'~"gamma['5']",
                          .variable == 's_stime_since_last_major_1[6]'~"gamma['6']",
                          .variable == 's_stime_since_last_major_1[7]'~"gamma['7']",
                          .variable == 's_stime_since_last_major_1[8]'~"gamma['8']",
                           .variable == 'sds_stime_since_last_major_1'~ "sigma['europe-anser']",
                          
                           .variable == 'sd_collection_year__Intercept'~ "sigma[year]",
                           
                           .variable == 'lprior' ~ 'prior')) %>%
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
        axis.text = element_text(size = 7),
        axis.title = element_text(size = 10),
        legend.text = element_text(size = 7))

ggsave('~/Downloads/flu_plots/ordinal_trace.jpeg',
       dpi = 360,
       height = 25,
       width = 17,
       units = 'cm')

# Plot ranked traces
# If chains are exploring the same space efficiently, the traces should be similar to one another 
# and largely overlapping.
as_draws_df(ordinal_model) %>% 
  
  # Selecting only beta coefficients
  mcmc_rank_overlay() %>%
  
  # Extract data
  .$data %>%
  rename(.variable = parameter) %>%
  mutate(label = case_when(.variable == 'b_Intercept[1]'~ "tau['minor']",
                           .variable == 'b_Intercept[2]'~ "tau['moderate']", 
                           
                           .variable == 'b_parent_classminor'~ "beta['minor']",
                           .variable == "b_parent_classmoderate"~ "beta['moderate']",
                           
                           .variable == 'b_cluster_regionasia'~ "beta['asia']",
                           .variable == 'b_cluster_regioncentral&northernamerica'~ "beta['americas']",
                           .variable == 'b_cluster_regioneurope'~ "beta['europe']",
                           .variable == 'b_segments_changed'~ "beta['segments-changed']",
                           .variable == 'bs_stime_since_last_major_1'~ "beta['smooth']",
                           .variable == 's_stime_since_last_major_1[1]'~"gamma['1']",
                           .variable == 's_stime_since_last_major_1[2]'~"gamma['2']",
                           .variable == 's_stime_since_last_major_1[3]'~"gamma['3']",
                           .variable == 's_stime_since_last_major_1[4]'~"gamma['4']",
                           .variable == 's_stime_since_last_major_1[5]'~"gamma['5']",
                           .variable == 's_stime_since_last_major_1[6]'~"gamma['6']",
                           .variable == 's_stime_since_last_major_1[7]'~"gamma['7']",
                           .variable == 's_stime_since_last_major_1[8]'~"gamma['8']",
                           .variable == 'sds_stime_since_last_major_1'~ "sigma['europe-anser']",
                           
                           .variable == 'sd_collection_year__Intercept'~ "sigma[year]",
                           
                           .variable == 'lprior' ~ 'prior')) %>%
  drop_na() %>%
  
  # Plot 
  ggplot() + 
  geom_step(aes(x = bin_start, 
                y = n,
                colour = chain),
            linewidth = 0.8) + 
  #scale_y_continuous(expand = c(0,0), limits = c(150, 225)) +
  scale_x_continuous(expand = c(0,0)) + 
  facet_wrap(~label,  ncol = 3, labeller = label_parsed) +
  scale_colour_brewer(palette = 'GnBu', 'Chains') +
  theme_minimal() + 
  theme(legend.position = 'bottom',
        axis.title = element_blank())


#MCMC chain resolution (ESS)
model_resolution <- list(tidy(ordinal_model, ess = T, rhat = T, effects = 'ran_vals') %>%
                           dplyr::select( ess, rhat, group, level) %>%
                           unite(term, group, level, sep = '_'),
                         tidy(ordinal_model, ess = T, rhat = T) %>%  dplyr::select(term, ess, rhat,term) ) %>%
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
ordinal_model %>% 
  mcmc_acf() %>% 
  .$data %>%
  as_tibble() %>%
  mutate(Parameter = case_when(Parameter == 'b_Intercept[1]'~ "tau['minor']",
                               Parameter == 'b_Intercept[2]'~ "tau['moderate']", 
                               
                               Parameter == 'b_parent_classminor'~ "beta['minor']",
                               Parameter == "b_parent_classmoderate"~ "beta['moderate']",
                               
                               Parameter == 'b_cluster_regionasia'~ "beta['asia']",
                               Parameter == 'b_cluster_regioncentral&northernamerica'~ "beta['americas']",
                               Parameter == 'b_cluster_regioneurope'~ "beta['europe']",
                               Parameter == 'b_segments_changed'~ "beta['segments-changed']",
                               #Parameter == 'bs_stime_since_last_major_1'~ "beta['smooth']",
                               Parameter == 's_stime_since_last_major_1[1]'~"gamma['1']",
                               Parameter == 's_stime_since_last_major_1[2]'~"gamma['2']",
                               Parameter == 's_stime_since_last_major_1[3]'~"gamma['3']",
                               Parameter == 's_stime_since_last_major_1[4]'~"gamma['4']",
                               Parameter == 's_stime_since_last_major_1[5]'~"gamma['5']",
                               Parameter == 's_stime_since_last_major_1[6]'~"gamma['6']",
                               Parameter == 's_stime_since_last_major_1[7]'~"gamma['7']",
                               Parameter == 's_stime_since_last_major_1[8]'~"gamma['8']",
                               Parameter == 'sds_stime_since_last_major_1'~ "sigma['europe-anser']",
                               
                               Parameter == 'sd_collection_year__Intercept'~ "sigma[year]",
                               
                               Parameter == 'lprior' ~ 'prior')) %>%
  drop_na(Parameter) %>%
  ggplot(aes(y = AC, x = Lag, colour = as.factor(Chain))) + 
  geom_path() + 
  facet_wrap(~Parameter, labeller = label_parsed, ncol = 4) + 
  theme_minimal()  + 
  scale_colour_brewer('Chain') +
  theme(legend.position = 'bottom',
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 10),
        legend.text = element_text(size = 8))

ggsave('~/Downloads/flu_plots/ordinal_autocorrelation.jpeg',
       dpi = 360,
       height = 20,
       width = 16,
       units = 'cm')


# Check Ratio of Effective Population Size to Total Sample Size 
# values <0.1 should raise concerns about autocorrelation
neff_ratio(ordinal_model) %>% 
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
pp_check(ordinal_model, 
         ndraws = 100, 
         type = 'bars_grouped', 
         group = 'cluster_region' ) %>%
  .$data %>%
  mutate(x = case_when(x == 1 ~ 'minor',
                       x == 2 ~ 'moderate',
                       x == 3 ~ 'major')) %>%
  #filter(!is_y) %>%
  ggplot() +
  geom_bar(aes(x = x,
               y = y_obs),
           stat = 'identity',
           fill = '#cbc9e2',
           colour = '#cbc9e2',
           alpha = 0.7) + 
  geom_pointinterval(aes(x = x,
                         y=m,
                         ymin = l,
                         ymax = h), orientation = 'x',
                     colour = '#54278f')+ 
  facet_wrap(~group, labeller  = as_labeller(str_to_title)) +
scale_x_discrete('Reassortants Class', labels = str_to_title)+
  scale_y_continuous('Count') +
  theme_classic() + 
  theme(strip.text = element_text(face = 'bold', size = 9),
        strip.background = element_blank(),
        legend.title = element_blank(),
        legend.position = 'bottom',
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 9),
        legend.text = element_text(size = 8))



ggsave('~/Downloads/flu_plots/ordinal_ppc.jpeg',
       dpi = 360,
       device = 'jpeg' ,
       height = 17,
       width = 17, 
       units = 'cm')

#Step 3B. Summarize posterior of variables
t <- get_variables(ordinal_model)

beta_draws <- ordinal_model %>%
  gather_draws(., !!!syms(t)) %>%
  mutate(type = 'posterior') %>%
  #filter(grepl('^b_', .variable)) %>%
  
  mutate(label = case_when(.variable == 'b_Intercept[1]'~ "tau['minor']",
                           .variable == 'b_Intercept[2]'~ "tau['moderate']", 
                           
                           .variable == 'b_parent_classminor'~ "beta['minor']",
                           .variable == "b_parent_classmoderate"~ "beta['moderate']",
                           
                           .variable == 'b_cluster_regionasia'~ "beta['asia']",
                           .variable == 'b_cluster_regioncentral&northernamerica'~ "beta['americas']",
                           .variable == 'b_cluster_regioneurope'~ "beta['europe']",
                           .variable == 'b_segments_changed'~ "beta['segments-changed']",
                           #.variable == 'bs_stime_since_last_major_1'~ "beta['smooth']",
                           .variable == 's_stime_since_last_major_1[1]'~"gamma['1']",
                           .variable == 's_stime_since_last_major_1[2]'~"gamma['2']",
                           .variable == 's_stime_since_last_major_1[3]'~"gamma['3']",
                           .variable == 's_stime_since_last_major_1[4]'~"gamma['4']",
                           .variable == 's_stime_since_last_major_1[5]'~"gamma['5']",
                           .variable == 's_stime_since_last_major_1[6]'~"gamma['6']",
                           .variable == 's_stime_since_last_major_1[7]'~"gamma['7']",
                           .variable == 's_stime_since_last_major_1[8]'~"gamma['8']",
                           .variable == 'sds_stime_since_last_major_1'~ "sigma['smooth']",
                           
                           .variable == 'sd_collection_year__Intercept'~ "sigma[year]")) %>%
  drop_na()



ggplot() + 
  geom_histogram(data = beta_draws %>% drop_na(label), 
                 aes(x = .value,
                     y = after_stat(density)),
                 inherit.aes = F, 
                 bins = 70,
                 fill = '#1b9e77') + 
  
  stat_function(fun = dstudent_t,
                data = tibble(label = "tau['minor']"),
                args = list(df = 3, mu = 0, sigma = 2.5),
                fill = '#d95f02',
                geom = 'area',
                alpha = 0.5) +
  
  stat_function(fun = dstudent_t,
                data = tibble(label = "tau['moderate']"),
                args = list(df = 3, mu = 0, sigma = 2.5),
                fill = '#d95f02',
                geom = 'area',
                alpha = 0.5) +
  
  stat_function(fun = dnorm,
                data =tibble(label = "beta['minor']"),
                args = list(mean = 0, sd = 2),,
                fill = '#d95f02',
                geom = 'area',
                alpha = 0.5) +
  
  stat_function(fun = dnorm,
                data = tibble(label = "beta['moderate']"),
                args = list(mean = 0, sd = 2),
                fill = '#d95f02',
                geom = 'area',
                alpha = 0.5) +   
  
  stat_function(fun = dnorm,
                data = tibble(label = "beta['asia']"),
                args = list(mean = 0, sd = 2),
                fill = '#d95f02',
                geom = 'area',
                alpha = 0.5) +
  
  stat_function(fun = dnorm,
                data = tibble(label = "beta['americas']"),
                args = list(mean = 0, sd = 2),
                fill = '#d95f02',
                geom = 'area',
                alpha = 0.5) +
  
  stat_function(fun = dnorm,
                data = tibble(label = "beta['europe']"),
                args = list(mean = 0, sd = 2),
                fill = '#d95f02',
                geom = 'area',
                alpha = 0.5) +
  
  
  stat_function(fun = dnorm,
                data = tibble(label = "beta['segments-changed']"),
                args = list(mean = 0, sd = 2),
                fill = '#d95f02',
                geom = 'area',
                alpha = 0.5) +
  
  stat_function(fun = dnorm,
                data = tibble(label = "gamma['1']"),
                args = list(mean = 0, sd = 2),
                fill = '#d95f02',
                geom = 'area',
                alpha = 0.5) +
  
  stat_function(fun = dnorm,
                data = tibble(label = "gamma['2']"),
                args = list(mean = 0, sd = 2),
                fill = '#d95f02',
                geom = 'area',
                alpha = 0.5) +
  
  stat_function(fun = dnorm,
                data = tibble(label = "gamma['3']"),
                args = list(mean = 0, sd = 2),
                fill = '#d95f02',
                geom = 'area',
                alpha = 0.5) +
  
  stat_function(fun = dnorm,
                data = tibble(label = "gamma['4']"),
                args = list(mean = 0, sd = 2),
                fill = '#d95f02',
                geom = 'area',
                alpha = 0.5) +
  
  
  stat_function(fun = dnorm,
                data = tibble(label = "gamma['5']"),
                args = list(mean = 0, sd = 1),
                fill = '#d95f02',
                geom = 'area',
                alpha = 0.5) +
  
  stat_function(fun = dnorm,
                data = tibble(label = "gamma['6']"),
                args = list(mean = 0, sd = 2),
                fill = '#d95f02',
                geom = 'area',
                alpha = 0.5) +
  
  stat_function(fun = dnorm,
                data = tibble(label = "gamma['7']"),
                args = list(mean = 0, sd = 2),
                fill = '#d95f02',
                geom = 'area',
                alpha = 0.5) +
  
  stat_function(fun = dnorm,
                data = tibble(label = "gamma['8']"),
                args = list(mean = 0, sd = 2),
                fill = '#d95f02',
                geom = 'area',
                alpha = 0.5) +
  
  
  stat_function(fun = dexp,
                data = tibble(label = "sigma['smooth']"),
                args = list(rate = 0.5),
                fill = '#d95f02',
                geom = 'area',
                alpha = 0.5) +
  
  # stat_function(fun = dstudent_t,
  #data = tibble(label = "rho[reassortant]"),
  # args = list(df = 3, mu= 0, sigma =0.25),
  #fill = '#d95f02',
  #geom = 'area',
  #alpha = 0.5) +
  
  scale_y_continuous('Probability Density') + 
  scale_x_continuous('Parameter Value') + 
  facet_wrap(~label, scales = 'free',  ncol = 3, labeller = label_parsed) +
  theme_minimal() +
  theme(axis.text = element_text(size = 8),
        axis.title = element_text(size = 10),
        legend.text = element_text(size = 8))

ggsave('~/Downloads/flu_plots/ordinal_identifiability.jpeg',
       dpi = 360,
       height = 29,
       width = 20,
       units = 'cm')



###################################### Residuals Checks ############################################
# Check Residuals using DHARMA 
# sample from the Posterior Predictive Distribution
simulatedResponse = ordinal_model %>% 
  predict(type="response") %>% 
  apply(1, \(x) sample(1:length(x), prob = x, size=1000, replace = TRUE)) %>% 
  t()
model.check <-

preds <- posterior_predict(ordinal_model, nsamples = 250, summary = FALSE)
preds <- t(preds)

res <- createDHARMa(
  simulatedResponse = simulatedResponse,
  observedResponse = as.integer(class_data %>% drop_na() %>% pull(cluster_class)),
  fittedPredictedResponse = ordinal_model %>% predict(), # the linear predictor
  integerResponse = !all(as.numeric(simulatedResponse)%%1!=0) # checks if prediction consists of integers
)


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

ggsave('~/Downloads/flu_plots/ordinal_qq.jpeg',
       dpi = 360,
       height = 12,
       width = 12,
       units = 'cm')



############################################## WRITE ###############################################




############################################## END #################################################
####################################################################################################
###################################################################################################