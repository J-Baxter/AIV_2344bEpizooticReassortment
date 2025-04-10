####################################################################################################
####################################################################################################
## Script name:
##
## Purpose of script:
##
## Date created: 2025-03-28
##
##
########################################## SYSTEM OPTIONS ##########################################
options(scipen = 6, digits = 7) 
memory.limit(30000000) 

  
########################################## DEPENDENCIES ############################################
# Packages
library(tidyverse)
library(magrittr)


# User functions

# User functions
scientific_10 <- function(x) {
  parse(text=gsub("e\\+*", " %*% 10^", scales::scientific_format()(x)))
  }

############################################## DATA ################################################
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
mcmc_class <- ggs(basic_model) # Warning message In custom.sort(D$Parameter) : NAs introduced by coercion

mcmc_class %>% 
  from_ggmcmc_names() %>%
  filter(.iteration > 400) %>%
  mutate(label = case_when(.variable == 'b_muminor_Intercept'~ "beta[{'minor'[0]}]",
                           .variable == 'b_mumoderate_Intercept'~ "beta[{'moderate'[0]}]",
                           
                           .variable == "b_muminor_previous_classmajor:collection_regionnameasia" ~ "beta[{'minor'[mj-asia]}]",
                           .variable == "b_muminor_previous_classminor:collection_regionnameasia" ~ "beta[{'minor'[mn-asia]}]",
                           .variable == "b_muminor_previous_classmoderate:collection_regionnameasia"~ "beta[{'minor'[md-asia]}]",
                           
                           .variable == "b_muminor_previous_classmajor:collection_regionnameeurope" ~ "beta[{'minor'[mj-europe]}]",
                           .variable == "b_muminor_previous_classminor:collection_regionnameeurope" ~ "beta[{'minor'[mn-europe]}]",
                           .variable == "b_muminor_previous_classmoderate:collection_regionnameeurope"~ "beta[{'minor'[md-europe]}]",
                           
                           .variable == "b_muminor_previous_classmajor:collection_regionnamecentral&northernamerica" ~ "beta[{'minor'[mj-americas]}]",
                           .variable == "b_muminor_previous_classminor:collection_regionnamecentral&northernamerica" ~ "beta[{'minor'[mn-americas]}]",
                           .variable == "b_muminor_previous_classmoderate:collection_regionnamecentral&northernamerica"~ "beta[{'minor'[md-americas]}]",
                           
                           .variable == "b_mumoderate_previous_classmajor:collection_regionnameasia" ~ "beta[{'moderate'[mj-asia]}]",
                           .variable == "b_mumoderate_previous_classminor:collection_regionnameasia" ~ "beta[{'moderate'[mn-asia]}]",
                           .variable == "b_mumoderate_previous_classmoderate:collection_regionnameasia"~ "beta[{'moderate'[md-asia]}]",
                           
                           .variable == "b_mumoderate_previous_classmajor:collection_regionnameeurope" ~ "beta[{'moderate'[mj-europe]}]",
                           .variable == "b_mumoderate_previous_classminor:collection_regionnameeurope" ~ "beta[{'moderate'[mn-europe]}]",
                           .variable == "b_mumoderate_previous_classmoderate:collection_regionnameeurope"~ "beta[{'moderate'[md-europe]}]",
                           
                           .variable == "b_mumoderate_previous_classmajor:collection_regionnamecentral&northernamerica" ~ "beta[{'moderate'[mj-americas]}]",
                           .variable == "b_mumoderate_previous_classminor:collection_regionnamecentral&northernamerica" ~ "beta[{'moderate'[mn-americas]}]",
                           .variable == "b_mumoderate_previous_classmoderate:collection_regionnamecentral&northernamerica"~ "beta[{'moderate'[md-americas]}]")) %>%
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
  theme_minimal(base_size = 8) + 
  theme(legend.position = 'bottom')


### Ranked Traces ###
trank <- as_draws_df(basic_model)  %>%  
  mcmc_rank_overlay(regex_pars = '^b_') 

trank$data %>%
  rename(.variable = parameter) %>%
  mutate(label = case_when(.variable == 'b_muminor_Intercept'~ "beta[{'minor'[0]}]",
                           .variable == 'b_mumoderate_Intercept'~ "beta[{'moderate'[0]}]",
                           
                           .variable == "b_muminor_previous_classmajor:collection_regionnameasia" ~ "beta[{'minor'[mj-asia]}]",
                           .variable == "b_muminor_previous_classminor:collection_regionnameasia" ~ "beta[{'minor'[mn-asia]}]",
                           .variable == "b_muminor_previous_classmoderate:collection_regionnameasia"~ "beta[{'minor'[md-asia]}]",
                           
                           .variable == "b_muminor_previous_classmajor:collection_regionnameeurope" ~ "beta[{'minor'[mj-europe]}]",
                           .variable == "b_muminor_previous_classminor:collection_regionnameeurope" ~ "beta[{'minor'[mn-europe]}]",
                           .variable == "b_muminor_previous_classmoderate:collection_regionnameeurope"~ "beta[{'minor'[md-europe]}]",
                           
                           .variable == "b_muminor_previous_classmajor:collection_regionnamecentral&northernamerica" ~ "beta[{'minor'[mj-americas]}]",
                           .variable == "b_muminor_previous_classminor:collection_regionnamecentral&northernamerica" ~ "beta[{'minor'[mn-americas]}]",
                           .variable == "b_muminor_previous_classmoderate:collection_regionnamecentral&northernamerica"~ "beta[{'minor'[md-americas]}]",
                           
                           .variable == "b_mumoderate_previous_classmajor:collection_regionnameasia" ~ "beta[{'moderate'[mj-asia]}]",
                           .variable == "b_mumoderate_previous_classminor:collection_regionnameasia" ~ "beta[{'moderate'[mn-asia]}]",
                           .variable == "b_mumoderate_previous_classmoderate:collection_regionnameasia"~ "beta[{'moderate'[md-asia]}]",
                           
                           .variable == "b_mumoderate_previous_classmajor:collection_regionnameeurope" ~ "beta[{'moderate'[mj-europe]}]",
                           .variable == "b_mumoderate_previous_classminor:collection_regionnameeurope" ~ "beta[{'moderate'[mn-europe]}]",
                           .variable == "b_mumoderate_previous_classmoderate:collection_regionnameeurope"~ "beta[{'moderate'[md-europe]}]",
                           
                           .variable == "b_mumoderate_previous_classmajor:collection_regionnamecentral&northernamerica" ~ "beta[{'moderate'[mj-americas]}]",
                           .variable == "b_mumoderate_previous_classminor:collection_regionnamecentral&northernamerica" ~ "beta[{'moderate'[mn-americas]}]",
                           .variable == "b_mumoderate_previous_classmoderate:collection_regionnamecentral&northernamerica"~ "beta[{'moderate'[md-americas]}]")) %>%
  drop_na(label) %>%
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


#### Check Ratio of Effective Population Size to Total Sample Size #### 
# values <0.1 should raise concerns about autocorrelation
neff_ratio(basic_model) %>% as_tibble(rownames = 'param') %>%
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


# Autocorrelation plot
basic_model %>% 
  mcmc_acf() %>% 
  .$data %>%
  as_tibble() %>%
  ggplot(aes(y = AC, x = Lag, colour = as.factor(Chain))) + 
  geom_path() + 
  facet_wrap(~Parameter) + 
  theme_classic() 



# Compare prior & posterior parameter distributions
t_class<- get_variables(basic_model)

beta_draws <- basic_model %>%
  gather_draws(., !!!syms(t_class)) %>%
  mutate(type = 'posterior') %>%
  filter(grepl('^b_', .variable)) %>%
  
  mutate(label = case_when(.variable == 'b_muminor_Intercept'~ "beta[{'minor'[0]}]",
                           .variable == 'b_mumoderate_Intercept'~ "beta[{'moderate'[0]}]",
                           
                           .variable == "b_muminor_previous_classmajor:collection_regionnameasia" ~ "beta[{'minor'[mj-asia]}]",
                           .variable == "b_muminor_previous_classminor:collection_regionnameasia" ~ "beta[{'minor'[mn-asia]}]",
                           .variable == "b_muminor_previous_classmoderate:collection_regionnameasia"~ "beta[{'minor'[md-asia]}]",
                           
                           .variable == "b_muminor_previous_classmajor:collection_regionnameeurope" ~ "beta[{'minor'[mj-europe]}]",
                           .variable == "b_muminor_previous_classminor:collection_regionnameeurope" ~ "beta[{'minor'[mn-europe]}]",
                           .variable == "b_muminor_previous_classmoderate:collection_regionnameeurope"~ "beta[{'minor'[md-europe]}]",
                           
                           .variable == "b_muminor_previous_classmajor:collection_regionnamecentral&northernamerica" ~ "beta[{'minor'[mj-americas]}]",
                           .variable == "b_muminor_previous_classminor:collection_regionnamecentral&northernamerica" ~ "beta[{'minor'[mn-americas]}]",
                           .variable == "b_muminor_previous_classmoderate:collection_regionnamecentral&northernamerica"~ "beta[{'minor'[md-americas]}]",
                           
                           .variable == "b_mumoderate_previous_classmajor:collection_regionnameasia" ~ "beta[{'moderate'[mj-asia]}]",
                           .variable == "b_mumoderate_previous_classminor:collection_regionnameasia" ~ "beta[{'moderate'[mn-asia]}]",
                           .variable == "b_mumoderate_previous_classmoderate:collection_regionnameasia"~ "beta[{'moderate'[md-asia]}]",
                           
                           .variable == "b_mumoderate_previous_classmajor:collection_regionnameeurope" ~ "beta[{'moderate'[mj-europe]}]",
                           .variable == "b_mumoderate_previous_classminor:collection_regionnameeurope" ~ "beta[{'moderate'[mn-europe]}]",
                           .variable == "b_mumoderate_previous_classmoderate:collection_regionnameeurope"~ "beta[{'moderate'[md-europe]}]",
                           
                           .variable == "b_mumoderate_previous_classmajor:collection_regionnamecentral&northernamerica" ~ "beta[{'moderate'[mj-americas]}]",
                           .variable == "b_mumoderate_previous_classminor:collection_regionnamecentral&northernamerica" ~ "beta[{'moderate'[mn-americas]}]",
                           .variable == "b_mumoderate_previous_classmoderate:collection_regionnamecentral&northernamerica"~ "beta[{'moderate'[md-americas]}]"))

 beta_draws %>% 
  drop_na(label) %>%
  ggplot() + 
  geom_histogram(aes(x = .value,
                     y = after_stat(density)),
                 inherit.aes = F, 
                 binwidth = 0.1, 
                 fill = '#1b9e77') + 
  
  stat_function(fun = dnorm,
                args = list(mean = 0, sd = 5),
                fill = '#d95f02',
                geom = 'area',
                alpha = 0.5) +
  
  xlim(c(-15,15)) + 
  facet_wrap(~label, scales = 'free_y',  ncol = 3, labeller = label_parsed) +
  theme_minimal()
 

#### Check Residuals using DHARMA #### # currently not working
# sample from the Posterior Predictive Distribution
preds <- posterior_predict(basic_model, nsamples = 250, summary = FALSE)
preds <- t(preds)

res <- createDHARMa(
  simulatedResponse = t(posterior_predict(basic_model)),
  observedResponse = class_data$class,
  fittedPredictedResponse = apply(t(posterior_epred(basic_model)), 1, mean),
  integerResponse = FALSE)

plot(res, quantreg = FALSE)
############################################## WRITE ###############################################




############################################## END #################################################
####################################################################################################
####################################################################################################