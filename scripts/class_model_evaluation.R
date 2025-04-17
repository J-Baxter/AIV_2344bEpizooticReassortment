####################################################################################################
####################################################################################################
## Script name: Class Model Evaluation

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
library(ggmcmc)
library(ggdist)


# User functions
library(rlang)
make_probability_residuals = function(data, prediction, y, y_upper = NA, n = 1) {
  .prediction = enquo(prediction)
  .y = enquo(y)
  .y_upper = enquo(y_upper)
  
  if (eval_tidy(expr(is.factor(!!.prediction) && !is.ordered(!!.prediction)), data)) {
    data = mutate(data, !!.prediction := ordered(!!.prediction, levels = levels(!!.prediction)))
  }
  
  if (is.na(get_expr(.y_upper))) {
    #no y_upper provided, use y as y_upper
    data = summarise(data,
                     .p_lower = mean(!!.prediction < !!.y),
                     .p_upper = mean(!!.prediction <= !!.y),
                     .groups = "drop_last"
    )
  } else {
    #y_upper should be a vector, and if an entry in it is NA, use the entry from y
    data = summarise(data,
                     .p_lower = mean(!!.prediction < !!.y),
                     .p_upper = mean(!!.prediction <= ifelse(is.na(!!.y_upper), !!.y, !!.y_upper)),
                     .groups = "drop_last"
    )
  }
  
  data %>%
    mutate(
      .p_residual = map2(.p_lower, .p_upper, runif, n = !!n),
      .residual_draw = map(.p_residual, seq_along)
    ) %>%
    unnest(c(.p_residual, .residual_draw)) %>%
    mutate(.z_residual = qnorm(.p_residual))
}
######################################### DATA & MODEL #############################################
class_model <- class_model

diffusion_data <- 
  ###################################### MCMC Diagnostics #############################################

#MCMC chain convergence (Visual Inspection)
# Plot Chains
ggs(class_model) %>% 
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
                           .variable == "b_mumoderate_previous_classmoderate:collection_regionnamecentral&northernamerica"~ "beta[{'moderate'[md-americas]}]",
                           
                           .variable == "b_muminor_time_since_previous"~ "beta[{'minor'[time]}]",
                           .variable == "b_mumoderate_time_since_previous"~ "beta[{'moderate'[time]}]",
                           
                           .variable == "Intercept_muminor"~ "sigma['minor']",
                           .variable == "Intercept_mumoderate"~ "sigma['moderate']",
                           
                           .variable == 'lprior' ~ 'prior',
                           .variable == 'lp__' ~ 'log~probability')) %>% 
  drop_na(label) %>%
  ggplot(aes(x = .iteration,
             y = .value, 
             col = as.factor(.chain)))+
  geom_line(alpha = 0.8)+
  facet_wrap(~label,  ncol = 2, labeller = label_parsed, scales = 'free_y') +
  scale_colour_brewer(palette = 'GnBu', 'Chains') +
  scale_x_continuous('Iteration') + 
  scale_y_continuous('Parameter Value') + 
  theme_minimal() + 
  theme(legend.position = 'bottom',
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 10),
        legend.text = element_text(size = 8))

ggsave('~/Downloads/flu_plots/class_trace.jpeg',
       dpi = 360,
       height = 29,
       width = 21,
       units = 'cm')


# Plot ranked traces
# If chains are exploring the same space efficiently, the traces should be similar to one another 
# and largely overlapping.
as_draws_df(class_model) %>% 
  
  # Selecting only beta coefficients
  mcmc_rank_overlay(regex_pars = '^b_') %>%
  
  # Extract data
  .$data %>%
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
                           .variable == "b_mumoderate_previous_classmoderate:collection_regionnamecentral&northernamerica"~ "beta[{'moderate'[md-americas]}]",
                           
                           .variable == "b_muminor_time_since_previous"~ "beta[{'minor'[time]}]",
                           .variable == "b_mumoderate_time_since_previous"~ "beta[{'moderate'[time]}]",
                           
                           .variable == "Intercept_muminor"~ "sigma['minor']",
                           .variable == "Intercept_mumoderate"~ "sigma['moderate']",
                           
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
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 10),
        legend.text = element_text(size = 8))


#MCMC chain resolution (ESS)
model_resolution <- tidy(class_model, ess = T, rhat = T) %>% 
  dplyr::select(term, ess, rhat,term) %>%
  filter(!grepl('none', term)) %>%
  mutate(term = case_when(term == 'muminor_(Intercept)'~ "beta[{'minor'[0]}]",
                          term == 'mumoderate_(Intercept)'~ "beta[{'moderate'[0]}]",
                          
                          term == "muminor_previous_classmajor:collection_regionnameasia" ~ "beta[{'minor'[mj-asia]}]",
                          term == "muminor_previous_classminor:collection_regionnameasia" ~ "beta[{'minor'[mn-asia]}]",
                          term == "muminor_previous_classmoderate:collection_regionnameasia"~ "beta[{'minor'[md-asia]}]",
                          
                          term == "muminor_previous_classmajor:collection_regionnameeurope" ~ "beta[{'minor'[mj-europe]}]",
                          term == "muminor_previous_classminor:collection_regionnameeurope" ~ "beta[{'minor'[mn-europe]}]",
                          term == "muminor_previous_classmoderate:collection_regionnameeurope"~ "beta[{'minor'[md-europe]}]",
                          
                          term == "muminor_previous_classmajor:collection_regionnamecentral&northernamerica" ~ "beta[{'minor'[mj-americas]}]",
                          term == "muminor_previous_classminor:collection_regionnamecentral&northernamerica" ~ "beta[{'minor'[mn-americas]}]",
                          term == "muminor_previous_classmoderate:collection_regionnamecentral&northernamerica"~ "beta[{'minor'[md-americas]}]",
                          
                          term == "mumoderate_previous_classmajor:collection_regionnameasia" ~ "beta[{'moderate'[mj-asia]}]",
                          term == "mumoderate_previous_classminor:collection_regionnameasia" ~ "beta[{'moderate'[mn-asia]}]",
                          term == "mumoderate_previous_classmoderate:collection_regionnameasia"~ "beta[{'moderate'[md-asia]}]",
                          
                          term == "mumoderate_previous_classmajor:collection_regionnameeurope" ~ "beta[{'moderate'[mj-europe]}]",
                          term == "mumoderate_previous_classminor:collection_regionnameeurope" ~ "beta[{'moderate'[mn-europe]}]",
                          term == "mumoderate_previous_classmoderate:collection_regionnameeurope"~ "beta[{'moderate'[md-europe]}]",
                          
                          term == "mumoderate_previous_classmajor:collection_regionnamecentral&northernamerica" ~ "beta[{'moderate'[mj-americas]}]",
                          term == "mumoderate_previous_classminor:collection_regionnamecentral&northernamerica" ~ "beta[{'moderate'[mn-americas]}]",
                          term == "mumoderate_previous_classmoderate:collection_regionnamecentral&northernamerica"~ "beta[{'moderate'[md-americas]}]",
                          
                          term == "muminor_time_since_previous"~ "beta[{'minor'[time]}]",
                          term == "mumoderate_time_since_previous"~ "beta[{'moderate'[time]}]")) 



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
class_model %>% 
  mcmc_acf() %>% 
  .$data %>%
  as_tibble() %>%
  mutate(Parameter = case_when(Parameter == 'b_muminor_Intercept'~ "beta[{'minor'[0]}]",
                               Parameter == 'b_mumoderate_Intercept'~ "beta[{'moderate'[0]}]",
                               
                               Parameter == "b_muminor_previous_classmajor:collection_regionnameasia" ~ "beta[{'minor'[mj-asia]}]",
                               Parameter == "b_muminor_previous_classminor:collection_regionnameasia" ~ "beta[{'minor'[mn-asia]}]",
                               Parameter == "b_muminor_previous_classmoderate:collection_regionnameasia"~ "beta[{'minor'[md-asia]}]",
                               
                               Parameter == "b_muminor_previous_classmajor:collection_regionnameeurope" ~ "beta[{'minor'[mj-europe]}]",
                               Parameter == "b_muminor_previous_classminor:collection_regionnameeurope" ~ "beta[{'minor'[mn-europe]}]",
                               Parameter == "b_muminor_previous_classmoderate:collection_regionnameeurope"~ "beta[{'minor'[md-europe]}]",
                               
                               Parameter == "b_muminor_previous_classmajor:collection_regionnamecentral&northernamerica" ~ "beta[{'minor'[mj-americas]}]",
                               Parameter == "b_muminor_previous_classminor:collection_regionnamecentral&northernamerica" ~ "beta[{'minor'[mn-americas]}]",
                               Parameter == "b_muminor_previous_classmoderate:collection_regionnamecentral&northernamerica"~ "beta[{'minor'[md-americas]}]",
                               
                               Parameter == "b_mumoderate_previous_classmajor:collection_regionnameasia" ~ "beta[{'moderate'[mj-asia]}]",
                               Parameter == "b_mumoderate_previous_classminor:collection_regionnameasia" ~ "beta[{'moderate'[mn-asia]}]",
                               Parameter == "b_mumoderate_previous_classmoderate:collection_regionnameasia"~ "beta[{'moderate'[md-asia]}]",
                               
                               Parameter == "b_mumoderate_previous_classmajor:collection_regionnameeurope" ~ "beta[{'moderate'[mj-europe]}]",
                               Parameter == "b_mumoderate_previous_classminor:collection_regionnameeurope" ~ "beta[{'moderate'[mn-europe]}]",
                               Parameter == "b_mumoderate_previous_classmoderate:collection_regionnameeurope"~ "beta[{'moderate'[md-europe]}]",
                               
                               Parameter == "b_mumoderate_previous_classmajor:collection_regionnamecentral&northernamerica" ~ "beta[{'moderate'[mj-americas]}]",
                               Parameter == "b_mumoderate_previous_classminor:collection_regionnamecentral&northernamerica" ~ "beta[{'moderate'[mn-americas]}]",
                               Parameter == "b_mumoderate_previous_classmoderate:collection_regionnamecentral&northernamerica"~ "beta[{'moderate'[md-americas]}]",
                               
                               Parameter == "b_muminor_time_since_previous"~ "beta[{'minor'[time]}]",
                               Parameter == "b_mumoderate_time_since_previous"~ "beta[{'moderate'[time]}]",
                               
                               Parameter == "Intercept_muminor"~ "sigma['minor']",
                               Parameter == "Intercept_mumoderate"~ "sigma['moderate']",
                               
                               Parameter == 'lprior' ~ 'prior',
                               Parameter == 'lp__' ~ 'log~probability')) %>%
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

ggsave('~/Downloads/flu_plots/class_autocorrelation.jpeg',
       dpi = 360,
       height = 29,
       width = 21,
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
pp_check(class_model, 
         ndraws = 100, 
         type = 'bars_grouped', 
         group = 'collection_regionname' ) %>%
  .$data %>%
  mutate(x = case_when(x == 1 ~ 'Major',
                       x == 2 ~ 'Minor',
                       x == 3 ~ 'Moderate') %>%
           factor(levels = c('Major', 'Moderate', 'Minor'))) %>%
  
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
                     colour = '#54278f') +

facet_wrap(~group, labeller  = as_labeller(str_to_title)) + 
  scale_x_discrete('Reassortant Class')+
  scale_y_continuous('Count', expand = c(0,0)) +
  theme_classic() + 
  theme(strip.text = element_text(face = 'bold', size = 10),
        strip.background = element_blank(),
        legend.title = element_blank(),
        legend.position = 'bottom',
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 10),
        legend.text = element_text(size = 8))


ggsave('~/Downloads/flu_plots/class_ppc.jpeg',
       dpi = 360,
       device = 'jpeg' ,
       height = 12,
       width = 17, 
       units = 'cm')




#Step 3B. Summarize posterior of variables
t <- get_variables(class_model)

beta_draws <- class_model %>%
  gather_draws(., !!!syms(t)) %>%
  mutate(type = 'posterior') %>%
  filter(grepl('^b_', .variable)) %>%
  filter(!grepl('none', .variable)) %>%
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
                           .variable == "b_mumoderate_previous_classmoderate:collection_regionnamecentral&northernamerica"~ "beta[{'moderate'[md-americas]}]",
                           
                           .variable == "b_muminor_time_since_previous"~ "beta[{'minor'[time]}]",
                           .variable == "b_mumoderate_time_since_previous"~ "beta[{'moderate'[time]}]",
                           
                           .variable == "Intercept_muminor"~ "sigma['minor']",
                           .variable == "Intercept_mumoderate"~ "sigma['moderate']",
                           
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
                #data = tibble(label = "beta[{'minor'[0]}]"),
                args = list(mean = 0, sd = 5),
                fill = '#d95f02',
                geom = 'area',
                alpha = 0.5) +
  scale_y_continuous('Probability Density') + 
  scale_x_continuous('Parameter Value', limits = c(-15,15)) + 
  facet_wrap(~label, scales = 'free_y',  ncol = 3, labeller = label_parsed) +
  theme_minimal() +
  theme(axis.text = element_text(size = 8),
        axis.title = element_text(size = 10),
        legend.text = element_text(size = 8))

ggsave('~/Downloads/flu_plots/class_identifiability.jpeg',
       dpi = 360,
       height = 29,
       width = 20,
       units = 'cm')


###################################### Residuals Checks ############################################
# Check Randommised Residuals 
# sample from the Posterior Predictive Distribution
class_data %>%
  add_predicted_draws(class_model) %>%
  make_probability_residuals(.prediction, class, n = 1) %>%
  ggplot(aes(sample = .p_residual)) +
  geom_qq(distribution = qunif) +
  geom_abline() +
  scale_y_continuous('Observed')+
  scale_x_continuous('Expected')+
  theme_minimal() +
  theme(axis.text = element_text(size = 8),
        axis.title = element_text(size = 10),
        legend.text = element_text(size = 8))

ggsave('~/Downloads/flu_plots/class_qq.jpeg',
       dpi = 360,
       height = 15,
       width = 15,
       units = 'cm')


############################################## END #################################################
####################################################################################################
####################################################################################################