####################################################################################################
####################################################################################################
## Script name: Number and Class of Reassortants Model
##
## Purpose of script: to determine whether reassortants emerging from some regions are more
## likely and more severe than those from others. Variables of interest interest include: region of
## origin, persistence time (persistence time increases uncertainty), number of sequences per region
##
## This script implements a workflow for a multivariate negative binomial model. 
##
## Date created: 2024-xx-xx
##
##
########################################## SYSTEM OPTIONS ##########################################
#options(scipen = 6, digits = 7) 
memory.limit(30000000) 


########################################## DEPENDENCIES ############################################
# Packages
library(tidyverse)
library(magrittr)
library(cmdstanr)
library(posterior)
library(broom)
library(broom.mixed)
library(tidybayes)
library(bayesplot)
library(emmeans)
library(marginaleffects)
library(magrittr)
library(ggmcmc)


# User functions

FormatContinent <- function(dataframe){
  dataframe %<>%
    mutate(collection_regionname = case_when(grepl('europe', 
                                                   collection_regionname, 
                                                   ignore.case = TRUE) ~ 'europe',
                                             
                                             grepl('africa',
                                                   collection_regionname, 
                                                   ignore.case = TRUE) ~ 'africa',
                                             
                                             grepl('asia', 
                                                   collection_regionname, 
                                                   ignore.case = TRUE) ~ 'asia',
                                             
                                             grepl('(central|northern) america', 
                                                   collection_regionname, 
                                                   ignore.case = TRUE) ~ 'central & northern america',
                                             
                                             .default = NA_character_ ))
  
  return(dataframe)
}

############################################## DATA ################################################
combined_data <- read_csv('./2024Aug18/treedata_extractions/2024-09-20_combined_data.csv')

data <- read_csv('countmodeldata_2025Apr23.csv')

############################################## MAIN ################################################
# Data pre-processing
data_processed_2 <- data %>%
  
  # remove reassortant classes
  dplyr::select(-c(ends_with('class'), time_since_last_dominant)) %>%
  
  # scaling
  mutate(across(c(collection_dec, collection_year), .fns = ~ subtract(.x, 2015))) %>%
  
  # set default NA values
  replace_na(list(woah_cases = 0, 
                  woah_susceptibles = 0,
                  woah_deaths = 0,
                  n_sequences = 0,
                  n_reassortants = 0,
                  time_since_last_dominant = 0, # this may be quite a strong assumption (accidentally)
                  minor = 0,
                  major = 0,
                  dominant = 0)) %>%
  mutate(across(starts_with('woah'), ~.x/6, .names = "{.col}_monthly")) %>%
  
  mutate(across(ends_with('_monthly'), ~log1p(.x), .names = "{.col}_log1p")) %>%
  mutate(n_sequences_log1p = log1p(n_sequences)) %>%
  #rename_with(~gsub('_', '-' ,.x)) %>%
  
  filter(collection_regionname != 'south america')


# Compile Stan Model
numbers_mod <- cmdstan_model('./scripts/stan_models/n_reassortants_single_raneff.stan')

numbers_data <- list(N = nrow(data_processed_2),
                    y = data_processed_2 %>% pull(n_reassortants),
                    K = data_processed_2 %>% pull(n_reassortants) %>% max(),
                    C = data_processed_2 %>% pull(collection_regionname) %>% n_distinct(),
                    Y = data_processed_2 %>% pull(collection_year) %>% n_distinct(),
                    continent_index = data_processed_2 %>% pull(collection_regionname) %>% as.factor() %>% as.numeric(),
                    year_index = data_processed_2 %>% pull(collection_year) %>% as.factor() %>% as.numeric(),
                    cases =  data_processed_2 %>% pull(woah_susceptibles_monthly_log1p),
                    sequences =  data_processed_2 %>% pull(n_sequences_log1p))


# Run the model
CHAINS <- 4
CORES <- 4
ITER <- 4000
BURNIN <- ITER/10 # Discard 10% burn in from each chain
SEED <- 4472

numbers_model_2 <- numbers_mod$sample(
  data = numbers_data,
  seed = SEED,
  chains = CHAINS,
  parallel_chains = CORES,
  iter_warmup = BURNIN,
  iter_sampling = ITER)
    

numbers_parms <- numbers_model_2$summary()

############################################## END #################################################
####################################################################################################
####################################################################################################
# Deprecated

# Marginal probability density of the number of unique reassortants/region
# ie, irrespective of class
count_prob <- jointmodel_temp %>% 
  epred_draws(newdata = count_data %>%
                select(collection_regionname) %>%
                distinct(),resp = "nreassortants",
              re_formula = NA)

count_prob %>% 
  median_hdi(.epred)

ggplot(count_prob, aes(x = .epred, fill = collection_regionname , y = collection_regionname)) +
  stat_halfeye() +
  theme_minimal() + 
  scale_x_continuous('Number of Unique Reassortants/Year') +
  
  # Scales
  scale_y_discrete('Region of Origin', labels = function(x) str_wrap(x, width = 20) %>% str_to_title())+
  scale_colour_manual(
    'Reassortant Class',
    values = region_colour %>% pull(Trait, name = Name))  +
  theme(legend.position = 'none') + 
  coord_cartesian(xlim = c(0,10))

# Marginal probability of reassortant class/region
# ie, irrespective of number
class_prob <- jointmodel_temp %>% 
  epred_draws(newdata = count_data %>%
                select(collection_regionname) %>%
                distinct(),
              resp = "group2",
              re_formula = NA)
class_prob %>% 
  median_hdi(.epred)

ggplot(class_prob, aes(x = collection_regionname, colour = .category, y = .epred)) +
  geom_boxplot() +
  theme_minimal() + 
  scale_y_continuous('Posterior Probability of Reassortant Class', labels = scales::percent) +
  
  # Scales
  scale_x_discrete('Region of Origin', labels = function(x) str_wrap(x, width = 20) %>% str_to_title())+
  scale_colour_manual(
    'Reassortant Class',
    values = riskgroup_colour %>% pull(Trait, name = group2)) +
  theme(legend.position = 'bottom')


#### Contrasts ####
fit_year %>% 
  emmeans(., ~ collection_regionname ,
          epred = TRUE,
          resp = 'group2') %>% 
  filter(rep.meas == 'nreassortants') %>%
  contrast(method = 'revpairwise')



joint_prob <- count_prob %>%
  left_join(class_prob, 
            by = c(".draw",  '.row', "collection_regionname"), 
            suffix = c("_count", "_class")) %>%
  mutate(joint_epred = .epred_count * .epred_class)