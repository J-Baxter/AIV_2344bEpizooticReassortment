####################################################################################################
####################################################################################################
## Script name: Class Model
##
## Purpose of script: to model the probability that any given reassortant belongs to class X,
## stratified by continent
##
## Date created: 2025-03-24
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

# variables required
# - class of reassortant /
# - immediately preceeding reassortant /
# - time since last ressortant /
# - continent /
# (- diversity)
# - proxy incidence
# - number of genomes


############################################## DATA ################################################
combined_data <- read_csv('./2024Aug18/treedata_extractions/2024-09-20_combined_data.csv')
summary_data <- read_csv('./2024Aug18/treedata_extractions/summary_reassortant_metadata_20240904.csv') %>%
  dplyr::select(-c(cluster_label,
                   clade)) 

final_clusters <- read_csv('./final_clusters_2025Mar21.csv')


############################################## MAIN ################################################
class_data <- combined_data %>% 
  dplyr::select(-c(group2, col2)) %>%
  left_join(final_clusters) %>%
  
  dplyr::select(cluster_profile, 
                segment,
                TMRCA,
                collection_regionname,
                class) %>%
  
  mutate(collection_regionname = case_when(grepl('europe', collection_regionname) ~ 'europe',
                                           grepl('africa', collection_regionname) ~ 'africa',
                                           grepl('asia', collection_regionname) ~ 'asia',
                                           grepl('(central|northern) america', collection_regionname) ~ 'central & northern america',
                                           .default = collection_regionname
  )) %>%
  
  mutate(TMRCA = date_decimal(TMRCA) %>%
           round_date(unit = 'day')%>%
           as_date(.)) %>%
  
  filter(segment== 'ha') %>%
  arrange(TMRCA) %>%
  group_by(collection_regionname) %>%
  
  # class of most recent prior reassortant
  mutate(previous_class = lag(class),
         previous_class = replace_na(previous_class, 'none')) %>%
  
  # time since last reassortant
  # calculate as interval, then convert to decimal of years
  mutate(time_since_previous = interval(lag(TMRCA), TMRCA) %>%
           as.numeric('years')) %>%
  
  # revert formatting of dates for model
  # TMRCA to decimal, then subtract 2016 for scale
  mutate(tmrca_decimal = decimal_date(TMRCA),
         tmrca_decimal_adjusted = tmrca_decimal - 2016) %>%
  
  # select vars of interest
  dplyr::select(collection_regionname,
                class, 
                previous_class,
                tmrca_decimal_adjusted,
                time_since_previous) %>%
  
  # exclude the first reassortant from each continent
  drop_na(time_since_previous) %>%
  drop_na(class) %>%
  
  filter(collection_regionname != 'africa')


# Plot Data
ggplot(class_data) + geom_bar(aes(x = class, fill = previous_class), position = 'stack') + facet_wrap( ~ collection_regionname)
ggplot(class_data, aes(x = time_since_previous)) + geom_histogram() + 
  #geom_smooth(method = "lm", se = FALSE) +
  facet_wrap( class ~ collection_regionname)


# Formula
class_formula_priors <- get_prior(class ~  previous_class:collection_regionname + time_since_previous,
                                  data = class_data,
                                  family = categorical(link ='logit')) 



# Priors
class_formula_priors$prior[c(1:14)] <- "normal(0,5)"
class_formula_priors$prior[c(16:29)] <- "normal(0,5)"


# Model
CHAINS <- 4
CORES <- 4
ITER <- 4000
BURNIN <- ITER/10 # Discard 10% burn in from each chain
SEED <- 4472


class_model <- brm(
  class ~  previous_class:collection_regionname + time_since_previous,
  prior = class_formula_priors,
  data = class_data,
  family = categorical(link ='logit'),
  chains = CHAINS,
  cores = CORES, 
  threads = 2, 
  backend = "cmdstanr",
  iter = ITER,
  warmup = BURNIN,
  seed = SEED,
  control = list(adapt_delta = 0.95)
)



############################################## WRITE ###############################################




############################################## END #################################################
####################################################################################################
####################################################################################################
