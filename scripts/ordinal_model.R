####################################################################################################
####################################################################################################
## Script name: Class Model (using cumulative distribution)
##
## Purpose of script:to model the probability that any given reassortant belongs to class X,
## stratified by continent
##
## Date created: 2025-05-13
##
##
########################################## SYSTEM OPTIONS ##########################################
options(scipen = 6, digits = 7) 
memory.limit(30000000) 

  
########################################## DEPENDENCIES ############################################
# Packages
library(tidyverse)
library(magrittr)
library(broom)
library(broom.mixed)
library(brms)
library(cmdstanr)


# User functions


############################################## DATA ################################################

reassortant_ancestral_changes <- read_csv('./reassortant_ancestral_changes.csv')

# Obtain 'nearest major reassortant'
library(igraph)
my_edgelist<- reassortant_ancestral_changes %>% 
  select(ends_with('label'), cluster_class) %>%
  drop_na() %>%
  relocate(parent_label, cluster_label) %>% 
  graph_from_data_frame(.)
plot(my_edgelist, vertex.size = 7)
most_recent_major <- distances(my_edgelist,
                               weights = reassortant_ancestral_changes %>%  select(ends_with('label'), segments_changed) %>% drop_na() %>% pull(segments_changed),
                               mode = 'in',
                               to =  c( "H5N8/2019/R7_Africa",
                                        "H5N1/2020/R1_Europe",
                                        "H5N1/2021/R1_Europe",
                                        "H5N1/2021/R3_Europe" , 
                                        "H5N1/2022/R7_NAmerica" ,
                                        "H5N1/2022/R12_Europe")) %>%
  as_tibble(rownames = 'cluster_label') %>%
  rowwise() %>%
  filter(cluster_label != "H5N8/2019/R7_Africa") %>%
  mutate(across(starts_with('H5'), .fns = ~ case_when(.x==0 ~ NaN, 
                                                      abs(.x) == Inf ~ NaN,
                                                      .default = .x))) %>%
  mutate(most_recent_major = list(names(.)[2:7][which(c_across(starts_with('H5')) == min(c_across(starts_with('H5')), na.rm = T))])) %>%
  mutate(most_recent_major = paste(most_recent_major, collapse = '')) %>%
  mutate(last_major_label = na_if(most_recent_major, '')) %>%
  as_tibble() %>%
  select(cluster_label, last_major_label) %>% 
  left_join(meta %>% select(last_major_profile = cluster_profile, last_major_label = cluster_label) %>% distinct())



####################################################################################################
# Update reassortant 'ancestry'
updated <- reassortant_ancestral_changes %>%
  
  # Add  origin continent of current reassortant
  left_join(combined_data %>% filter(segment == 'ha') %>% 
              select(cluster_profile,
                     cluster_region = collection_regionname) %>% 
              group_by(cluster_profile) %>% 
              #slice_min(cluster_tmrca, n = 1) %>% 
              distinct()) %>%
  
  # Add tmrca and origin continent of 'parent' reassortant (with respect to HA)
  left_join(combined_data %>% 
              filter(segment == 'ha') %>% 
              select(cluster_profile, 
                     parent_region = collection_regionname) %>% 
              group_by(cluster_profile) %>% 
              distinct(),
            by = join_by(parent_profile == cluster_profile), 
            relationship = 'many-to-many') %>%
  
  rename(cluster_tmrca = height_median,
         parent_tmrca = parent_height_median) %>%
  
  # Add 'last' (with respect to HA) major reassortant
  left_join(most_recent_major) %>%
  left_join(., select(., last_major_profile = cluster_profile, 
                      last_major_region = cluster_region, 
                      last_major_tmrca = cluster_tmrca)) %>%
  
  mutate(across(ends_with('tmrca'), .fns = ~subtract(2024.21, .x))) %>%
  
  # calculate time to most recent reassortant
  rowwise() %>%
  mutate(time_since_last_major =  cluster_tmrca-last_major_tmrca) %>%
  as_tibble() %>% 
  
  # group continents
  mutate(across(ends_with('region'),
                .fns = ~case_when(grepl('europe', .x) ~ 'europe',
                                  grepl('africa', .x) ~ 'africa',
                                  grepl('asia', .x) ~ 'asia',
                                  grepl('(central|northern) america', .x) ~ 'central & northern america',
                                  .default = .x))) %>%
  
  # Does reassortant coincide with region change
  mutate(region_changed_from_previous = if_else(parent_region != cluster_region, '1', '0')) 


class_data <- updated %>%
  select(cluster_class,
         parent_class,
         cluster_region, 
         segments_changed,
         region_changed_from_previous,
         time_since_last_major,
         time_since_parent) %>%
  
  mutate(cluster_class = ordered(cluster_class, levels = c('minor', 'moderate', 'major')))

############################################## MAIN ################################################
# Formula
ordinal_formula_priors <- get_prior(cluster_class~1 + parent_class + cluster_region + s(time_since_last_major) + segments_changed,
                                  data = class_data,
                                  family = cumulative("probit")) 



# Priors
ordinal_formula_priors$prior[c(1:8)] <- "normal(0,5)"


CHAINS <- 4
CORES <- 4
ITER <- 2000
BURNIN <- ITER/10 # Discard 10% burn in from each chain
SEED <- 4472

ordinal_model <- brm(cluster_class~1 + parent_class + cluster_region + s(time_since_last_major) + segments_changed,
            data=class_data,
            prior = ordinal_formula_priors,
            family =cumulative("probit"),
            chains = CHAINS,
            cores = CORES, 
            threads = 2, 
            backend = "cmdstanr",
            iter = ITER,
            warmup = BURNIN,
            seed = SEED)

#marginal_effects(temp, "parent_class", categorical = TRUE) #Graph
#marginal_effects(temp, "time_since_last_major", categorical = TRUE)
#marginal_effects(temp, "cluster_region", categorical = TRUE)
#marginal_effects(temp, "region_changed_from_previous", categorical = TRUE) # no difference

############################################## WRITE ###############################################


simulatedResponse = ordinal_model %>% 
  predict(type="response") %>% 
  apply(1, \(x) sample(1:length(x), prob = x, size=1000, replace = TRUE)) %>% 
  t()
model.check <-createDHARMa(
  simulatedResponse = simulatedResponse,
  observedResponse = as.integer(class_data %>% drop_na() %>% pull(cluster_class)),
  fittedPredictedResponse = ordinal_model %>% predict(), # the linear predictor
  integerResponse = !all(as.numeric(simulatedResponse)%%1!=0) # checks if prediction consists of integers
)

plot(model.check)

############################################## END #################################################
####################################################################################################
####################################################################################################