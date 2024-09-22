####################################################################################################
####################################################################################################
# The aim of this model is to determine associations between variables obtained from our phylogenetic
# analysis and difusion coefficient for each reassortant

# variables of interest include: 
# 1. evolutionary rates
# 2. persistence time
# 3. region of origin
# 4. persistence time in wild birds
# 5. persistence time in domestic birds
# 6. total number of species jumps

# This script implements a workflow for a zero-inflated lognormal model. Initial exploratory analysis 
# is conducted using OLS regression, thereafter progressing to Bayesian analysis using BRMS.

########################################## DEPENDENCIES ############################################
library(broom)
library(broom.mixed)
library(tidybayes)
library(bayesplot)


########################################### IMPORT DATA ############################################

combined_data <- read_csv('./2024Aug18/treedata_extractions/2024-09-20_combined_data.csv')


########################################### FORMAT DATA ############################################

model_data <- combined_data %>%
  
  # select variables of interes
  select(
    cluster_profile,
    group2,
    weighted_diff_coeff,
    original_diff_coeff,
    evoRate,
    persist.time,
    collection_regionname,
    host_simplifiedhost,
    count_cross_species,
    starts_with('median'),
    starts_with('max')) %>%
  
  # Substitute NA values in diffusion coefficient with 0
  replace_na(list(original_diff_coeff = 0,
                  weighted_diff_coeff = 0,
                  `median_galliformes-domestic` = 0,
                  `median_other-bird` = 0,
                  `median_anseriformes-wild` = 0,
                  `median_galliformes-wild` = 0,
                  median_environment = 0,
                  median_mammal = 0,
                  `median_charadriiformes-wild` = 0,
                  median_human = 0,
                  `median_anseriformes-domestic` = 0)) %>%
  rename_with(~gsub('-', '_', .x))



################################### Ancestral host ~ Region #####################################

combined_data %>%
  mutate(collection_regionname = case_when(grepl('europe', collection_regionname) ~ 'europe',
                                           grepl('africa', collection_regionname) ~ 'africa',
                                           grepl('asia', collection_regionname) ~ 'asia',
                                           grepl('(central|northern) america', collection_regionname) ~ 'central & northern america',
                                           .default = collection_regionname )) %>%
  filter(!is.na(collection_regionname)) %>% 
  filter(!grepl('\\+', host_simplifiedhost)) %>%
  ggplot(aes(x = collection_regionname,
             fill = host_simplifiedhost)) +
  
  # Geom objects
  geom_bar(position = "fill") +
  scale_y_continuous('Number of Reassortants', labels = scales::percent) +
  
  # Scales
  scale_x_discrete('Region of Origin', labels = function(x) str_wrap(x, width = 20) %>% str_to_title())+
  scale_fill_brewer('Reassortant Class', direction = -1) +
  
  
  
  # Graphical
  facet_grid(rows = vars(segment)) + 
  theme_minimal() 
  


  
props <- combined_data %>%
  select(c(collection_regionname,
           segment,
           host_simplifiedhost,
           group2)) %>%
  filter(!grepl('\\+', host_simplifiedhost)) %>%
  mutate(collection_regionname = case_when(grepl('europe', collection_regionname) ~ 'europe',
                                           grepl('africa', collection_regionname) ~ 'africa',
                                           grepl('asia', collection_regionname) ~ 'asia',
                                           grepl('(central|northern) america', collection_regionname) ~ 'central & northern america',
                                           .default = collection_regionname
  )) %>%
  group_by(collection_regionname, segment) %>%
  group_by(host_simplifiedhost, .add = TRUE) %>%
  summarise(prop = n())  %>%
  filter(segment == 'ha') %>%
  ungroup()%>%
  select(-segment) %>%
  mutate(prop = as.numeric(prop)) %>%
  pivot_wider(values_from = prop, names_from = host_simplifiedhost) %>%  # Mood to columns
  replace_na(., list(major = 0)) %>%
  pivot_longer(., -collection_regionname, names_to = 'host_simplifiedhost', values_to = 'prop')


MASS::loglm(prop ~ host_simplifiedhost + collection_regionname, props)
#######################################  Reassortant Class ~ Region  ########################################
combined_data %>%
  mutate(collection_regionname = case_when(grepl('europe', collection_regionname) ~ 'europe',
                                           grepl('africa', collection_regionname) ~ 'africa',
                                           grepl('asia', collection_regionname) ~ 'asia',
                                           grepl('(central|northern) america', collection_regionname) ~ 'central & northern america',
                                           .default = collection_regionname
  )) %>%
  filter(!is.na(collection_regionname)) %>%
  ggplot(aes(x = collection_regionname,
             fill = group2)) +
  
  # Geom objects
  geom_bar(position = "fill") +
  scale_y_continuous('Number of Reassortants', labels = scales::percent) +
  
  # Scales
  scale_x_discrete('Region of Origin', labels = function(x) str_wrap(x, width = 20) %>% str_to_title())+
  scale_fill_brewer('Reassortant Class', direction = -1) +
  
  
  
  # Graphical
  facet_grid(rows = vars(segment)) + 
  theme_minimal() 




model_data <- combined_data %>%
    select(c(collection_regionname,
             segment,
             host_simplifiedhost,
             group2)) %>%
  filter(!grepl('\\+', host_simplifiedhost)) %>%
  mutate(collection_regionname = case_when(grepl('europe', collection_regionname) ~ 'europe',
                                           grepl('africa', collection_regionname) ~ 'africa',
                                           grepl('asia', collection_regionname) ~ 'asia',
                                           grepl('(central|northern) america', collection_regionname) ~ 'central & northern america',
                                           .default = collection_regionname
  )) %>%
  mutate(dominant_binary = if_else(group2 == 'dominant', 1, 0),
         segment_binary = if_else(segment %in%  c('ha', 'nx'), 1, 0))

props <- combined_data %>%
  select(c(collection_regionname,
           segment,
           host_simplifiedhost,
           group2)) %>%
  filter(!grepl('\\+', host_simplifiedhost)) %>%
  mutate(collection_regionname = case_when(grepl('europe', collection_regionname) ~ 'europe',
                                           grepl('africa', collection_regionname) ~ 'africa',
                                           grepl('asia', collection_regionname) ~ 'asia',
                                           grepl('(central|northern) america', collection_regionname) ~ 'central & northern america',
                                           .default = collection_regionname
  )) %>%
  group_by(collection_regionname, segment) %>%
  group_by(group2, .add = TRUE) %>%
  summarise(prop = n())  %>%
  filter(segment == 'ha') %>%
  ungroup()%>%
  select(-segment) %>%
  mutate(prop = as.numeric(prop)) %>%
  pivot_wider(values_from = prop, names_from = group2) %>%  # Mood to columns
  replace_na(., list(major = 0)) %>%
  pivot_longer(., -collection_regionname, names_to = 'group2', values_to = 'prop')


MASS::loglm(prop ~ group2 + collection_regionname, props)




####################################################################################################
####################################################################################################