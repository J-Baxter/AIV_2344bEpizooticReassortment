####################################################################################################
####################################################################################################
## Script name:
##
## Purpose of script: Figure 2 - plot numnber of reassortments over time alongside Will's diversity
## metrics
##
## Date created: 2025-02-28
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


############################################## DATA ################################################

h5_diversity_sliding_window <- read_csv('./h5_diversity_sliding_window.csv')
combined_data <- read_csv('./2024Aug18/treedata_extractions/2024-09-20_combined_data.csv')

############################################## MAIN ################################################
region_colours <- c('europe' = '#1b9e77',
                    'asia' ='#d95f02',
                    'africa' ='#7570b3',
                    'australasia' = '#e7298a',
                    'central & northern america' ='#66a61e',
                    'south america' ='#e6ab02')

# Cumulative number of reassortants
combined_data %>% 
  select(cluster_profile, collection_regionname, TMRCA) %>%
  mutate(tmrca_date = date_decimal(TMRCA) %>% as_date()) %>%
  filter(tmrca_date >= as_date('2019-01-01')) %>%
  arrange(tmrca_date) %>%
  mutate(n = 1) %>%
  mutate(collection_regionname = case_when(grepl('europe', collection_regionname) ~ 'europe',
                                           grepl('africa', collection_regionname) ~ 'africa',
                                           grepl('asia', collection_regionname) ~ 'asia',
                                           grepl('(central|northern) america|caribbean', collection_regionname) ~ 'central & northern america',
                                           grepl('south america|southern ocean', collection_regionname) ~ 'south america',
                                           grepl('australia|melanesia', collection_regionname) ~ 'australasia',
                                           .default = collection_regionname)) %>%
  group_by(collection_regionname) %>%
  mutate(cum_sum = cumsum(n))  %>%
  ggplot() + 
  geom_line(aes(x = tmrca_date, y = cum_sum, colour = collection_regionname)) + 
  scale_colour_manual('Continent', values = region_colours, labels = str_to_title) +
  scale_y_continuous('Reassortans (Cumulative)', expand = c(0.01, 0)) + 
  scale_x_date('Date', expand = c(0.01, 0)) +
  theme_classic() + 
  theme(legend.position = 'none')


combined_data# 'Diversity' -  the exponent of Shannon entropy 
h5_diversity_sliding_window %>%
  filter(! continent %in%  c('South America', 'Antarctica')) %>%
  mutate(continent = str_to_lower(continent) %>%
           case_when(grepl('north america', .) ~ 'central & northern america',
                     .default = .)) %>%
  mutate(midpoint = (start_point + end_point)/2) %>%
  mutate(midpoint_date = date_decimal(midpoint)) %>%
  ggplot(aes(y = diversity, x = midpoint_date, colour = continent, fill = continent)) +
  #geom_point() + 
  geom_smooth(method = 'gam' , formula = y ~ s(x, bs = 'cs'), se = TRUE, alpha = 0.1)  +
  scale_fill_manual('Continent', values = region_colours, labels = str_to_title) + 
  scale_colour_manual('Continent', values = region_colours, labels = str_to_title) +
  scale_y_continuous('Numbers-Equivalent Shanon Entropy', expand = c(0.01, 0)) + 
  scale_x_datetime('Date', expand = c(0.01, 0)) +
  theme_classic() + 
  theme(legend.position = 'none')





############################################## WRITE ###############################################




############################################## END #################################################
####################################################################################################
####################################################################################################