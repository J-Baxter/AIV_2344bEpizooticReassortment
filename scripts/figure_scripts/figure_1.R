####################################################################################################
####################################################################################################
## Script name: Paper Figure 1
##
## Purpose of script: To plot paper figures 1
##
## Date created: 2025-01-14
##
##
########################################## SYSTEM OPTIONS ##########################################
options(scipen = 6, digits = 7) 
memory.limit(30000000) 

  
########################################## DEPENDENCIES ############################################
# Packages
library(tidyverse)
library(magrittr)
library(RColorBrewer)
library(ggtree)
library(ggtreeExtra)
library(scales)
library(cowplot)

# User functions


############################################## DATA ################################################
combined_data <- read_csv('./2024Aug18/treedata_extractions/2024-09-20_combined_data.csv')
summary_data <- read_csv('./2024Aug18/treedata_extractions/summary_reassortant_metadata_20240904.csv') %>%
  select(-c(cluster_label,
            clade)) 

meta <- read_csv('./2024-09-09_meta.csv') 


############################################## MAIN ################################################

# Number of Sequences ~ Continent
meta %>%
  mutate(collection_regionname = case_when(grepl('europe', collection_regionname) ~ 'europe',
                                           grepl('africa', collection_regionname) ~ 'africa',
                                           grepl('asia', collection_regionname) ~ 'asia',
                                           grepl('(central|northern) america', collection_regionname) ~ 'central & northern america',
                                           grepl('south america|southern ocean', collection_regionname) ~ 'south america',
                                           grepl('australia', collection_regionname) ~ 'australasia',
                                           .default = collection_regionname
  )) %>%
  drop_na(collection_regionname) %>%
  summarise(n = n(), .by = c(collection_regionname)) %>%
  ggplot() + 
  geom_bar(aes(x = collection_regionname, y = n), stat = 'identity') 


# Number of Sequences ~ Continent (cumulatively)
all <- meta %>%
  mutate(collection_regionname = case_when(grepl('europe', collection_regionname) ~ 'europe',
                                           grepl('africa', collection_regionname) ~ 'africa',
                                           grepl('asia', collection_regionname) ~ 'asia',
                                           grepl('(central|northern) america', collection_regionname) ~ 'central & northern america',
                                           grepl('south america|southern ocean', collection_regionname) ~ 'south america',
                                           grepl('australia', collection_regionname) ~ 'australasia',
                                           .default = collection_regionname)) %>%
  drop_na(collection_regionname) %>%
  select(collection_regionname, collection_dateyear) %>%
  expand(collection_regionname, collection_dateyear = full_seq(collection_dateyear,1))

  
plt_1b <- meta %>%
  mutate(collection_regionname = case_when(grepl('europe', collection_regionname) ~ 'europe',
                                           grepl('africa', collection_regionname) ~ 'africa',
                                           grepl('asia', collection_regionname) ~ 'asia',
                                           grepl('(central|northern) america', collection_regionname) ~ 'central & northern america',
                                           grepl('south america|southern ocean', collection_regionname) ~ 'south america',
                                           grepl('australia', collection_regionname) ~ 'australasia',
                                           .default = collection_regionname)) %>%
  drop_na(collection_regionname) %>%
  summarise(n = n(), .by = c(collection_regionname, collection_dateyear)) %>%
  right_join(all) %>%
  mutate(n = case_when(is.na(n) ~ 0, .default = n)) %>%
  arrange(collection_dateyear) %>%
  group_by(collection_regionname) %>%
  mutate(cum_sum = cumsum(n)) %>%
  ggplot() + 
  geom_area(aes(x = collection_dateyear, y = cum_sum, fill = collection_regionname)) +
  theme(legend.position = 'inside',
        legend.position.inside = c(0.2, 0.8))



# Number of Sequences ~ Host Type
meta %>%
  mutate(collection_regionname = case_when(grepl('europe', collection_regionname) ~ 'europe',
                                           grepl('africa', collection_regionname) ~ 'africa',
                                           grepl('asia', collection_regionname) ~ 'asia',
                                           grepl('(central|northern) america', collection_regionname) ~ 'central & northern america',
                                           grepl('south america|southern ocean', collection_regionname) ~ 'south america',
                                           grepl('australia', collection_regionname) ~ 'australasia',
                                           .default = collection_regionname
  )) %>%
  drop_na(host_simplifiedhost) %>%
  summarise(n = n(), .by = c(host_simplifiedhost)) %>%
  ggplot() + 
  geom_bar(aes(x = host_simplifiedhost, y = n), stat = 'identity') 


# Number of Sequences ~ Host Type (cumulatively)
all <- meta %>%
  drop_na(host_simplifiedhost) %>%
  select(host_simplifiedhost, collection_dateyear) %>%
  expand(host_simplifiedhost, collection_dateyear = full_seq(collection_dateyear,1))


plt_1c <- meta %>%
  drop_na(host_simplifiedhost) %>%
  summarise(n = n(), .by = c(host_simplifiedhost, collection_dateyear)) %>%
  right_join(all) %>%
  mutate(n = case_when(is.na(n) ~ 0, .default = n)) %>%
  arrange(collection_dateyear) %>%
  group_by(host_simplifiedhost) %>%
  mutate(cum_sum = cumsum(n)) %>%
  ggplot() + 
  geom_area(aes(x = collection_dateyear, y = cum_sum, fill = host_simplifiedhost)) +
  theme(legend.position = 'inside',
        legend.position.inside = c(0.2, 0.8))


# Reassortants ~ Time
plt_1d <- combined_data %>% 
  filter(segment == 'ha') %>%
  select(TMRCA) %>%
  drop_na(TMRCA) %>%
  mutate(tmrca_date = date_decimal(TMRCA)) %>%
  mutate(tmrca_yearmonth = round_date(tmrca_date, 'month')) %>%
  ggplot() + 
  geom_bar(aes(x = tmrca_yearmonth)) 


# Reassortants ~ Continent (TMRCA)
plt_1e <- summary_data %>% 
  mutate(Continent_of_Earliest_Date = case_when(grepl('europe', Continent_of_Earliest_Date) ~ 'europe',
                                           grepl('africa', Continent_of_Earliest_Date) ~ 'africa',
                                           grepl('asia', Continent_of_Earliest_Date) ~ 'asia',
                                           grepl('(central|northern) america', Continent_of_Earliest_Date) ~ 'central & northern america',
                                           grepl('south america|southern ocean', Continent_of_Earliest_Date) ~ 'south america',
                                           grepl('australia', Continent_of_Earliest_Date) ~ 'australasia',
                                           .default = Continent_of_Earliest_Date
  )) %>%
  drop_na(Continent_of_Earliest_Date) %>%
  summarise(n = n(), .by = c(Continent_of_Earliest_Date)) %>%
  ggplot() + 
  geom_bar(aes(x = Continent_of_Earliest_Date, y = n), stat = 'identity') 
  


# Reassortants ~ Host Type (TMRCA)
plt_1f <- summary_data %>% 
  mutate(Host_of_Earliest_Sample_Date = case_when(grepl('europe', Host_of_Earliest_Sample_Date) ~ 'europe',
                                                grepl('africa', Host_of_Earliest_Sample_Date) ~ 'africa',
                                                grepl('asia', Host_of_Earliest_Sample_Date) ~ 'asia',
                                                grepl('(central|northern) america', Host_of_Earliest_Sample_Date) ~ 'central & northern america',
                                                grepl('south america|southern ocean', Host_of_Earliest_Sample_Date) ~ 'south america',
                                                grepl('australia', Host_of_Earliest_Sample_Date) ~ 'australasia',
                                                .default = Host_of_Earliest_Sample_Date
  )) %>%
  drop_na(Host_of_Earliest_Sample_Date) %>%
  summarise(n = n(), .by = c(Host_of_Earliest_Sample_Date)) %>%
  ggplot() + 
  geom_bar(aes(x = Host_of_Earliest_Sample_Date, y = n), stat = 'identity') 

# Combine to make panel
plt_1rightplots <- align_plots(plt_1d, plt_1e, plt_1b, align = 'v', axis = 'l')

plt_1rightupper <- plot_grid(plt_1rightplots[[3]],  plt_1c, 
                             labels = c('B', 'C'), 
                             label_size = 12)

plt_1rightlower <- plot_grid(plt_1rightplots[[2]], plt_1f,
                             labels = c('E', 'F'),
                             label_size = 12)

plt_1right <- plot_grid(plt_1rightupper, plt_1rightplots[[1]], plt_1rightlower,
                        labels = c('', 'D', ''), 
                        label_size = 12,
                        ncol = 1)

plt_1right
############################################## WRITE ###############################################




############################################## END #################################################
####################################################################################################
####################################################################################################