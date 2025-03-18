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
library(geomtextpath)

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

ggplot(cdf01, aes(Sepal.Width)) + 
  stat_ecdf(geom = "step", color="purple")



t <- combined_data %>% 
  filter(segment == 'ha') %>%
  select(cluster_profile, collection_regionname, TMRCA) %>%
  mutate(tmrca_date = date_decimal(TMRCA) %>% as_date()) %>%
  filter(tmrca_date >= as_date('2019-01-01')) %>%
  arrange(tmrca_date) %>%
  mutate(collection_regionname = case_when(grepl('europe', collection_regionname) ~ 'europe',
                                           grepl('africa', collection_regionname) ~ 'africa',
                                           grepl('asia', collection_regionname) ~ 'asia',
                                           grepl('(central|northern) america|caribbean', collection_regionname) ~ 'central & northern america',
                                           grepl('south america|southern ocean', collection_regionname) ~ 'south america',
                                           grepl('australia|melanesia', collection_regionname) ~ 'australasia',
                                           .default = collection_regionname)) %>%
  #dplyr::select(-cluster_profile) %>%
  summarise(n = n(), .by = c(collection_regionname, tmrca_date)) %>%
  group_by(collection_regionname) %>%
  mutate(cs = cumsum(n)) %>%
  ungroup()  %>%
  select(-n)


plt_2a <- expand_grid(tmrca_date = seq(ymd('2019-01-01'),ymd('2024-05-01'),by='day'),
                      collection_regionname = c('europe', 'africa', 'asia', 'central & northern america'),
                      cs = NA_real_) %>%
  rows_patch(t, by = c('collection_regionname', 'tmrca_date')) %>%
  group_by(collection_regionname) %>%
  fill(cs) %>%
  ungroup() %>%
  mutate(cs = replace_na(cs, 0)) %>%
  
  ggplot() + 
  geom_step(aes(x=tmrca_date, y=cs, color=collection_regionname), lwd = 0.8) +
  scale_colour_manual('Continent', values = region_colours, labels = str_to_title) +
  scale_y_continuous('Reassortans (Cumulative)', expand = c(0.01, 0)) + 
  scale_x_date(limits = as_date(c('2019-01-01', '2024-02-01')),
                   breaks = '1 year', 
                   date_labels = "%Y", 'Date (Years)',
                   expand = c(0,0)) + 
  theme_classic() + 
  theme(legend.position = 'inside',
        legend.position.inside = c(0,1),
        legend.justification.inside = c(0,1),
        legend.background = element_blank(),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 10),
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 10))


# 'Diversity' -  the exponent of Shannon entropy 
plt_2b <- h5_diversity_sliding_window %>%
  filter(! continent %in%  c('South America', 'Antarctica')) %>%
  mutate(continent = str_to_lower(continent) %>%
           case_when(grepl('north america', .) ~ 'central & northern america',
                     .default = .)) %>%
  mutate(midpoint = (start_point + end_point)/2) %>%
  mutate(midpoint_date = date_decimal(midpoint)) %>%
  ggplot(aes(y = diversity, x = midpoint_date, colour = continent, fill = continent,label = str_to_title(continent))) +
  #geom_point() + 
  geom_smooth(method = 'gam' , formula = y ~ s(x, bs = 'cs'), se = TRUE, alpha = 0.1)  +
  #geom_labelsmooth(text_smoothing = 30, 
                  # fill = "white",
                   #formula = y ~ s(x, bs = 'cs'), 
                   #method = "gam",
                   #size = 4, 
                   #linewidth = 0.1,
                   #boxlinewidth = 0.3) +
  scale_fill_manual('Continent', values = region_colours, labels = str_to_title) + 
  scale_colour_manual('Continent', values = region_colours, labels = str_to_title) +
  scale_y_continuous('Numbers-Equivalent Shanon Entropy', expand = c(0.01, 0)) + 
  scale_x_datetime(limits = as_datetime(c('2020-06-01', '2024-02-01')),
               breaks = '1 year', 
               date_labels = "%Y", 'Date (Years)',
               expand = c(0,0)) + 
  
  theme_classic()+ 
  theme(legend.position = 'none',
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 10))


# Genetic distance between reassortants

# Combine

plt_2 <- cowplot::plot_grid(plt_2a, plt_2b, nrow = 1, labels = 'AUTO', align = 'v', axis = 'tb')
plt_2
############################################## WRITE ###############################################
ggsave('~/Downloads/flu_plots/figure2.jpeg', height = 12.5, width = 25, units = 'cm', dpi = 360)




############################################## END #################################################
####################################################################################################
####################################################################################################