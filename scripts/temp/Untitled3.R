####################################################################################################
####################################################################################################
## Script name: 
##
## Purpose of script: 
##
## Date created: 2025-05-23
##
##
########################################## SYSTEM OPTIONS ##########################################
options(scipen = 6, digits = 7) 
memory.limit(30000000) 


########################################## DEPENDENCIES ############################################
# Packages
library(tidyverse)
library(magrittr)

############################################## DATA ################################################
data0 <- read_csv("~/Downloads/OneDrive_1_21-05-2025/RegionalTre_mrca_stats_disp_jumps_combined.csv")


############################################## MAIN ################################################

data1 <- data0 %>% 
  select(data_name, count_cross_species, count_to_mammal) %>% 
  separate(data_name, into = c("segment", "continent"), sep = "_", remove = FALSE)

unique(data1$continent)


count.cross.species <- data1 %>% 
  filter(!is.na(count_cross_species))

count.to.mammal <- data1 %>% 
  filter(!is.na(count_to_mammal))


## cross species frequency
stat.cross.species <- count.cross.species %>% 
  #filter(segment == 'ha') %>%
  group_by(continent, segment) %>% 
  summarise(stat.cross.species = sum(count_cross_species)) %>% 
  ungroup()


## spread into mammal species
stat.to.mammal <- count.to.mammal %>% 
  group_by(continent, segment) %>% 
  summarise(stat.to.mammal = sum(count_to_mammal)) %>% 
  ungroup()



join0 <- stat.cross.species %>% 
  left_join(stat.to.mammal, by = c("continent", "segment")) 


total_mammals <- meta %>%
  
  mutate(collection_regionname = case_when(grepl('europe', collection_regionname) ~ 'europe',
                                           grepl('africa', collection_regionname) ~ 'africa',
                                           grepl('asia', collection_regionname) ~ 'asia',
                                           grepl('(central|northern) america', collection_regionname) ~ 'central & northern america',
                                           .default = collection_regionname)) %>%
  drop_na(collection_regionname) %>%
  group_by(collection_regionname) %>%
  mutate(total_seqs = n()) %>% 
  select(collection_regionname, total_seqs, host_class) %>%
  mutate(host_class = factor(host_class, levels = c("mammalia" ,   "aves"  ,      "environment" ,"unknown"))) %>%
  count(host_class, total_seqs, .drop = F) %>%
  fill(total_seqs, .direction = 'updown') %>%
  filter(host_class == 'mammalia') %>%
  rename(continent = collection_regionname, 
         total_mammal = n) 

plt_a <- join0 %>%
  mutate(continent = case_when(grepl('southamerica', continent) ~ 'south america',
                               grepl('northamerica', continent) ~ 'central & northern america',
                               .default = continent)) %>%
  #slice_max(tibble(stat.cross.species, stat.to.mammal), by = continent) %>%
  summarise(stat.cross.species = list(quantile(stat.cross.species, na.rm = T, probs = c(0.25, 0.5, 0.75))),
            stat.to.mammal = list(quantile(stat.to.mammal, na.rm = T, probs = c(0.25, 0.5, 0.75))),
            .by = continent) %>%
  unnest_longer(c('stat.cross.species', 'stat.to.mammal')) %>%
  select(-stat.to.mammal_id) %>%
  pivot_longer(c('stat.cross.species', 'stat.to.mammal'), names_to = 'name', values_to = 'jumps') %>%
  pivot_wider(names_from = stat.cross.species_id, values_from = jumps) %>%
  
  left_join(total_mammals) %>%
  rowwise() %>%
  #mutate(scaled_species = stat.cross.species/total_seqs,
  # scaled_mammal = stat.to.mammal/total_seqs) %>%
  mutate(across(ends_with('%'), .fns = ~.x/total_seqs*100)) %>%
  # pivot_longer(starts_with('scaled'), names_to = 'group', values_to = 'frequency') %>%
  as_tibble() %>%
  mutate(across(3:5, ~replace_na(.x, 0))) %>%
  
  ggplot(aes(x = continent, fill = name)) + 
  
  geom_bar(aes(y = `50%`),
           stat = "identity",
           color = "grey20", 
           position = position_dodge(0.7), 
           width = 0.65) + 
  geom_errorbar(aes(ymin = `25%`, ymax =`75%`),  position = position_dodge(0.7), 
                width = 0.4) +
  scale_fill_manual('', values=c( 'stat.to.mammal' ='#FDB462',
                                  'stat.cross.species' = "#80B1D3"), 
                    labels = c('stat.to.mammal' = 'Transitions to mammal',
                               'stat.cross.species' = 'All transitions')) + 
  scale_y_continuous('Transitions per 100 Sequences', expand = c(0.001,0),limits = c(0,6), breaks = seq(0,6,by = 2),  labels = seq(0,6,by = 2)) +
  scale_x_discrete('Continent',labels = ~str_to_title(.x) %>% str_wrap(width = 20)) +
  global_theme + 
  theme(legend.position  = 'inside',
        legend.position.inside = c(0.25,0.9),
        strip.text = element_text(face = 'bold', size = 10),
        strip.background = element_blank(),
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 9),
        legend.text = element_text(size = 8))


plt_b <- join0 %>%
  mutate(continent = case_when(grepl('southamerica', continent) ~ 'south america',
                               grepl('northamerica', continent) ~ 'central & northern america',
                               .default = continent)) %>%
  #slice_max(tibble(stat.cross.species, stat.to.mammal), by = continent) %>%
  summarise(#stat.cross.species = list(quantile(stat.cross.species, na.rm = T, probs = c(0.25, 0.5, 0.75))),
    stat.to.mammal = list(quantile(stat.to.mammal, na.rm = T, probs = c(0.25, 0.5, 0.75))),
    .by = continent) %>%
  unnest_longer('stat.to.mammal') %>%
  #select(-stat.to.mammal_id)
  left_join(total_mammals) %>%
  rowwise() %>%
  mutate(scaled_mammal = stat.to.mammal/total_mammal*100) %>%
  select(-c(total_seqs,stat.to.mammal, total_mammal)) %>%
  pivot_wider( names_from = 'stat.to.mammal_id', values_from = 'scaled_mammal') %>%
  mutate(across(3:5, ~replace_na(.x, 0))) %>%
  mutate(name = 'stat.to.mammal') %>%
  #filter(group == 'scaled_mammal') %>%
  ggplot(aes(x = continent, fill = name)) + 
  
  geom_bar(aes(y = `50%`),
           stat = "identity",
           color = "grey20", 
           position = position_dodge(0.7), 
           width = 0.65) + 
  geom_errorbar(aes(ymin = `25%`, ymax =`75%`),  position = position_dodge(0.7), 
                width = 0.4) +
  scale_fill_manual('', values=c( 'stat.to.mammal' ='#FDB462',
                                  'stat.cross.species' = "#80B1D3"), 
                    labels = c('stat.to.mammal' = 'Transitions to mammal',
                               'stat.cross.species' = 'All transitions')) + 
  scale_y_continuous('Transitions per 100 Sequences', expand = c(0.001,0),limits = c(0,8), breaks = seq(0,8,by = 2),  labels = seq(0,8,by = 2)) +
  scale_x_discrete('Continent',labels = ~str_to_title(.x) %>% str_wrap(width = 20)) +
  global_theme + 
  theme(legend.position  = 'inside',
        legend.position.inside = c(0.25,0.9),
        strip.text = element_text(face = 'bold', size = 10),
        strip.background = element_blank(),
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 9),
        legend.text = element_text(size = 8))

plot_grid(plt_a, plt_b, nrow = 1, labels = 'AUTO', label_size = 9, align = 'hv')

############################################## WRITE ###############################################
ggsave('~/Downloads/flu_plots/host-jumps.jpeg', height = 12, width = 22, units = 'cm', dpi = 360)

############################################## END #################################################
####################################################################################################
####################################################################################################