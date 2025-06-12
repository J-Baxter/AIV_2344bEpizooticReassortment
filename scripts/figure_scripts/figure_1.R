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
library(treeio)
library(TreeTools)
library(ggnewscale)
library(rnaturalearth)
library(rnaturalearthdata)
library(ggstream)
# User functions


############################################## DATA ################################################
combined_data <- read_csv('./2024Aug18/treedata_extractions/2024-09-20_combined_data.csv')
summary_data <- read_csv('./2024Aug18/treedata_extractions/summary_reassortant_metadata_20240904.csv') %>%
  select(-c(cluster_label,
            clade)) 

meta <- read_csv('./2024-09-09_meta.csv') 

new_tree <- read.beast('./2025Feb26/globalsubsample/ha_global_SRD06_relaxLn_constant_mcc.tree')

hpai_cases <- read_csv('~/Downloads/overview-raw-data_202502241440.csv') 
woah_hpai <- read_csv('~/Downloads/Quantitative data 2025-02-25.csv')

ref <- ne_countries() %>%
  as_tibble()  %>% 
  select(name, continent)

############################################## MAIN ################################################
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


# Host Palette
host_colours <- c(
 'anseriformes-domestic' = '#a6cee3',
 'anseriformes-wild' = '#1f78b4',
 'galliformes-domestic' = '#b2df8a',
 'galliformes' = '#33a02c',
 'mammal' = '#fb9a99',
 'human' = '#e31a1c',
 'charadriiformes-wild' = '#fdbf6f',
 'other-bird' = '#ff7f00',
 'unknown' = '#cab2d6',
 'environment' = '#6a3d9a')


#E8E1E9FF #C0A5AAFF #4D3944FF #7083A4FF #B3A2B4FF #C9CCEAFF #3B3960FF  

#host_colours <- c(
 # 'anseriformes-domestic' = '#E8E1E9FF',
 # 'anseriformes-wild' = '#C0A5AAFF',
 # 'galliformes-domestic' = '#4D3944FF',
 # 'galliformes-wild' = '#7083A4FF',
 # 'mammal' = '#B3A2B4FF',
 # 'human' = '#C9CCEAFF',
 # 'charadriiformes-wild' = '#3B3960FF',
 # 'other-bird' = '#1E2142FF',
 # 'unknown' = '#586085FF',
 # 'environment' = '#F6E0D2FF')


# Region Palette
region_colours <- c('europe' = '#1b9e77',
            'asia' ='#d95f02',
            'africa' ='#7570b3',
            'australasia' = '#e7298a',
            'central & northern america' ='#66a61e',
            'south america' ='#e6ab02')


#region_colours <- c('europe' = '#EFCBCBFF',
#                      'asia' ='#B47880FF',
#                      'africa' ='#824B51FF',
#                      'australasia' = '#635761FF',
#                      'central & northern america' ='#AD616CFF',
#                      'south america' ='#D1A391FF')


#F6E0D2FF #DFA398FF #9C6755FF #659794FF #EA967CFF #F5C98EFF #D65B5AFF #586085FF 


# Format WOAH data and estimate the minimum number of cases

  # set groupings
  #group_by(collection_regionname, interval) %>%
 # mutate(scale_factor = woah_cases/sum(fao_cases)) %>%
 # ungroup() %>%
 # rowwise() %>%
 # mutate(scaled_estimate = fao_cases*scale_factor) %>%
 # filter(!grepl('australasia', collection_regionname)) %>%
 # filter(collection_datemonth >= as_date('2019-01-01'))


# Plot
plt_1a <- new_tree %>%
  mutate(isolate_id = str_extract(label, "EPI_ISL_(china_){0,1}\\d+[^.|]*")) %>%
  # mutate(date = str_extract(label, "\\d{4}-.*")) %>%
  left_join(reassortant_profiles) %>%
  left_join(continent) %>%
  #left_join(mrca_nodes) %>%
  
  ggtree(mrsd = "2024-03-18") + 
  theme_tree2(#base_family = "LM Sans 10",
    #plot.margin = unit(c(1,1,1,1), units = "cm"),
    axis.text.x = element_text(size = 8),
    axis.title.x = element_text(size = 10)
  ) +
  
  #theme_classic() + 
  
  scale_x_continuous(
    # limits = c(2014.9, Inf),
    'Time (Years)',
    breaks = seq(2014, 2024, 1)) +
  
  # scale_y_reverse() + 
  
  scale_y_continuous(expand = c(0.01,0.01))  +
  scale_fill_manual('Continent' ,values = region_colours,
                    labels = str_to_title) +
  
  annotate("rect", 
           xmin = shading_intervals, 
           xmax = shading_intervals_end, 
           ymin = -Inf, ymax = Inf, alpha = 0.2, fill = "grey") +
  
  geom_tippoint(aes(fill = collection_regionname),
                shape = 21,
                size = 2,
                size = 0.75) +
  
  scale_fill_manual(values = region_colours) +
  scale_colour_manual(values = region_colours) +
  
  
  new_scale_fill()+
  new_scale_colour()+
  
  
  geom_fruit(geom = geom_tile,
             mapping = aes(fill = PB2),
             #width = 4,
             #colour = "white",
             #pwidth = 1.2,
             #offset = 0.03
  ) + 
  scale_fill_distiller(palette = 'Purples', direction  = 1, limits = c(1, 30)) + 
  
  new_scale_fill()+
  geom_fruit(geom = geom_tile,
             mapping = aes(fill = PB1),
             #width = 4,
             #colour = "white",
             #pwidth = 1.2,
             offset = 0.03
  ) + 
  #scale_fill_paletteer_c("ggthemes::Orange")+
  scale_fill_distiller(palette = 'Purples', direction  = 1, limits = c(1, 30)) + 
  
  
  
  new_scale_fill()+
  geom_fruit(geom = geom_tile,
             mapping = aes(fill = PA),
             #width = 4,
             # colour = "white",
             #pwidth = 1.2,
             offset = 0.03
  ) + 
  #scale_fill_paletteer_c("ggthemes::Orange-Gold")+
  scale_fill_distiller(palette = 'Purples', direction  = 1, limits = c(1, 30)) + 
  
  new_scale_fill()+
  geom_fruit(geom = geom_tile,
             mapping = aes(fill = HA),
             #width = 4,
             # colour = "white",
             #pwidth = 1.2,
             offset = 0.03
  ) + 
  #scale_fill_paletteer_c("ggthemes::Green-Gold")+
  scale_fill_distiller(palette = 'Purples', direction  = 1, limits = c(1, 30)) + 
  
  
  new_scale_fill()+
  geom_fruit(geom = geom_tile,
             mapping = aes(fill = N),
             #width = 4,
             #colour = "white",
             #pwidth = 1.2,
             offset = 0.03
  ) + 
  #scale_fill_paletteer_c("ggthemes::Green")+
  scale_fill_distiller(palette = 'Purples', direction  = 1, limits = c(1, 30)) + 
  
  new_scale_fill()+
  geom_fruit(geom = geom_tile,
             mapping = aes(fill = NP),
             #width = 4,
             # colour = "white",
             #pwidth = 1.2,
             offset = 0.03
  ) + 
  #scale_fill_paletteer_c("ggthemes::Blue-Green Sequential")+
  scale_fill_distiller(palette = 'Purples', direction  = 1, limits = c(1, 30)) + 
  new_scale_fill()+
  geom_fruit(geom = geom_tile,
             mapping = aes(fill = N),
             #width = 4,
             # colour = "white",
             #pwidth = 1.2,
             offset = 0.03
  ) + 
  #scale_fill_paletteer_c("ggthemes::Blue")+
  scale_fill_distiller(palette = 'Purples', direction  = 1, limits = c(1, 30)) + 
  
  new_scale_fill()+
  geom_fruit(geom = geom_tile,
             mapping = aes(fill = NS),
             #width = 4,
             # colour = "white",
             #pwidth = 1.2,
             offset = 0.03
  ) + 
  #scale_fill_paletteer_c("ggthemes::Purple")+
  scale_fill_distiller(palette = 'Purples', direction  = 1, limits = c(1, 30)) + 
  
  theme(legend.position = 'none',
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 9))


plt_1b <- combined_data %>%
  filter(segment == 'ha') %>%
  mutate(collection_regionname = case_when(grepl('europe', collection_regionname) ~ 'europe',
                                           grepl('africa', collection_regionname) ~ 'africa',
                                           grepl('asia', collection_regionname) ~ 'asia',
                                           grepl('(central|northern) america|caribbean', collection_regionname) ~ 'central & northern america',
                                           grepl('south america|southern ocean', collection_regionname) ~ 'south america',
                                           grepl('australia|melanesia', collection_regionname) ~ 'australasia',
                                           .default = collection_regionname)) %>%
  summarise(n = n(), .by = collection_regionname) %>%
  mutate(collection_regionname = reorder(collection_regionname, n)) %>%
  ggplot() +
  geom_bar(aes(x = collection_regionname,
               y = n, 
               fill = collection_regionname), 
           stat = 'identity',
           colour = 'black') + 
  scale_fill_manual(values = region_colours) +
  scale_y_continuous(expand = c(0,0), 'Reassortants (N)') + 
  scale_x_discrete(expand = c(0.2,0), labels = function(x) str_to_title(x) %>% str_wrap(., width = 10), '') + 
  
  theme_classic() + 
  theme(
    panel.background = element_rect(fill = "transparent", colour = NA), 
    plot.background = element_rect(fill = "transparent", colour = NA),
    legend.position = 'none',
    strip.background = element_blank(),
    strip.text = element_blank(),
    axis.text = element_text(size = 8),
    axis.title = element_text(size = 9),
  )


plt_1c <- combined_data %>% 
  filter(segment == 'ha') %>%
  mutate(host_simplifiedhost = case_when(.default = host_simplifiedhost,
                                         host_simplifiedhost == 'anseriformes-wild+other-bird' ~ 'unknown',
                                         host_simplifiedhost == 'other-bird' ~ 'unknown',
                                         host_simplifiedhost == 'environment' ~ 'unknown',
                                         grepl('galliformes', host_simplifiedhost) ~ 'galliformes')) %>%
  summarise(n = n(), .by = host_simplifiedhost) %>%
  mutate(host_simplifiedhost = reorder(host_simplifiedhost, n)) %>%
  ggplot() +
  geom_bar(aes(x = host_simplifiedhost,y = n, fill = host_simplifiedhost), colour= 'black', stat = 'identity') + 
  scale_fill_manual(values = host_colours) +
  scale_y_continuous(expand = c(0,0), 'Reassortants (N)') + 
  scale_x_discrete(expand = c(0.2,0), labels = c('unknown' = 'Unknown',
                                                 'anseriformes-wild' = 'Anser. Wild' ,
                                                 'charadriiformes-wild' = 'Charad. Wild',
                                                 'galliformes' = 'Gall.',
                                                 'mammal' = 'Mammal' ) %>% 
                     str_wrap(., width = 10), '') + 
  
  theme_classic() + 
  theme(
    panel.background = element_rect(fill = "transparent", colour = NA), 
    plot.background = element_rect(fill = "transparent", colour = NA),
    legend.position = 'none',
    strip.background = element_blank(),
    strip.text = element_blank(),
    axis.text = element_text(size = 8),
    axis.title = element_text(size = 9),
  )



plt_1d <- combined_data %>% 
  filter(segment == 'ha') %>%
  
  mutate( 
    collection_month = date_decimal(TMRCA) %>% format(., "%m") %>% as.integer(),
    collection_season = case_when(collection_month %in% c(12,1,2) ~ 'overwintering', 
                                  collection_month %in% c(3,4,5)  ~ 'spring migration', # Rename to spring migration
                                  collection_month %in% c(6,7,8)  ~ 'breeding', 
                                  collection_month %in% c(9,10,11)  ~ 'autumn migration' )) %>%
  summarise(n = n(), .by = collection_season) %>%
  mutate(collection_season = reorder(collection_season, n)) %>%
  ggplot() +
  geom_bar(aes(x = collection_season,y = n, fill = collection_season), colour= 'black', stat = 'identity') + 
  scale_fill_brewer(palette = 'Accent')+
  scale_y_continuous(expand = c(0,0), 'Reassortants (N)') + 
  scale_x_discrete(expand = c(0.2,0), labels = function(x) str_to_title(x) %>% str_wrap(., width = 10), '') + 
  
  theme_classic() + 
  theme(
    panel.background = element_rect(fill = "transparent", colour = NA), 
    plot.background = element_rect(fill = "transparent", colour = NA),
    legend.position = 'none',
    strip.background = element_blank(),
    strip.text = element_blank(),
    axis.text = element_text(size = 8),
    axis.title = element_text(size = 9),
  )


# N~ Reassortants over time (inset)
plot_data <- count_data %>%
  #drop_na(n_reassortants) %>%
  mutate(collection_year_quarter = quarter(collection_datemonth, type = 'year.quarter' )) %>%
  replace_na(list(n_reassortants = 0 )) %>%
  drop_na(n_reassortants) %>%
  summarise(n_reassortants = sum(n_reassortants), .by = c(collection_year_quarter, collection_regionname)) 

years <- 2014:2024
quarters <- c(0.1, 0.2, 0.3, 0.4)

plot_data <- expand_grid(collection_year_quarter =  as.vector(outer(years, quarters, `+`)),
                         collection_regionname = unique(plot_data$collection_regionname)) %>%
  left_join(plot_data) %>%
  mutate(collection_year_dec = case_when(grepl('1$', as.numeric(collection_year_quarter)) ~ floor(as.numeric(collection_year_quarter)) + 0.125,
                                         grepl('2$', as.numeric(collection_year_quarter)) ~ floor(as.numeric(collection_year_quarter)) + 0.375,
                                         grepl('3$', as.numeric(collection_year_quarter)) ~ floor(as.numeric(collection_year_quarter)) + 0.625,
                                         grepl('4$', as.numeric(collection_year_quarter)) ~ floor(as.numeric(collection_year_quarter)) + 0.875)) %>%
  
  filter(between(collection_year_quarter, 2014.3, 2024.4)) %>%
  drop_na(collection_regionname) %>%
  mutate(collection_year_quarter = as.factor(collection_year_quarter))

inset_count <- plot_data %>%
  ggplot(aes(x = collection_year_dec, 
             y = n_reassortants,
             fill = collection_regionname, 
             colour = collection_regionname)) + 
  geom_bar(position = 'stack',
           stat = 'identity',
           alpha = 0.7)  +
  scale_x_continuous(
    limits = c(2016, 2025),
    'Time (Years)',
    breaks = seq(2014, 2024, 1)) +
  #scale_x_discrete(breaks = levels(plot_data$collection_year_quarter)[c(T, rep(F, 3))],
  #  labels =  gsub('\\.\\d+$', '', levels(plot_data$collection_year_quarter)[c(T, rep(F, 3))] ),
  # expand = c(0,0),
  #'Date') + 
  #annotate("rect", 
  #  xmin = shading_intervals, 
  #  xmax = shading_intervals_end, 
  #  ymin = -Inf, ymax = Inf, alpha = 0.2, fill = "grey") +
  scale_y_continuous('Count',expand = c(0,0)) +
  scale_fill_manual('Continent' ,values = region_colours,
                    labels = str_to_title) +
  scale_colour_manual('Continent' , values = region_colours,
                      labels = str_to_title) +
  theme_classic()+
  theme(axis.text = element_text(size = 8),
        axis.title = element_text(size = 9),
        panel.background = element_rect(fill = "transparent", colour = NA), 
        plot.background = element_rect(fill = "transparent", colour = NA),
        
        legend.position = 'inside',
        legend.justification.inside = c(0.1, 1),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 9),
        legend.position.inside = c(0.1,0.95),
        legend.background = element_rect(fill = "transparent", colour = NA)
  )



plot_1a_with_inset <-ggdraw() +
  draw_plot(plt_1a) +
  draw_plot(inset_count, x = 0.01, y = .7, width = .5, height = .3) 

plt_lower <- cowplot::plot_grid(plt_1b,
                                plt_1c, 
                                plt_1d, 
                                ncol = 3, 
                                align = 'hv',
                                axis = 'tblr', 
                                labels = c('B', 'C', 'D'), 
                                label_size = 9)

cowplot::plot_grid(plot_a_with_inset, 
                   plt_lower, 
                   align = 'h', 
                   axis = 'lr',
                   nrow = 2, 
                   rel_heights = c(0.8, 0.2),
                   labels = c('A', ''), 
                   label_size = 9)


############################################## WRITE ###############################################

ggsave('~/Downloads/flu_plots/figure1_new.jpeg', height = 30, width = 25, units = 'cm', dpi = 360)



############################################## END #################################################
####################################################################################################
####################################################################################################