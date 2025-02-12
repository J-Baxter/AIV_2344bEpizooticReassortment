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
# User functions


############################################## DATA ################################################
combined_data <- read_csv('./2024Aug18/treedata_extractions/2024-09-20_combined_data.csv')
summary_data <- read_csv('./2024Aug18/treedata_extractions/summary_reassortant_metadata_20240904.csv') %>%
  select(-c(cluster_label,
            clade)) 

meta <- read_csv('./2024-09-09_meta.csv') 

will_tree <- read.beast('./global_subsample/h52344b_ha_s1.mcc.trees')

new_tree <- read.beast('./2025Jan06/globalsubsample/ha_global_subsample_traits_mcc.tree')

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
 'galliformes-wild' = '#33a02c',
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

all <- meta %>%
  mutate(collection_regionname = case_when(grepl('europe', collection_regionname) ~ 'europe',
                                           grepl('africa', collection_regionname) ~ 'africa',
                                           grepl('asia', collection_regionname) ~ 'asia',
                                           grepl('(central|northern) america', collection_regionname) ~ 'central & northern america',
                                           grepl('south america|southern ocean', collection_regionname) ~ 'south america',
                                           grepl('australia', collection_regionname) ~ 'australasia',
                                           .default = collection_regionname)) %>%
  drop_na(collection_regionname) %>%
  dplyr::select(collection_regionname, collection_dateyear) %>%
  tidyr::expand(collection_regionname, collection_dateyear = full_seq(collection_dateyear,1))

  
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
  scale_y_continuous('Cumulative Frequency', expand = c(0,0)) + 
  scale_x_continuous('Collection Year', expand = c(0,0)) + 
  scale_fill_manual(values = region_colours) + 
  global_theme +
  theme(legend.position = 'inside',
        legend.position.inside = c(0.05, 0.95),
        legend.title = element_blank(),
        legend.justification = c("left", "top")) 
  


# Number of Sequences ~ Host Type (cumulatively)
all <- meta %>%
  drop_na(host_simplifiedhost) %>%
  dplyr::select(host_simplifiedhost, collection_dateyear) %>%
  tidyr::expand(host_simplifiedhost, collection_dateyear = full_seq(collection_dateyear,1))


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
  scale_y_continuous('Cumulative Frequency', expand = c(0,0)) + 
  scale_x_continuous('Collection Year', expand = c(0,0)) + 
  scale_fill_manual(values = host_colours) +
  global_theme +
  theme(legend.position = 'inside',
        legend.position.inside = c(0.05, 0.95),
        legend.title = element_blank(),
        legend.justification = c("left", "top"))


# Reassortants ~ Time
plt_1d <- combined_data %>% 
  filter(segment == 'ha') %>%
  select(TMRCA) %>%
  drop_na(TMRCA) %>%
  mutate(tmrca_date = date_decimal(TMRCA)) %>%
  mutate(tmrca_yearmonth = round_date(tmrca_date, 'month') %>% as_date(.)) %>%
  ggplot() + 
  geom_bar(aes(x = tmrca_yearmonth)) +
  scale_x_date('Collection Year', expand = c(0,0), limits = as_date(c('2016-01-01','2024-06-01'))) + 
  scale_y_continuous('Frequency', expand = c(0,0)) + 
  global_theme


# Reassortants ~ Continent (TMRCA)
plt_1e <- combined_data %>% 
  filter(segment == 'ha') %>%
  select(collection_regionname) %>%
  mutate(collection_regionname = case_when(grepl('europe', collection_regionname) ~ 'europe',
                                           grepl('africa', collection_regionname) ~ 'africa',
                                           grepl('asia', collection_regionname) ~ 'asia',
                                           grepl('(central|northern) america', collection_regionname) ~ 'central & northern america',
                                           grepl('south america|southern ocean', collection_regionname) ~ 'south america',
                                           grepl('australia', collection_regionname) ~ 'australasia',
                                           .default = collection_regionname)) %>%
  drop_na(collection_regionname) %>%
  summarise(n = n(), .by = c(collection_regionname)) %>%
  ggplot() + 
  geom_bar(aes(x = collection_regionname, y = n, fill = collection_regionname), stat = 'identity') +
  scale_fill_manual(values = region_colours) + 
  scale_x_discrete('Earliest Region', expand = c(0.15,0.1),labels = c('europe' = 'EUR',
                                                                      'africa' = 'AFR',
                                                                      'asia' = 'ASIA',
                                                                      'central & northern america' = 'AMR')) + 
  scale_y_continuous('Frequency', expand = c(0,0)) + 
  global_theme+ 
  theme(legend.position = 'none')
  


# Reassortants ~ Host Type (TMRCA)
plt_1f <- summary_data %>% 
  mutate(Host_of_Earliest_Sample_Date = case_when(grepl('Anseriformes domestic', Host_of_Earliest_Sample_Date) ~ 'anseriformes-domestic',
                                                  grepl('Anseriformes wild', Host_of_Earliest_Sample_Date) ~ 'anseriformes-wild',
                                                  grepl('Avian other', Host_of_Earliest_Sample_Date) ~ 'other-bird',
                                                grepl('Charadriiformes', Host_of_Earliest_Sample_Date) ~ 'charadriiformes-wild',
                                                grepl('Environment', Host_of_Earliest_Sample_Date) ~ 'environment',
                                                grepl('Galliformes', Host_of_Earliest_Sample_Date) ~ 'galliformes-domestic',
                                                grepl('Mammal', Host_of_Earliest_Sample_Date) ~ 'mammal',
                                                .default = Host_of_Earliest_Sample_Date
  )) %>%
  drop_na(Host_of_Earliest_Sample_Date) %>%
  summarise(n = n(), .by = c(Host_of_Earliest_Sample_Date)) %>%
  ggplot() + 
  geom_bar(aes(x = Host_of_Earliest_Sample_Date, y = n, fill = Host_of_Earliest_Sample_Date), stat = 'identity') +
  scale_fill_manual(values = host_colours) + 
  scale_x_discrete('Earliest Host', expand = c(0.1,0.1) ,labels = c('anseriformes-domestic' = 'ANS-dom',
                                                                    'anseriformes-wild' = 'ANS-wild',
                                                                    'charadriiformes-wild' = 'CHAR-wild',
                                                                    'galliformes-domestic' = 'GAL-dom',
                                                                    'galliformes-wild' = 'GAL-wild',
                                                                    'environment' = 'ENV',
                                                                    'mammal' = 'MAM',
                                                                    'other-bird' = 'OTHER')) + 
  scale_y_continuous('Frequency', expand = c(0,0)) + 
  global_theme + 
  theme(legend.position = 'none')

# Combine to make right panel
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





# Tree Panel
# A time-scaled MCC tree of HA
# column geoms for reassortant profile
# node = red for inferred changes
reassortant_profiles <- meta %>% 
  dplyr::select(isolate_id, cluster_profile) %>%
  drop_na(cluster_profile) %>%
  separate_wider_delim(cluster_profile, '_', 
                       names = c('PB2', 'PB1', 'PA', 'HA', 'NP', 'N', 'M', 'NS')) %>%
  mutate(across(-1, .fns = ~ as.double(.x)))

# Use Will's for now, to be replaced by thorney tree
plt_1left <- new_tree %>%
  mutate(isolate_id = gsub('\\|.*','', label)) %>%
  left_join(reassortant_profiles) %>%
  
  ggtree(mrsd = '2024-04') + ################################# NEEDS UPDATING #####################
  theme_tree2(#base_family = "LM Sans 10",
              #plot.margin = unit(c(1,1,1,1), units = "cm"),
              axis.text.x = element_text(family = "LM Sans", size = 8),
              axis.title.x = element_text(family = "LM Sans", size = 9)
              ) +
  
  scale_x_continuous(
    #limits = c(2000, 2023),
    'Time',
    breaks = seq(2016, 2024, 1)) +

  
  geom_fruit(geom = geom_tile,
             mapping = aes(fill = PB2),
             #width = 4,
             #colour = "white",
            #pwidth = 1.2,
             #offset = 0.03
             ) + 
  scale_fill_paletteer_c("ggthemes::Red")+
  #scale_fill_distiller(palette = 'Greys', direction  = 1) + 
  
  new_scale_fill()+
  geom_fruit(geom = geom_tile,
             mapping = aes(fill = PB1),
             #width = 4,
             #colour = "white",
             #pwidth = 1.2,
              offset = 0.03
  ) + 
  scale_fill_paletteer_c("ggthemes::Orange")+
  #scale_fill_distiller(palette = 'Greys', direction  = 1) + 
  
  new_scale_fill()+
  geom_fruit(geom = geom_tile,
             mapping = aes(fill = PA),
             #width = 4,
            # colour = "white",
             #pwidth = 1.2,
             offset = 0.03
  ) + 
  scale_fill_paletteer_c("ggthemes::Orange-Gold")+
  #scale_fill_distiller(palette = 'Greys', direction  = 1) + 
  
  new_scale_fill()+
  geom_fruit(geom = geom_tile,
             mapping = aes(fill = HA),
             #width = 4,
            # colour = "white",
             #pwidth = 1.2,
             offset = 0.03
  ) + 
  scale_fill_paletteer_c("ggthemes::Green-Gold")+
  #scale_fill_distiller(palette = 'Greys', direction  = 1) + 
  
  new_scale_fill()+
  geom_fruit(geom = geom_tile,
             mapping = aes(fill = N),
             #width = 4,
             #colour = "white",
             #pwidth = 1.2,
             offset = 0.03
  ) + 
  scale_fill_paletteer_c("ggthemes::Green")+
  #scale_fill_distiller(palette = 'Greys', direction  = 1) + 
  
  new_scale_fill()+
  geom_fruit(geom = geom_tile,
             mapping = aes(fill = NP),
             #width = 4,
            # colour = "white",
             #pwidth = 1.2,
             offset = 0.03
  ) + 
  scale_fill_paletteer_c("ggthemes::Blue-Green Sequential")+
  #scale_fill_distiller(palette = 'Greys', direction  = 1) + 
  new_scale_fill()+
  geom_fruit(geom = geom_tile,
             mapping = aes(fill = N),
             #width = 4,
            # colour = "white",
             #pwidth = 1.2,
             offset = 0.03
  ) + 
  scale_fill_paletteer_c("ggthemes::Blue")+
  #scale_fill_distiller(palette = 'Greys', direction  = 1) + 
  
  new_scale_fill()+
  geom_fruit(geom = geom_tile,
             mapping = aes(fill = NS),
             #width = 4,
            # colour = "white",
             #pwidth = 1.2,
             offset = 0.03
  ) + 
  scale_fill_paletteer_c("ggthemes::Purple")+
  #scale_fill_distiller(palette = 'Greys', direction  = 1) + 
  theme(legend.position = 'none' )


plt_1 <- plot_grid(plt_1left,
                   plt_1right, 
                        labels = c('A', ''), 
                   label_size = 12,
                   rel_widths = c(0.8, 1),
                        ncol = 2)
plt_1
############################################## WRITE ###############################################

ggsave('~/Downloads/figure1.jpeg', height = 30, width = 40, units = 'cm', dpi = 360)



############################################## END #################################################
####################################################################################################
####################################################################################################