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
source('./scripts/figure_scripts/global_theme.R')
# User functions


############################################## DATA ################################################
combined_data <- read_csv('./2024Aug18/treedata_extractions/2024-09-20_combined_data.csv')
summary_data <- read_csv('./2024Aug18/treedata_extractions/summary_reassortant_metadata_20240904.csv') %>%
  dplyr::select(-c(cluster_label,
            clade)) 

meta <- read_csv('./2024-09-09_meta.csv') 

new_tree <- read.beast('./2025Feb26/globalsubsample/ha_global_SRD06_relaxLn_constant_mcc.tree')

hpai_cases <- read_csv('~/Downloads/overview-raw-data_202502241440.csv') 
woah_hpai <- read_csv('~/Downloads/Quantitative data 2025-02-25.csv')

ref <- ne_countries() %>%
  as_tibble()  %>% 
  dplyr::select(name, continent)

############################################## MAIN ################################################
# Plot
plt_1a <- new_tree %>%
  p <- new_tree %>%
  mutate(isolate_id = str_extract(label, "EPI_ISL_(china_){0,1}\\d+[^.|]*")) %>%
  # mutate(date = str_extract(label, "\\d{4}-.*")) %>%
  left_join(reassortant_profiles) %>%
  left_join(continent) %>%
  #left_join(mrca_nodes) %>%
  
  ggtree(mrsd = "2024-03-18") + 
  theme_tree2(#base_family = "LM Sans 10",
    #plot.margin = unit(c(1,1,1,1), units = "cm"),
    axis.text.x = element_text(size = 7),
    axis.title.x = element_text(size = 8)
  ) +
  
  #theme_classic() + 
  
  #scale_x_continuous(
  # limits = c(2014.9, Inf),
  # 'Time (Years)',
  # breaks = seq(2014, 2024, 1)) +
  
  # scale_y_reverse() + 
  
  
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
  theme(legend.position = 'none',
        axis.text = element_text(size = 7),
        axis.title = element_text(size = 8))


p_data <- reassortant_profiles %>% left_join(meta %>% dplyr::select(isolate_id, tipnames)) %>% dplyr::select(-isolate_id)  %>% as.data.frame()
rownames(p_data) <- p_data$tipnames

all_colors <- brewer.pal(9, "Purples")
color_subset <- colorRampPalette(all_colors)(100)[30:100]  # middle 50%

plt_1a <- gheatmap(p, p_data %>% dplyr::select(-tipnames),  width=0.4, offset = 0.1,
                   color = NA,
         colnames=FALSE, legend_title="genotype") +
  scale_x_ggtree() +
  #scale_fill_gradientn(colours = color_subset, na.value = 'white') + 
  scale_fill_distiller(palette = 'Purples', direction  = 1, na.value = 'white', transform = 'log2') + 
  scale_y_continuous(expand = c(0.01,0.01), limits = c(0, 1000))  +
  theme(legend.position = 'none',
        axis.text = element_text(size = 7),
        axis.title = element_text(size = 8))


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
    axis.text = element_text(size = 7),
    axis.title = element_text(size = 8),
  )

plt1d_colours <- c(
  'anseriformes-domestic' = '#a6cee3',
  'anser. wild' = '#1f78b4',
  'galliformes-domestic' = '#b2df8a',
  'gall.' = '#33a02c',
  'mammal' = '#fb9a99',
  'human' = '#e31a1c',
  'charad. wild' = '#fdbf6f',
  'other-bird' = '#ff7f00',
  'unknown' = '#cab2d6',
  'environment' = '#6a3d9a')

plt_1c <- combined_data %>% 
  filter(segment == 'ha') %>%
  mutate(host_simplifiedhost = case_when(.default = host_simplifiedhost,
                                         host_simplifiedhost == 'anseriformes-wild+other-bird' ~ 'unknown',
                                         host_simplifiedhost == 'other-bird' ~ 'unknown',
                                         host_simplifiedhost == 'environment' ~ 'unknown',
                                         host_simplifiedhost == 'anseriformes-wild' ~ 'anser. wild',
                                         host_simplifiedhost == 'charadriiformes-wild' ~ 'charad. wild',
                                         grepl('galliformes', host_simplifiedhost) ~ 'gall.')) %>%
  summarise(n = n(), .by = host_simplifiedhost) %>%
  mutate(host_simplifiedhost = reorder(host_simplifiedhost, n)) %>%
  ggplot() +
  geom_bar(aes(x = host_simplifiedhost,y = n, fill = host_simplifiedhost), colour= 'black', stat = 'identity') + 
  scale_fill_manual(values = plt1d_colours) +
  scale_y_continuous(expand = c(0,0), 'Reassortants (N)') + 
  scale_x_discrete(expand = c(0.2,0), labels = function(x) str_to_title(x) %>% str_wrap(., width = 10), '') + 
  
  theme_classic() + 
  theme(
    panel.background = element_rect(fill = "transparent", colour = NA), 
    plot.background = element_rect(fill = "transparent", colour = NA),
    legend.position = 'none',
    strip.background = element_blank(),
    strip.text = element_blank(),
    axis.text = element_text(size = 7),
    axis.title = element_text(size = 8),
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
    axis.text = element_text(size = 7),
    axis.title = element_text(size = 8),
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
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 8),
        legend.position.inside = c(0.1,0.95),
        legend.background = element_rect(fill = "transparent", colour = NA)
  )


inset_seqs <-  ggplot(data = sequences_month, 
                      aes(x = collection_datemonth, 
                          y = n_sequences,
                          coloru = collection_regionname,
                          fill = collection_regionname)) +
  geom_area(position = 'stack') + 
  scale_y_continuous('Count',expand = c(0,0)) +
  xlab('Sampling Date') + 
  scale_colour_manual('Continent' , values = region_colours,
                    labels = str_to_title) +
  scale_fill_manual('Continent' , values = region_colours,
                      labels = str_to_title) +
  theme_classic()+
  theme(axis.text = element_text(size = 7),
        axis.title = element_text(size = 8),
        panel.background = element_rect(fill = "transparent", colour = NA), 
        plot.background = element_rect(fill = "transparent", colour = NA),
        
        legend.position = 'inside',
        legend.justification.inside = c(0, 1),
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 8),
        legend.position.inside = c(0,0.99),
        legend.key.size = unit(0.75,"line"),
        legend.background = element_rect(fill = "transparent", colour = NA)
  )

  
  
plot_1a_with_inset <- ggdraw() +
  draw_plot(plt_1a) +
  draw_plot(inset_seqs, x = 0, y = .75, width = .45, height = .25) 

plt_lower <- cowplot::plot_grid(plt_1b,
                                plt_1c, 
                                plt_1d, 
                                ncol = 3, 
                                align = 'hv',
                                axis = 'tblr', 
                                labels = c('B', 'C', 'D'), 
                                label_size = 9)

cowplot::plot_grid(plot_1a_with_inset, 
                   plt_lower, 
                   align = 'h', 
                   axis = 'lr',
                   nrow = 2, 
                   rel_heights = c(0.8, 0.2),
                   labels = c('A', ''), 
                   label_size = 9)


############################################## WRITE ###############################################
#210x297 ()
ggsave('~/Downloads/flu_plots/figure1_new.pdf', height = 257, width = 206, units = 'mm', dpi = 360)



############################################## END #################################################
####################################################################################################
####################################################################################################