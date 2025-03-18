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
plt_1a <- data %>%
  as_tibble() %>%  
  filter(collection_datemonth > as_date('2018-06-01'))  %>%
  mutate(across(where(is.numeric), .fns = ~replace_na(.x, 0))) %>%
  
  ggplot() + 
  
  geom_bar(data = sequences_month, 
           aes(x = collection_datemonth, 
               y = n_sequences*20000, 
               colour = collection_regionname,
               fill = collection_regionname), 
           alpha = 0.7, 
           stat = 'identity',
           position = 'stack') +
  
  geom_line(aes(x = collection_datemonth, 
                y = woah_susceptibles, 
                colour= collection_regionname),
            linetype = 'dashed') +
  
  scale_fill_manual(values = region_colours) +
  scale_colour_manual(values = region_colours) +
  theme_classic() + 
  scale_y_continuous(expand = c(0.01, 0),
                     'WOAH Susceptibles (n)',
                     sec.axis = sec_axis(transform = ~ ./20000, 
                                         'GISAID Whole Genomes (n)')) + 
  scale_x_date(limits = as_date(c('2019-01-01', '2024-05-01')),
               breaks = '2 year', 
               date_labels = "%Y", 'Date',
               expand = c(0,0)) + 
  
  facet_wrap(~collection_regionname,  
             ncol = 6) + 
  
  geom_text(data = cbind.data.frame(collection_regionname = unique(all_casedata_monthly$collection_regionname)), 
            aes(label = str_wrap(str_to_title(collection_regionname), width = 20)),
            y = Inf, 
            x = as_date('2019-01-01'), 
            size = 3,
            fontface = 'bold',
            vjust = 'top',
            hjust = 'left') +
  theme(
    legend.position = 'none',
    strip.background = element_blank(),
    strip.text = element_blank(),
    axis.text = element_text(size = 8),
    axis.title = element_text(size = 10)
  )


plt_1a


##################################
# Plot 1C
sequences_host_month <- meta %>%  
  drop_na(cluster_profile) %>%
  select(starts_with('collection_date'),
         host_simplifiedhost,
         collection_regionname) %>%
  mutate(collection_regionname = case_when(grepl('europe', collection_regionname) ~ 'europe',
                                           grepl('africa', collection_regionname) ~ 'africa',
                                           grepl('asia', collection_regionname) ~ 'asia',
                                           grepl('(central|northern) america|caribbean', collection_regionname) ~ 'central & northern america',
                                           grepl('south america|southern ocean', collection_regionname) ~ 'south america',
                                           grepl('australia|melanesia', collection_regionname) ~ 'australasia',
                                           .default = collection_regionname)) %>%
  group_by(collection_datemonth, host_simplifiedhost) %>%
  summarise(n_sequences = n()) %>%
  ungroup() %>%
  mutate(collection_datemonth = ymd(paste0(collection_datemonth, '-01'))) %>%
  drop_na(host_simplifiedhost,collection_datemonth)


plt_1c <-ggplot(sequences_host_month) + 
  geom_stream(aes(x = collection_datemonth, 
                  y = n_sequences, 
                  colour = host_simplifiedhost,
                  fill = host_simplifiedhost),
              alpha = 0.7) + 
  scale_x_date(limits = as_date(c('2016-01-01', '2024-05-01')),
               breaks = '1 year', 
               date_labels = "%Y", 'Date',
               expand = c(0,0)) + 
  scale_fill_manual('Host Order', values = host_colours) +
  scale_colour_manual('Host Order', values = host_colours) +
  scale_y_continuous('GISAID Whole Genomes (n)' ,
                     labels = abs,
                     breaks = seq(-225, 225, by = 75),
                     limits = c(-230, 230)) + 
  theme_classic() + 
  theme(axis.text = element_text(size = 8),
        axis.title = element_text(size = 10),
        legend.position = 'inside',
        legend.justification.inside = c(0, 1),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 10),
        legend.position.inside = c(0.01,1.05),
        legend.background = element_blank())



##################################
# Plot 1D
sequences_subtype_month <- meta %>%  
  drop_na(cluster_profile) %>%
  select(starts_with('collection_date'),
         virus_subtype,
         collection_regionname) %>%
  mutate(collection_regionname = case_when(grepl('europe', collection_regionname) ~ 'europe',
                                           grepl('africa', collection_regionname) ~ 'africa',
                                           grepl('asia', collection_regionname) ~ 'asia',
                                           grepl('(central|northern) america|caribbean', collection_regionname) ~ 'central & northern america',
                                           grepl('south america|southern ocean', collection_regionname) ~ 'south america',
                                           grepl('australia|melanesia', collection_regionname) ~ 'australasia',
                                           .default = collection_regionname)) %>%
  group_by(collection_datemonth, virus_subtype) %>%
  summarise(n_sequences = n()) %>%
  ungroup() %>%
  mutate(collection_datemonth = ymd(paste0(collection_datemonth, '-01'))) %>%
  drop_na(virus_subtype,collection_datemonth)


plt_1d <- ggplot(sequences_subtype_month) + 
  geom_stream(aes(x = collection_datemonth, 
                  y = n_sequences, 
                  colour = virus_subtype,
                  fill = virus_subtype),
              alpha = 0.7) + 
  scale_x_date(limits = as_date(c('2016-01-01', '2024-05-01')),
               breaks = '1 year', 
               date_labels = "%Y", 'Date',
               expand = c(0,0)) + 
  scale_fill_brewer('Subtype', palette = 'OrRd', direction = -1) +
  scale_colour_brewer('Subtype', palette = 'OrRd', direction = -1) +
  scale_y_continuous('GISAID Whole Genomes (n)' ,
                     labels = abs,
                     breaks = seq(-225, 225, by = 75),
                     limits = c(-230, 230)) + 
  theme_classic() + 
  theme(
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 10),
        legend.position = 'inside',
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 10),
        legend.justification.inside = c(0, 1),
        legend.position.inside = c(0.01,1.05),
        legend.background = element_blank())


plt_1e <- data %>%
drop_na(n_reassortants) %>%
  ggplot() + 
  geom_histogram(aes(x = collection_datemonth, 
                     fill = collection_regionname, 
                     colour = collection_regionname),
                 bins = 45,
                 alpha = 0.7)  +
  scale_x_date(limits = as_date(c('2016-01-01', '2024-05-01')),
               breaks = '1 year', 
               date_labels = "%Y", 'Date',
               expand = c(0,0)) + 
  
  scale_y_continuous('Count',expand = c(0,0)) +
  scale_fill_manual('Continent' ,values = region_colours,
                    labels = str_to_title) +
  scale_colour_manual('Continent' , values = region_colours,
                      labels = str_to_title) +
  global_theme + 
  theme(axis.text = element_text(size = 8),
        axis.title = element_text(size = 10),
        legend.position = 'inside',
        legend.justification.inside = c(0, 1),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 10),
        legend.position.inside = c(0.01,1.05),
        legend.background = element_blank())


##################################
# reassortant MRCA (will run for 10-20 mins)
traits <- new_tree@phylo$tip.label %>%
  str_extract(., "(\\d_)+[^.|]*") 

new_tree@phylo$edge.length <- ifelse(new_tree@phylo$edge.length < 0, 0, new_tree@phylo$edge.length)

ace_test <- ace(traits,
                new_tree@phylo,
                type = 'discrete',
                method = 'ML',
                model = 'ARD')

ace_test$lik.anc


# extract likelihood for each state and each node
anc_lik <- ace_test$lik.anc

# infer most likely state at each node and 
anc_state <- apply(anc_lik,
                   1,
                   function(x) names(x)[which.max(x)]) %>%
  #format output
  enframe(., 
          name = 'node', 
          value = 'cluster_profile') %>%
  mutate(node = as.integer(node))


mrca_nodes <- anc_state %>%
  left_join(as_tibble(new_tree) %>% dplyr::select(node, height_median)) %>%
  group_by(cluster_profile) %>%
  slice_max(height_median, n =1, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(tmrca_node = '1') %>%
  dplyr::select(-height_median)



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

shading_intervals <- seq(2015, 2023, by = 2)
shading_intervals_end <- c(shading_intervals[-1] - 1, 2024)

continent <- meta %>% 
  dplyr::select(isolate_id, collection_regionname) %>%
  mutate(collection_regionname = case_when(grepl('europe', collection_regionname) ~ 'europe',
                                           grepl('africa', collection_regionname) ~ 'africa',
                                           grepl('asia', collection_regionname) ~ 'asia',
                                           grepl('(central|northern) america|caribbean', collection_regionname) ~ 'central & northern america',
                                           grepl('south america|southern ocean', collection_regionname) ~ 'south america',
                                           grepl('australia|melanesia', collection_regionname) ~ 'australasia',
                                           .default = collection_regionname))
plt_1left <- new_tree %>%
  mutate(isolate_id = str_extract(label, "EPI_ISL_(china_){0,1}\\d+[^.|]*")) %>%
 # mutate(date = str_extract(label, "\\d{4}-.*")) %>%
  left_join(reassortant_profiles) %>%
  left_join(continent) %>%
  left_join(mrca_nodes) %>%
  
  ggtree(mrsd = "2024-03-18") + 
  theme_tree2(#base_family = "LM Sans 10",
              #plot.margin = unit(c(1,1,1,1), units = "cm"),
              axis.text.x = element_text(size = 8),
              axis.title.x = element_text(size = 10)
              ) +
  
  scale_x_continuous(
   # limits = c(2014.9, Inf),
    'Time (Years)',
    breaks = seq(2015, 2024, 1)) +
  
  
  annotate("rect", 
           xmin = shading_intervals, 
           xmax = shading_intervals_end, 
           ymin = -Inf, ymax = Inf, alpha = 0.2, fill = "grey") +
  
  #geom_nodepoint(aes(colour = tmrca_node,
                   #  shape = tmrca_node)) +
 # scale_shape_manual(values = c("1" = 18),
                     #guide = 'none') +
  #scale_colour_manual(values = c("1" = 'red'), 
                     # guide = 'none') + 
  #new_scale_colour()+
  
  
  geom_tippoint(aes(fill = collection_regionname,
                    colour = collection_regionname),
                shape = 21,
                size = 0.75, 
                alpha = 0.9) +
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
  #scale_fill_distiller(palette = 'Purples', direction = -1)
  #scale_fill_paletteer_c("ggthemes::Red") +
  scale_fill_distiller(palette = 'Purples', direction  = 1) + 
  
  new_scale_fill()+
  geom_fruit(geom = geom_tile,
             mapping = aes(fill = PB1),
             #width = 4,
             #colour = "white",
             #pwidth = 1.2,
              offset = 0.03
  ) + 
  #scale_fill_paletteer_c("ggthemes::Orange")+
  scale_fill_distiller(palette = 'Purples', direction  = 1) + 
  
  new_scale_fill()+
  geom_fruit(geom = geom_tile,
             mapping = aes(fill = PA),
             #width = 4,
            # colour = "white",
             #pwidth = 1.2,
             offset = 0.03
  ) + 
  #scale_fill_paletteer_c("ggthemes::Orange-Gold")+
  scale_fill_distiller(palette = 'Purples', direction  = 1) + 
  
  new_scale_fill()+
  geom_fruit(geom = geom_tile,
             mapping = aes(fill = HA),
             #width = 4,
            # colour = "white",
             #pwidth = 1.2,
             offset = 0.03
  ) + 
  #scale_fill_paletteer_c("ggthemes::Green-Gold")+
  scale_fill_distiller(palette = 'Purples', direction  = 1) + 
  
  new_scale_fill()+
  geom_fruit(geom = geom_tile,
             mapping = aes(fill = N),
             #width = 4,
             #colour = "white",
             #pwidth = 1.2,
             offset = 0.03
  ) + 
  #scale_fill_paletteer_c("ggthemes::Green")+
  scale_fill_distiller(palette = 'Purples', direction  = 1) + 
  
  new_scale_fill()+
  geom_fruit(geom = geom_tile,
             mapping = aes(fill = NP),
             #width = 4,
            # colour = "white",
             #pwidth = 1.2,
             offset = 0.03
  ) + 
  #scale_fill_paletteer_c("ggthemes::Blue-Green Sequential")+
  scale_fill_distiller(palette = 'Purples', direction  = 1) + 
  new_scale_fill()+
  geom_fruit(geom = geom_tile,
             mapping = aes(fill = N),
             #width = 4,
            # colour = "white",
             #pwidth = 1.2,
             offset = 0.03
  ) + 
  #scale_fill_paletteer_c("ggthemes::Blue")+
  scale_fill_distiller(palette = 'Purples', direction  = 1) + 
  
  new_scale_fill()+
  geom_fruit(geom = geom_tile,
             mapping = aes(fill = NS),
             #width = 4,
            # colour = "white",
             #pwidth = 1.2,
             offset = 0.03
  ) + 
  #scale_fill_paletteer_c("ggthemes::Purple")+
  scale_fill_distiller(palette = 'Purples', direction  = 1) + 
  theme(legend.position = 'none' )


# Combine to make right panel
legend  <- cowplot::get_legend(plt_1c + theme(legend.position = 'right',
                                              legend.justification = c(0.5,1),
                                              legend.text = element_text(size = 8),
                                              legend.title = element_text(size = 10)))
legend1 <- cowplot::get_legend(plt_1d + theme(legend.position = 'right',
                                              legend.justification = c(0.5,1),
                                              legend.text = element_text(size = 8),
                                              legend.title = element_text(size = 10)))
combineLegend <- cowplot::plot_grid(
  legend,
  legend1,
  align = 'v',
  axis = 't',
  ncol= 2)


###### COMBINE
plt_1rightplots <- align_plots(plt_1a, plt_1c, plt_1d, plt_1e, align = 'h', axis = 'r')

plt_1rightpanel <- plot_grid(
  plt_1rightplots[[2]],  
  plt_1rightplots[[3]], 
  #combineLegend,
  plt_1rightplots[[4]], 
  labels = c( 'C', 'D', 'E'), 
  nrow = 3,
  label_size = 10)

plt_1lower <- plot_grid(plt_1left, plt_1rightpanel,
                        labels = c('B', ''),
                        ncol = 2,
                        label_size = 10)

plt_1 <- plot_grid(plt_1rightplots[[1]], plt_1lower,
                   labels = c('A', ''), 
                   rel_heights = c(0.25, 1),
                   label_size = 10,
                   nrow = 2)

plt_1
############################################## WRITE ###############################################

ggsave('~/Downloads/flu_plots/figure1.jpeg', height = 30, width = 25, units = 'cm', dpi = 360)



############################################## END #################################################
####################################################################################################
####################################################################################################