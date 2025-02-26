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
# User functions


############################################## DATA ################################################
combined_data <- read_csv('./2024Aug18/treedata_extractions/2024-09-20_combined_data.csv')
summary_data <- read_csv('./2024Aug18/treedata_extractions/summary_reassortant_metadata_20240904.csv') %>%
  select(-c(cluster_label,
            clade)) 

meta <- read_csv('./2024-09-09_meta.csv') 

will_tree <- read.beast('./global_subsample/h52344b_ha_s1.mcc.trees')

new_tree <- read.beast('./2025Jan06/globalsubsample/ha_global_subsample_traits_mcc.tree')
timetree <- read.newick('./2025Feb26/globalsubsample/2025-02-26_clock/rerooted.newick')

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


# Sequences
sequences_month <- meta %>%  
  drop_na(cluster_profile) %>%
  select(starts_with('collection_date'),
         collection_regionname) %>%
  mutate(collection_regionname = case_when(grepl('europe', collection_regionname) ~ 'europe',
                                           grepl('africa', collection_regionname) ~ 'africa',
                                           grepl('asia', collection_regionname) ~ 'asia',
                                           grepl('(central|northern) america|caribbean', collection_regionname) ~ 'central & northern america',
                                           grepl('south america|southern ocean', collection_regionname) ~ 'south america',
                                           grepl('australia|melanesia', collection_regionname) ~ 'australasia',
                                           .default = collection_regionname)) %>%
  group_by(collection_datemonth, collection_regionname) %>%
  summarise(n_sequences = n()) %>%
  ungroup() %>%
  mutate(collection_datemonth = ymd(paste0(collection_datemonth, '-01'))) %>%
  drop_na(collection_regionname,collection_datemonth)


# Format WOAH data and estimate the minimum number of cases
woah_minimuminferredcases <- woah_hpai %>%
  select(Year,
         Semester,
         `World region`,
         Country,
         `Animal Category`,
         `Serotype/Subtype/Genotype`,
         `New outbreaks`, 
         Susceptible, 
         `Measuring units`,
         Cases,
         Deaths) %>%
  
  # replace dashes with NA
  mutate(across(!Year, ~ case_when(grepl('^-$', .x) ~ NA_character_,
                                   .default = .x))) %>%
  
  # format biannual periods
  mutate(date_start = case_when(grepl('Jan-Jun', Semester) ~ paste0(Year, '-01-01'),
                                grepl('Jul-Dec', Semester) ~ paste0(Year, '-07-01'),
                                .default = NA_character_),
         date_end = case_when(grepl('Jan-Jun', Semester) ~ paste0(Year, '-06-30'),
                              grepl('Jul-Dec', Semester) ~ paste0(Year, '-12-31'),
                              .default = NA_character_)) %>%
  
  # format country for continent assignment
  mutate(Country = case_when(Country ==  "Cote D'Ivoire" ~ "Côte d'Ivoire",
                             Country == "Congo (Dem. Rep. of the)" ~ "Dem. Rep. Congo",
                             Country == "China (People's Rep. of)"~ "China",
                             Country == "Chinese Taipei" ~ "Taiwan",
                             Country == "Korea (Dem People's Rep. of)" ~"North Korea",
                             Country == "Korea (Rep. of)" ~"South Korea",
                             Country == "Dominican (Rep.)" ~ "Dominican Rep.",
                             grepl( "Falkland Islands \\(Malvinas\\)|South Georgia and the South Sandwich Islands", Country) ~ "Falkland Is.",
                             Country ==  "Türkiye (Rep. of)" ~ "Turkey",
                             Country == "Hong Kong" ~ 'China',                                 
                             Country == "Bosnia and Herzegovina" ~  "Bosnia and Herz." ,                 
                             Country == "Czech Republic" ~ "Czechia",                               
                             Country == "Faeroe Islands" ~ "Denmark",                              
                             Country == "Reunion" ~ 'Madagascar',                                      
                             .default = Country))%>%
  # Group by continent
  left_join(ref, by = join_by(Country == name)) %>%
  mutate(animal_category  = str_to_lower(`Animal Category`)) %>%   
  mutate(across(any_of(c('New outbreaks', 'Susceptible', 'Cases', 'Deaths')), ~as.numeric(.x))) %>%
  replace_na(list('New outbreaks' = 0,
                  'Susceptible' = 0,
                  'Cases' = 0,
                  'Deaths' = 0)) %>%
  
  # Filter only H5NX
  filter(grepl('^[Hh]5', `Serotype/Subtype/Genotype`)) %>%
  
  # Infer the minimum number of 'cases' This is mostly due to USA seemingly incapable of reporting
  # their data in a similar manner to the rest of the world.
  mutate(inferred_minimum_cases = case_when(Cases >= Deaths ~ Cases, 
                                            Cases < Deaths ~ Deaths,
                                            .default = Cases)) %>%
  
  rename(collection_regionname = continent) %>%
  mutate(collection_regionname = str_to_lower(collection_regionname)) %>%
  mutate(collection_regionname = case_when(grepl('europe', collection_regionname) ~ 'europe',
                                           grepl('africa', collection_regionname) ~ 'africa',
                                           grepl('asia', collection_regionname) ~ 'asia',
                                           grepl('(central|northern|north) america|caribbean', collection_regionname) ~  'central & northern america',
                                           grepl('south america|southern ocean|antarctica', collection_regionname) ~ 'south america',
                                           grepl('australia|melanesia|oceania', collection_regionname) ~ 'australasia',
                                           .default = collection_regionname))


# Group WOAH minimum cases by month
woah_minimuminferredcases_monthly <- woah_minimuminferredcases %>%
  
  summarise(sum_IMC = sum(inferred_minimum_cases, na.rm = TRUE), .by = c(date_start, 
                                                                         date_end,
                                                                         collection_regionname))  %>%
  mutate(across(starts_with('date'), ~ymd(.x))) %>%
  mutate(interval = interval(date_start, date_end)) %>%
  rename(woah_cases = sum_IMC)


# Group FAO cases by month
fao_hpai_monthly <- fao_hpai %>% 
  filter(pathogenicity == 'HPAI') %>%
  filter(grepl('H5', serotype)) %>%
  mutate(observation_month = format.Date(observation_date, 
                                         format = '%Y-%m')) %>%
  group_by(observation_month, 
           collection_regionname) %>%
  summarise(n_cases = n()) %>%
  ungroup() %>%
  mutate(observation_month = ymd(paste0(observation_month, '-01'))) %>%
  filter(observation_month> ymd('2015-06-01')) %>% 
  rename(collection_datemonth = observation_month)  %>%
  mutate(date_start = floor_date(collection_datemonth,
                                 unit = 'month'),
         date_end = ceiling_date(collection_datemonth, 
                                 unit = 'month') - days(1)) %>%
  mutate(collection_regionname = case_when(collection_regionname == 'north and central america' ~ 'central & northern america', 
                                           .default = collection_regionname)) %>%
  rename(fao_cases = n_cases)


# Join everything together and calculated scaled estimates (not used)
all_casedata_monthly <- woah_minimuminferredcases_monthly %>%
  interval_left_join(fao_hpai_monthly, by = c('date_start', 'date_end')) %>%
  filter(collection_regionname.x == collection_regionname.y) %>%
  select(-ends_with('y')) %>%
  rename_with(~gsub('\\.x', '', .x)) %>%
  
  # set groupings
  group_by(collection_regionname, interval) %>%
  mutate(scale_factor = woah_cases/sum(fao_cases)) %>%
  ungroup() %>%
  rowwise() %>%
  mutate(scaled_estimate = fao_cases*scale_factor) %>%
  filter(!grepl('australasia', collection_regionname)) %>%
  filter(collection_datemonth >= as_date('2019-01-01'))


# Plot
plt_1a <- all_casedata_monthly %>%
  as_tibble() %>%  filter(date_start < as_date('2024-06-01'))  %>%
  ggplot() + 

  geom_bar(data = sequences_month, 
           aes(x = collection_datemonth, 
               y = n_sequences*7000, 
               colour = collection_regionname,
               fill = collection_regionname), 
           alpha = 0.7, 
           stat = 'identity',
           position = 'stack') +
  
  geom_line(aes(x = collection_datemonth, 
                y = woah_cases, 
                colour= collection_regionname),
            linetype = 'dashed') +
  
  scale_fill_manual(values = region_colours) +
  scale_colour_manual(values = region_colours) +
  theme_classic() + 
  scale_y_continuous(expand = c(0.01, 0),
                     'WOAH Minimum Cases (n)',
                     sec.axis = sec_axis(transform = ~ ./7000, 
                                         'GISAID Whole Genomes (n)')) + 
  scale_x_date('Date') + 
  
  facet_wrap(~collection_regionname,  
             ncol = 6) + 
  
  geom_text(data = cbind.data.frame(collection_regionname = unique(test_3$collection_regionname)), 
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
  scale_x_date(limits = as_date(c('2019-01-01', '2024-05-01')), breaks = '1 year', date_labels = "%Y", 'Date') + 
  scale_fill_manual(values = host_colours) +
  scale_colour_manual(values = host_colours) +
  scale_y_continuous('GISAID Whole Genomes (n)' ,
                     labels = abs,
                     breaks = seq(-250, 250, by = 50),
                     limits = c(-250, 250)) + 
  theme_classic() + 
  theme(legend.position = 'none',
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 10))



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
  scale_x_date(limits = as_date(c('2019-01-01', '2024-05-01')), breaks = '1 year', date_labels = "%Y", 'Date') + 
  scale_fill_brewer(palette = 'OrRd', direction = -1) +
  scale_colour_brewer(palette = 'OrRd', direction = -1) +
  scale_y_continuous('GISAID Whole Genomes (n)' ,
                     labels = abs,
                     breaks = seq(-250, 250, by = 50),
                     limits = c(-250, 250)) + 
  theme_classic() + 
  theme(legend.position = 'none',
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 10))



##################################
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



plt_1rightplots <- align_plots(plt_1a, plt_1c, plt_1d, align = 'h', axis = 'r')

plt_1rightpanel <- plot_grid(
                             plt_1rightplots[[2]],  
                             plt_1rightplots[[3]], 
                             combineLegend,
                             labels = c( 'C', 'D', ''), 
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
plt_1left <- timetree %>%
  as.treedata() %>%
  mutate(isolate_id = str_extract(label, "EPI_ISL_(china_){0,1}\\d+[^.|]*")) %>%
 # mutate(date = str_extract(label, "\\d{4}-.*")) %>%
  left_join(reassortant_profiles) %>%
  
  ggtree(mrsd = "2024-03-18") + ################################# NEEDS UPDATING #####################
  theme_tree2(#base_family = "LM Sans 10",
              #plot.margin = unit(c(1,1,1,1), units = "cm"),
              axis.text = element_text(size = 8),
              axis.title = element_text(size = 10)
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
  scale_fill_paletteer_c("ggthemes::Red") +
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