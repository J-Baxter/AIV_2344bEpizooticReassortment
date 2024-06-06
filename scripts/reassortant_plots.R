# Reassortant summaries
require(extrafont)
library(tidyverse)
library(sf)
library(rnaturalearth)
#https://www.fontsquirrel.com/fonts/latin-modern-sans
font_import(pattern = "lmsans10*") 


my_theme <- theme_classic(base_family = "LM Sans 10")+
  theme(
    #text = element_text(size=10),
    axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0), size = 8),
    axis.title.y = element_text(margin = margin(t = 0, r = 5, b = 0, l = 0), size = 8),
    axis.text = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0), size = 7),
    strip.text  = element_text(size = 8),
    legend.text = element_text(size = 7),
    legend.title = element_text(size = 8),
    legend.position = 'none', 
    #panel.spacing = unit(2, "lines"), 
    strip.background = element_blank()
  )

summary_tbl <- read_csv('summary_reassortant_metadata.csv') %>%
  filter(!clade %in% c("EA_nonGsGD","Am_nonGsGD")) %>%
  rename(reassortant_class = group) %>%
  mutate(reassortant_class = case_when(cluster_profile == '1_1_1_1_1_1_1_1' ~ 'major', .default = reassortant_class)) 
  
  #mutate(reassortant_class = case_when(
   # Length_Between_First_Last_Sample > 1 & number_Conti > 1 & number_hostType >1 ~ 'dominant',
   # Length_Between_First_Last_Sample >= 0.5 & Length_Between_First_Last_Sample <= 1.4 & number_Conti >= 1 & number_hostType >=1 & Num_Sequence >= 3 ~ 'major', 
   # Length_Between_First_Last_Sample < 0.5 & number_Conti == 1~ 'minor',
   # .default = NA)) 

plt_1a <- summary_tbl %>%
  select(c(Continent_of_Earliest_Date)) %>%
  ggplot() +
  geom_bar(aes(x = Continent_of_Earliest_Date),
           stat = "count") +
  scale_y_continuous(expand = c(0,0),
                     'Number of uniqe reassortants') +
  scale_x_discrete(expand = c(0.2,0.2),
                   'Region of earliest sequence date') +
  my_theme 

plt_1b <- summary_tbl %>%
  select(c(Continent_of_Earliest_Date, reassortant_class, cluster_profile)) %>%
  mutate(reassortant_class = factor(reassortant_class, levels = c('minor', 'major', 'dominant'))) %>%
  ggplot() +
  geom_bar(aes(x = Continent_of_Earliest_Date,
               fill = reassortant_class),
           stat = "count", 
           position = 'fill') +
  scale_y_continuous(expand = c(0,0),
                     'Proportion of uniqe reassortants') +
  scale_x_discrete(expand = c(0.2,0.2),
                   'Region of earliest sequence date') +
  scale_fill_brewer(palette = 'Dark2',
                    'Reassortant class') +
  my_theme +
  theme(legend.position = 'right')


frequencyplots <- cowplot::plot_grid(plt_1a, plt_1b, nrow = 1, rel_widths = c(0.4, 0.6), align = 'hv', labels='AUTO', label_size =10)

ggsave('frequencyplots.svg', device= 'svg',  height = 110, width = 240, units = 'mm')
Sys.sleep(0.5)
frequencyplots
dev.off()

# Frequency of reassortants 


reassortant_classification <- summary_tbl%>% 
  filter(reassortant_class == 'dominant') %>%
  mutate(reassortant_class = case_when(cluster_profile == '1_1_1_1_1_1_1_1' ~ 'major', .default = reassortant_class)) %>%
  select(reassortant_class, cluster_profile, colour) %>%
  distinct()

# import formmated reassortant data, 
formatted_reassortantdata <- read_csv('./2023Dec01/metadata/h5_metadata_global_6280_update_formatted.csv') %>%
  left_join(., reassortant_classification, by = join_by(cluster.profile == cluster_profile)) %>%
  filter(reassortant_class == 'dominant') 
# left join with classification

reassortant_counts <- formatted_reassortantdata %>%
  filter(!is.na(collection.date)) %>%
  mutate(flu_season = ifelse(month(collection.date) < 9, paste0(year(collection.date) - 1, "-", year(collection.date)), paste0(year(collection.date), "-", year(collection.date) + 1))) %>%
  summarise(n = n(), .by = c(cluster.profile, collection.country.code, flu_season, colour)) %>%
  mutate(collection.country.code = toupper(collection.country.code))

reassortants <- ne_countries(returnclass = "sf") %>% select(adm0_iso) %>%
  expand_grid(., reassortant_counts %>% select(cluster.profile) %>% distinct()) %>%
  expand_grid(., reassortant_counts %>% select(flu_season) %>% distinct()) %>%
  left_join(reassortant_counts, by = join_by(adm0_iso ==collection.country.code, cluster.profile == cluster.profile, flu_season == flu_season))


create_sequential_palette <- function(base_colour, n) {
  # Create a function to generate a sequential color palette
  ramp_func <- colorRampPalette(c("white", base_colour, 'black'))
  
  # Generate the color palette
  palette <- ramp_func(n)
  
  return(palette)
}

MapPlot <- function(data){
  col <- base::unique(data$colour[!is.na(data$colour)])

  map <- ne_countries(returnclass = "sf") %>%
    left_join(data, by = join_by(adm0_iso), multiple = "all")
   

  plot <- ggplot() +
    geom_sf(data = map, aes(fill = n), colour = 'white', linewidth = 0.1) + 
    coord_sf(ylim = c(-55,80), expand = FALSE) +
    my_theme + 
    scale_fill_gradientn(colors = create_sequential_palette(col, 11)[3:10] ,
                         limits = c(0,300),
                         breaks = seq(0,300,by=75),
                         na.value=create_sequential_palette(col, 10)[2]) +
    facet_grid(cols = vars(flu_season),
               drop = F)+
    my_theme +
    theme(legend.position = 'right',
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.line = element_blank()) + 
    guides(fill = guide_coloursteps(barwidth = 0.76, barheight = 4))
  
  return(plot)
}

mapplotlist <- reassortants %>%
  group_split(cluster.profile) %>%
  lapply(., MapPlot)

profile_list <- reassortants %>%
  group_split(cluster.profile) %>%
  lapply(. ,function(x) x %>% 
           pull(cluster.profile) %>%
           unique()) %>%
  unlist() %>%
  paste(collapse = ',')


ggsave('reassortantmapsnew.svg', device= 'svg',  height = 250, width = 320, units = 'mm')
Sys.sleep(0.5)
cowplot::plot_grid(plotlist =  mapplotlist, 
                   nrow = 6,
                   labels = c("1_1_2_1_1_1_1_1",
                              "2_1_1_1_1_1_1_1","3_2_3_1_1_2_1_2",
                              "4_1_1_1_2_1_1_8","5_3_1_1_2_1_1_4", 
                              "6_1_5_2_1_3_1_2"), 
                   label_size = 9)
Sys.sleep(0.5)
dev.off()
