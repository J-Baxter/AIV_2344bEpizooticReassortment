library(tidyverse)

summary_data <- read_csv('2024Aug18/treedata_extractions/2024-09-20_combined_data.csv')
metadata_date <- lapply(metadatafiles, read_csv, col_types = cols(collection_tipdate = col_character())) %>%
  lapply(. , function(x) x %>% filter(!is.na(cluster_profile)) %>% select(isolate_id, collection_date, cluster_profile)) %>%
  bind_rows() %>%
  distinct() %>%
  filter(!is.na(collection_date)) %>%
  summarise(first_date = min(collection_date),
            laste_date = max(collection_date),
            .by = cluster_profile)


summary_dataplus <- summary_data %>%
  left_join(metadata_date)


ggplot(summary_dataplus) +
  geom_segment(aes(
    y = cluster_profile, 
    yend = cluster_profile,
                x = first_date,
                xend = laste_date+5,
    colour = group)) +
  facet_grid(rows = vars(Continent_of_Earliest_Date),
             scales = 'free_y',
             space = 'free_y',
             switch = 'y')+ 
  theme_classic() +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank())

ggplot(summary_dataplus) +
  geom_histogram(aes(x = Length_Between_First_Last_Sample)) +
  facet_grid(rows = vars(Continent_of_Earliest_Date),
             scales = 'free_y') 



read_csv('2024Aug18/treedata_extractions/2024-09-20_combined_data.csv') %>%
  select(c(evoRate, persist.time, group2, collection_regionname, segment)) %>%
  #filter(segment %in% c('ha', 'pb2')) %>%
  mutate(collection_regionname = case_when(grepl('europe', 
                                                 collection_regionname) ~ 'europe',
                                           grepl('africa',
                                                 collection_regionname) ~ 'africa',
                                           grepl('asia', 
                                                 collection_regionname) ~ 'asia',
                                           grepl('(central|northern) america', 
                                                 collection_regionname) ~ 'central & northern america',
                                           .default = collection_regionname)) %>%
  ggplot()+
  geom_point(aes(x = evoRate, y = persist.time, colour = group2)) + 
  #facet_wrap(~ collection_regionname) + 
  facet_wrap(~ segment) + 
  theme_minimal(base_size = 8)+
  scale_colour_brewer(palette = 'Set1')+
  scale_x_continuous('Evolutionary Rate', limits = c(0, 0.01)) +
  scale_y_continuous('Persistence Time')