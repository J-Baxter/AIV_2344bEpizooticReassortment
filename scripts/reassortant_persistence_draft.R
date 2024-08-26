library(tidyverse)

summary_data <- read_csv('2024-06-05_reassortant_summary.csv')
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