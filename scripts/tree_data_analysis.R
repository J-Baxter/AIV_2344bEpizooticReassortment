# Dependencies
library(tidyverse)

# Import Tree Data Extractions
csv_paths <- list.files('./2024Aug18/treedata_extractions',
                        pattern = '.csv',
                        full.names = TRUE)

csv_names <- gsub('.*treedata_extractions/|.csv$',
                  '',
                  csv_paths)

csv_data <- lapply(csv_paths, read_csv) %>%
  setNames(csv_names)

# Expand to separate objects
for(i in 1:length(csv_data)){
  assign(names(csv_data)[i],as_tibble(csv_data[[i]]))
}


################### Combined Tibble ###################

combined_data <- RegionalTre_mrca_stats_disp_jumps_combined %>%
  
  #remove dominant cluster profiles
  filter(!Rcode %in% c("3_2_3_1_3_2_1_2",
                       "2_1_1_1_1_1_1_1",
                       "1_1_1_1_1_1_1_1",
                       "2_1_2_1_1_1_1_1",
                       "1_1_2_1_1_1_1_1",
                       "1_6_2_1_1_1_1_1",
                       "1_1_4_1_4_1_1_4",
                       "2_6_1_1_6_1_1_1",
                       "2_1_6_1_1_4_1_1",
                       "7_1_5_2_1_3_1_2",
                       "4_3_1_1_2_1_1_3",
                       "5_1_1_1_2_1_1_3",
                       "5_4_9_1_2_1_1_1")) %>%
  
  # remove deprecated columns
  select(-c(cluster_number, 
            collection_countryname,
            virus_subtype,
            prob_cluster_number, 
            prob_collection_countryname,
            prob_virus_subtype)) %>%
  
  # bind results from dominant trees
  bind_rows(domainantTre_mrca_stats_disp_jumps_combined) %>%
  
  # bind summary table                             
  left_join(summary_reassortant_metadata_20240902 %>%
              select(c(cluster_profile, 
                       col2,
                       group2,
                       Length_Between_First_Last_Sample)),
            by = join_by(Rcode == cluster_profile)) %>%
  
  # split labels to extract segments
  separate_wider_delim(data_name, '_', names = c('segment', 'deprecated')) %>%
  mutate(cluster_profile = case_when(deprecated == '11111111A'& Rcode == '1_1_1_1_1_1_1_1' ~ '1_1_1_1_1_1_1_1A',
                                     .default = Rcode)) %>%
  select(- c(deprecated,Rcode)) %>%
  mutate(segment = gsub('n[:0-9:]', 'nx', segment)) %>%
  relocate(cluster_profile, segment)


################### Plot 1: 'Persistence Plot' for Dominant Reassortants ################### 

combined_data %>%
  filter(group2 == 'dominant') %>%
  ggplot() +
  
  # Geom objects
  geom_linerange(aes(x = cluster_profile,
                     ymin = youngestTip.time-Length_Between_First_Last_Sample, 
                     ymax = youngestTip.time,
                     colour = segment), 
                 linewidth = 1.5, 
                 position = position_dodge(width = 0.75)) +
  
  geom_linerange(aes(x = cluster_profile,
                     ymin = TMRCA, 
                     ymax = youngestTip.time,
                     colour = segment), 
                 linewidth = 2, 
                 alpha = 0.5,
                 position = position_dodge(width = 0.75)) +
  # Scales
  scale_x_discrete('Reassortant Profile')+
  scale_y_continuous()+
  scale_colour_uchicago(palette = 'light') +
  coord_flip() + 
  
  # Graphical
  #facet_grid(rows = vars(segment)) + 
  theme_minimal() 



#### persistence time distributions by country ####
region_persistencetime <- RegionalTre_mrca_stats_disp_jumps_combined %>%
  select(c(data_name,
           youngestTip.time,
           TMRCA,
           Rcode)) %>%
  left_join(summary_reassortant_metadata_20240902,
            by = join_by(Rcode == cluster_profile)) %>%
  select(c(data_name,
           TMRCA,
           Rcode,
           Length_Between_First_Last_Sample,
           youngestTip.time,
           group2)) %>%
  rename(cluster_group = group2) %>%
  separate_wider_delim(data_name, '_', names = c('segment', 'region')) %>%
  mutate(segment = gsub('n[:0-9:]', 'nx', segment)) %>%
  filter(cluster_group != 'dominant')


ggplot(dominant_persistencetime) +
  # Geom objects
  geom_linerange(aes(x = cluster_profile,
                     ymin = youngestTip.time-Length_Between_First_Last_Sample, 
                     ymax = youngestTip.time,
                     colour = segment), 
                 linewidth = 1.5, 
                 position = position_dodge(width = 0.75)) +
  
  geom_linerange(aes(x = cluster_profile,
                     ymin = TMRCA, 
                     ymax = youngestTip.time,
                     colour = segment), 
                 linewidth = 2, 
                 alpha = 0.5,
                 position = position_dodge(width = 0.75)) +
  # Scales
  scale_x_discrete('Reassortant Profile')+
  scale_y_continuous()+
  scale_colour_uchicago(palette = 'light') +
  coord_flip() + 
  #facet_grid(rows = vars(segment)) + 
  theme_minimal() 