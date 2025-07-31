################### Plot 1: 'Persistence Plot' for Dominant Reassortants ################### 

combined_data %>%
  filter(group2 == 'dominant') %>%
  mutate(cluster_profile = gsub('_', '', cluster_profile)) %>%
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
  scale_colour_brewer(palette = 'Dark2') +
  coord_flip() + 
  
  # Graphical
  #facet_grid(rows = vars(segment)) + 
  theme_minimal() 


################### Plot 2: Region MRCA ################### 
combined_data %>%
  mutate(collection_regionname = case_when(grepl('europe', collection_regionname) ~ 'europe',
                                           grepl('africa', collection_regionname) ~ 'africa',
                                           grepl('asia', collection_regionname) ~ 'asia',
                                           grepl('(central|northern) america', collection_regionname) ~ 'central & northern america',
                                           .default = collection_regionname
  )) %>%
  filter(!is.na(collection_regionname)) %>%
  ggplot(aes(x = collection_regionname,
             fill = group2)) +
  
  # Geom objects
  geom_bar(position = "fill") +
  scale_y_continuous('Number of Reassortants', labels = scales::percent) +
  
  # Scales
  scale_x_discrete('Region of Origin', labels = function(x) str_wrap(x, width = 20) %>% str_to_title())+
  scale_fill_brewer('Reassortant Class', direction = -1) +
  
  
  
  # Graphical
  facet_grid(rows = vars(segment)) + 
  theme_minimal() 


################### Plot 2: Host MRCA - in progress ################### 
combined_data %>%
  filter(!is.na(host_simplifiedhost)) %>%
  ggplot(aes(x = host_simplifiedhost,
             fill = group2)) +
  
  # Geom objects
  geom_bar(position = "fill") +
  scale_y_continuous('Proportion of Reassortants') +
  
  # Scales
  scale_x_discrete('Host of Origin',
                   labels = function(x) str_wrap(x, width = 20) %>% str_to_title())+
  scale_fill_brewer('Host' ,
                    palette = 'Blues') +
  
  # Graphical
  facet_grid(rows = vars(segment)) + 
  theme_minimal() 


################### Plot 3:  Evo-rate density plot (all) ################### 

combined_data %>%
  ggplot(aes(x = evoRate,
             fill = group2,
             y = after_stat(scaled),
             colour = group2)) +
  
  # Geom objects
  geom_density(alpha = 0.7) +
  scale_y_continuous('Probability Density') +
  
  # Scales
  scale_x_continuous('Evolutionary Rate')+
  scale_fill_brewer('Reassortant Class', direction = -1) +
  scale_colour_brewer('Reassortant Class', direction = -1) +
  
  
  
  
  # Graphical
  #facet_grid(rows = vars(segment)) + 
  theme_minimal() 


################### Plot 4:  Evo-rate plot (dominant) ################### 

combined_data %>%
  filter(cluster_profile %in% c("3_2_3_1_3_2_1_2",
                                "2_1_1_1_1_1_1_1",
                                "1_1_1_1_1_1_1_1",
                                "1_1_1_1_1_1_1_1A",
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
  mutate(cluster_profile = gsub('_', '', cluster_profile)) %>%
  
  ggplot(aes(x = cluster_profile,
             ymin = evoRate_95lower,
             ymax = evoRate_95upper,
             y = evoRate,
             group = cluster_profile,
             colour = segment)) +
  
  # Geom objects
  geom_point(position = position_dodge2(
    width = 0.6,
    preserve = "total",
    padding = 0.1,
    reverse = FALSE
  )) + 
  
  geom_linerange(position = position_dodge2(
    width = 0.6,
    preserve = "total",
    padding = 0.1,
    reverse = FALSE
  )) + 
  # Scales
  scale_y_continuous('Evolutionary Rate')+
  scale_x_discrete('Cluster Profile') + 
  scale_color_brewer(palette = 'Dark2') +
  coord_cartesian(ylim = c(0, 0.02)) +
  
  
  # Graphical
  #facet_grid(rows = vars(segment)) + 
  theme_minimal() 
