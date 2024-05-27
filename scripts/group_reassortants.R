GroupReassortants <- function(data){
  
  summary <- data %>%
    group_by(cluster_profile) %>%
    
    # Count unique continents/host type (wild, domestic , mammal, other) for each cluster
    mutate(number_Conti = length(unique(location_1)),
           number_hostType = length(unique(host_class))) %>%
    group_by(cluster_label,clade,number_Conti,number_hostType, .add = T) %>%
    
    # Summarise across all data
    dplyr::summarise(
      Num_Sequence = n(),
      Range_of_Host = paste(unique(host_order), collapse = ", "),
      Range_of_HostType = paste(unique(host_class), collapse = ", "),
      Host_Range_Per_Continent = paste(unique(paste(host_class, location_1, sep = " - ")), collapse = ", "),
      Continent_of_Earliest_Date = location_1[which.min(date_frac)],
      Host_of_Earliest_Sample_Date = host_class[which.min(date_frac)],
      Range_of_Continent = paste(unique(location_1), collapse = ", "),
      Length_Between_First_Last_Sample = max(date_frac) - min(date_frac)
    ) %>%
    ungroup() %>%
    
    # Infer profile-specific groups
    mutate(group = case_when(
      Num_Sequence < 3 & Length_Between_First_Last_Sample <0.01 ~ 'minor',
      Length_Between_First_Last_Sample > 1 & number_Conti > 1 & number_hostType > 1 ~ 'dominant',
      TRUE ~ 'major'
    )) 
  
  return(summary)
}