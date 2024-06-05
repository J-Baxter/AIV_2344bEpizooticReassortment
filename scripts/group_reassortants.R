GroupReassortants <- function(data){
  
  summary <- data %>% 
    as_tibble() %>%
    mutate(date_year = as.double(date_year),
           location_1 = case_when(location_1 == 'Antarctica' ~ 'South America', 
                                  .default = location_1))  %>%
    # Format and filter dates
    mutate(date_format = case_when(
      str_count(date, "-") < 1 ~'yyyy',
      str_count(date, "-") == 1 ~'yyyy-mm',
      str_count(date, "-") == 2 ~'yyyy-mm-dd'
    )) %>%
    filter(date_format != 'yyyy') %>%
    split(~date_format) %>% 
    map_at("yyyy-mm-dd",  
           ~ mutate(.x,
                    date_parsed = ymd(date),
                    date_ymd = format(date_parsed, '%Y-%m-%d'), 
                    date_dec = decimal_date(date_parsed),
                    date_ym = format(date_parsed, '%Y-%m'),
                    date_y = format(date_parsed, '%Y')))  %>%
    map_at("yyyy-mm",
           ~ mutate(.x,
                    date_parsed = ym(date),
                    date_ymd = NA_character_, 
                    date_dec = decimal_date(date_parsed),
                    date_ym = format(date_parsed,'%Y-%m') ,
                    date_y = format(date_parsed, '%Y'))) %>%
    list_rbind() %>%
    
    # Count unique continents/host type (wild, domestic , mammal, other) for each cluster
    group_by(cluster_profile) %>%
    mutate(number_Conti = length(unique(location_1)),
           number_hostType = length(unique(host_class))) %>%
    
    # Summarise across all data
    group_by(cluster_label,clade,number_Conti,number_hostType, .add = T) %>%
    dplyr::summarise(
      Num_Sequence = n(),
      Range_of_Host = paste(unique(host_order), collapse = ", "),
      Range_of_HostType = paste(unique(host_class), collapse = ", "),
      Host_Range_Per_Continent = paste(unique(paste(host_class, location_1, sep = " - ")), collapse = ", "),
      Continent_of_Earliest_Date = location_1[which.min(date_dec)],
      Host_of_Earliest_Sample_Date = host_class[which.min(date_dec)],
      Range_of_Continent = paste(unique(location_1), collapse = ", "),
      Length_Between_First_Last_Sample = max(date_dec) - min(date_dec),
      number_Mammal = sum(host_class == "Mammal"),
      prop_Mammal = sum(host_class == "Mammal")/n()
    ) %>%
    ungroup() %>%
    
    # Infer profile-specific groups
    mutate(group = case_when(
      Length_Between_First_Last_Sample <0.5 & prop_Mammal < 0.01 & number_Mammal <=2~ 'minor',
      Length_Between_First_Last_Sample > 1 & number_Conti > 1 & number_hostType > 1 ~ 'dominant',
      TRUE ~ 'major'
    )) 
  
  return(summary)
}