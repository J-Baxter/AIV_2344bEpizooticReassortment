count_data <- combined_data %>%
  filter(segment == 'ha') %>%
  mutate(collection_monthyear = date_decimal(TMRCA) %>% 
           format(., "%Y-%m")) %>%
  mutate(collection_regionname = case_when(grepl('europe', 
                                                 collection_regionname) ~ 'europe',
                                           grepl('africa',
                                                 collection_regionname) ~ 'africa',
                                           grepl('asia', 
                                                 collection_regionname) ~ 'asia',
                                           grepl('(central|northern) america', 
                                                 collection_regionname) ~ 'central & northern america',
                                           .default = collection_regionname)) %>%
  drop_na(collection_regionname) %>%
  select(collection_regionname,
         collection_monthyear, 
         cluster_profile,
         group2) %>% # include all month-years over collection period to generate zer countr
  summarise(n_reassortants = n_distinct(cluster_profile), 
            .by = c(collection_monthyear,
                    collection_regionname,
                    group2))  %>%
  #mutate(collection_year = collection_monthyear %>%
  #  paste0('-01') %>%
  # ymd() %>%
  #decimal_date() %>%
  #round(., digits = 1))  %>%
  summarise(n_reassortants = sum(n_reassortants), 
            .by = c(collection_monthyear,collection_regionname, group2))



# for each region, consider a transition matrix of 4 states (none, minor, major dominant)
# what are the transition probabilities fitted to these. data (ie. what is the probability of major folowing
# minor and vice versa)

library(markovchains)


# give equal prior transition probabilities, what information do our posterior distributions tell us? 