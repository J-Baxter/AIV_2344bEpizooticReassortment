combined_data %>% select(cluster_profile,
                         segment,
                         count_cross_species,
                         collection_regionname,
                         count_to_mammal) %>%
  mutate(segment = factor(segment, levels = c('pb2', 'pb1', 'pa', 'ha', 'np', 'nx', 'mp', 'ns'))) %>%
  mutate(collection_regionname = case_when(grepl('europe', collection_regionname) ~ 'europe',
                                           grepl('africa', collection_regionname) ~ 'africa',
                                           grepl('asia', collection_regionname) ~ 'asia',
                                           grepl('(central|northern) america', collection_regionname) ~ 'central & northern america',
                                           grepl('south america|southern ocean', collection_regionname) ~ 'south america',
                                           grepl('australia', collection_regionname) ~ 'australasia',
                                           .default = collection_regionname)) %>%
  summarise(count_cross_species = sum(count_cross_species, na.rm = T), 
            count_to_mammal = sum(count_to_mammal, na.rm = T), 
            .by = c(collection_regionname, segment)) %>%
  
  #replace_na(list(count_to_mammal = 0, count_cross_species = 0)) %>% 
  #filter(segment == 'ha') %>%
  pivot_longer(cols = starts_with('count'), values_to = 'frequency', names_to = 'jumps') %>%
  ggplot(aes(x = collection_regionname, y = frequency, fill = jumps)) + 
  geom_bar(stat = 'identity', position = 'dodge') + 
  facet_wrap(~ segment)
  





# Stratified by region (for supplementary)
counts %>%
  filter(collection_regionname %in% c('europe', 'asia','central & northern america','africa')) %>%
  pivot_longer(starts_with('n_'), names_to = 'var', values_to = 'value') %>%
  mutate(collection_monthyear = ym(collection_monthyear)) %>%
  filter(collection_monthyear > '2017-01-01') %>%
  ggplot(aes(x=collection_monthyear, y= value, colour = var)) + 
  geom_bar(aes(x = ym(collection_monthyear), y = n_reassortants*30), stat = 'identity', data = count_data %>% filter(collection_monthyear > '2017-01-01'), inherit.aes = FALSE) + 
  geom_point(size =1) +
  geom_line(aes(group = var)) +
  facet_grid(rows = vars(collection_regionname)) + 
  scale_x_date()  +
  scale_y_continuous(sec.axis = sec_axis( trans=~./30, name="Reassortment")) +
  global_theme

# Global (for supplementary)
counts %>%
  filter(collection_regionname %in% c('europe', 'asia','central & northern america','africa')) %>%
  pivot_longer(starts_with('n_'), names_to = 'var', values_to = 'value') %>%
  summarise(value = sum(value), .by =  c(collection_monthyear, var)) %>%
  mutate(collection_monthyear = ym(collection_monthyear)) %>%
  filter(collection_monthyear > '2017-01-01') %>%
  ggplot(aes(x=collection_monthyear, y= value, colour = var)) + 
  geom_bar(aes(x = ym(collection_monthyear), y = n_reassortants*30), stat = 'identity', data = count_data %>% filter(collection_monthyear > '2017-01-01'), inherit.aes = FALSE) + 
  geom_point(size =1) +
  geom_line(aes(group = var)) +
  scale_x_date()  +
  scale_y_continuous(sec.axis = sec_axis( trans=~./30, name="Reassortment"))+
  global_theme