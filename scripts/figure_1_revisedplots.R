
hpai_cases <- read_csv('~/Downloads/overview-raw-data_202502241440.csv') %>%
  mutate(region = str_to_lower(Region),
         sub_region = str_to_lower(Subregion),
         observation_date = dmy(Observation.date),
         report_date = dmy(Report.date)) %>%
  separate_wider_delim(Serotype, delim = ' ', names = c('serotype', 'pathogenicity'),
                       too_few = 'align_end') %>%
  select(sub_region,
         serotype,
         pathogenicity,
         observation_date,
         report_date) %>%
  mutate(observation_date = coalesce(observation_date, report_date)) %>%
  drop_na(observation_date) %>%
  
  mutate(collection_regionname = case_when(grepl('europe', sub_region) ~ 'europe',
                                           grepl('africa', sub_region) ~ 'africa',
                                           grepl('asia', sub_region) ~ 'asia',
                                           grepl('(central|northern) america|caribbean', sub_region) ~ 'north and central america',
                                           grepl('south america|southern ocean', sub_region) ~ 'south america',
                                           grepl('australia|melanesia', sub_region) ~ 'australasia',
                                           .default = sub_region)) %>%
  select(-c(sub_region, report_date))
  

ref <- ne_countries() %>% as_tibble()  %>% select(name, continent)
woah_hpai <- read_csv('~/Downloads/Quantitative data 2025-02-24.csv') %>%
  select(Year,
         Semester,
         `World region`,
         Country,
         `Animal Category`,
         `Serotype/Subtype/Genotype`,
         `New outbreaks`, Susceptible, `Measuring units`, Cases) %>%
  mutate(across(!Year, ~ case_when(grepl('^-$', .x) ~ NA_character_,
                                          .default = .x))) %>%
  mutate(date_start = case_when(grepl('Jan-Jun', Semester) ~ paste0(Year, '-01-01'),
                                grepl('Jul-Dec', Semester) ~ paste0(Year, '-07-01'),
                                .default = NA_character_),
         date_end = case_when(grepl('Jan-Jun', Semester) ~ paste0(Year, '-06-30'),
                              grepl('Jul-Dec', Semester) ~ paste0(Year, '-12-31'),
                              .default = NA_character_)) %>%
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
  left_join(ref, by = join_by(Country == name)) %>%
  mutate(animal_category  = str_to_lower(`Animal Category`)) %>%
  mutate(across(any_of(c('New outbreaks', 'Susceptible', 'Cases')), ~as.numeric(.x))) %>%
  drop_na(Cases) %>%
  summarise(n_cases = sum(Cases, na.rm = TRUE), .by = c(date_start, continent, animal_category)) %>%
  rename(collection_regionname = continent) %>%
  mutate(collection_regionname = str_to_lower(collection_regionname)) %>%
  mutate(collection_regionname = case_when(grepl('europe', collection_regionname) ~ 'europe',
                                           grepl('africa', collection_regionname) ~ 'africa',
                                           grepl('asia', collection_regionname) ~ 'asia',
                                           grepl('(central|northern|north) america|caribbean', collection_regionname) ~  'central & northern america',
                                           grepl('south america|southern ocean|antarctica', collection_regionname) ~ 'south america',
                                           grepl('australia|melanesia|oceania', collection_regionname) ~ 'australasia',
                                           .default = collection_regionname))

# HPAI 
monthly_h5nx_cases <- hpai_cases %>% 
  filter(pathogenicity == 'HPAI') %>%
  filter(grepl('H5', serotype)) %>%
  mutate(observation_month = format.Date(observation_date, format = '%Y-%m')) %>%
  group_by(observation_month, collection_regionname) %>%
  summarise(n_cases = n()) %>%
  ungroup() %>%
  mutate(observation_month = ymd(paste0(observation_month, '-01'))) %>%
  filter(observation_month> ymd('2015-06-01')) %>% 
  rename(collection_datemonth = observation_month) 


daily_h5nx_cases <- hpai_cases %>% 
  filter(pathogenicity == 'HPAI') %>%
  filter(grepl('H5', serotype)) %>%
  #mutate(observation_month = format.Date(observation_date, format = '%Y-%m')) %>%
  group_by(observation_date, collection_regionname) %>%
  summarise(n_cases = n()) %>%
  ungroup() %>%
  filter(observation_date> ymd('2015-06-01'))  

  #group_by(region) %>%
  #mutate(n_averagecases = rollapply(n_cases, 28, mean, na.rm = TRUE, fill = NA))  



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


temp <- monthly_h5nx_cases %>% 
  left_join(sequences_month) %>%
  filter(!collection_regionname %in% c('antarctica', 'oceania', 'australasia')) %>%
  drop_na(collection_regionname,collection_datemonth) %>%
  mutate(n_sequences = case_when(is.na(n_sequences) ~ 0, .default = n_sequences)) %>%
 
  
plt_1b <-  ggplot() + 
 # geom_line(aes(x = collection_datemonth, y = n_cases, colour = collection_regionname)) + 
 
  
 # geom_bar(data = monthly_h5nx_cases, 
           #alpha = 0.7,
           #aes(x = collection_datemonth, y = n_cases*3000, colour  = collection_regionname, fill = collection_regionname), stat = 'identity') + 
  geom_line(data = woah_hpai, 
            aes(x = ymd(date_start), y = n_cases,  colour  = collection_regionname)) + 
  scale_x_date(limits = as_date(c('2018-05-01', '2024-05-01')), breaks = '6 month', date_labels = "%Y %b", 'Date') + 
  scale_y_continuous( expand = c(0,0), sec.axis = sec_axis(transform = ~ ./3000), 'WOAH Cases (N)') + 
  
  scale_fill_manual(values = region_colours) +
  scale_colour_manual(values = region_colours) +
  theme_classic() + 
  theme(legend.position = 'none')


plt_1c <-   ggplot() + 
    # geom_line(aes(x = collection_datemonth, y = n_cases, colour = collection_regionname)) + 
    
    geom_bar(data = sequences_month, 
             alpha = 0.7,
             aes(x = collection_datemonth, y = n_sequences, colour  = collection_regionname, fill = collection_regionname), stat = 'identity') + 
    scale_y_continuous( expand = c(0,0), 'GISAID Whole Genomes (n)') + 
  
    scale_x_date(limits = as_date(c('2018-05-01', '2024-05-01')), breaks = '6 month', date_labels = "%Y %b", 'Date') + 
    scale_fill_manual(values = region_colours) +
    scale_colour_manual(values = region_colours) +
    theme_classic() +
  theme(legend.position = 'none')


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
  group_by(collection_datemonth,  host_simplifiedhost,collection_regionname) %>%
  summarise(n_sequences = n()) %>%
  ungroup() %>%
  mutate(collection_datemonth = ymd(paste0(collection_datemonth, '-01'))) %>%
  drop_na(collection_regionname,collection_datemonth)


plt_1d <-   ggplot() + 
  # geom_line(aes(x = collection_datemonth, y = n_cases, colour = collection_regionname)) + 
  
  geom_bar(data = sequences_host_month, 
           alpha = 0.7,
           aes(x = collection_datemonth, y = n_sequences, colour  = host_simplifiedhost, fill = host_simplifiedhost), stat = 'identity', position = "fill") + 
  scale_y_continuous( expand = c(0,0), 'GISAID Whole Genomes (n)') + 
  
  scale_x_date(limits = as_date(c('2018-05-01', '2024-05-01')), breaks = '6 month', date_labels = "%Y %b", 'Date') + 
  scale_fill_manual(values = host_colours) +
  scale_colour_manual(values = host_colours) +
  theme_classic() +
  theme(legend.position = 'none')

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
  group_by(collection_datemonth,  host_simplifiedhost,collection_regionname) %>%
  summarise(n_sequences = n()) %>%
  ungroup() %>%
  mutate(collection_datemonth = ymd(paste0(collection_datemonth, '-01'))) %>%
  drop_na(collection_regionname,collection_datemonth)


plt_1d <-   ggplot() + 
  # geom_line(aes(x = collection_datemonth, y = n_cases, colour = collection_regionname)) + 
  
  geom_bar(data = sequences_host_month, 
           alpha = 0.7,
           aes(x = collection_datemonth, y = n_sequences, colour  = host_simplifiedhost, fill = host_simplifiedhost), stat = 'identity', position = "fill") + 
  scale_y_continuous( expand = c(0,0), 'GISAID Whole Genomes (n)') + 
  
  scale_x_date(limits = as_date(c('2018-05-01', '2024-05-01')), breaks = '6 month', date_labels = "%Y %b", 'Date') + 
  scale_fill_manual(values = host_colours) +
  scale_colour_manual(values = host_colours) +
  theme_classic() +
  theme(legend.position = 'none')



cowplot::plot_grid(plt_1b, plt_1c, plt_1d, nrow = 3, labels = c('B', 'C', 'D'), align = c('h', 'v'))


  facet_grid(rows = vars(collection_regionname), scales = "free_y")  +
    theme(legend.position = 'none')

# Sequences/Host



# Sequences/Subtype


combined_data %>% 
  filter(segment == 'ha') %>%
  select(TMRCA) %>%
  drop_na(TMRCA) %>%
  mutate(tmrca_date = date_decimal(TMRCA)) %>%
  mutate(tmrca_yearmonth = round_date(tmrca_date, 'month') %>% as_date(.)) %>%
  ggplot() + 
  geom_bar(aes(x = tmrca_yearmonth)) +
  scale_x_date('Collection Year', expand = c(0,0), limits = as_date(c('2016-01-01','2024-06-01'))) + 
  scale_y_continuous('Frequency', expand = c(0,0)) + 
  global_theme
