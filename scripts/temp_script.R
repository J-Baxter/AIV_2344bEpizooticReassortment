################### Import Evolutionary Rates from Reassortant log files ###################
reassortant_log_path <- list.files("./2024Nov18/plain_tree",
                                   pattern = '.log',
                                   full.names = T)



imported_logs <- lapply(reassortant_log_path, readLog) %>%
  lapply(., getHPDMedian) %>%
  lapply(., function(x) as_tibble(x, rownames = 'parameter') %>% filter(parameter %in% c("ucld.mean", "age.root."))) %>%
  setNames(gsub(".*plain_tree/|_subsampled.*", "", reassortant_log_path)) %>%
  bind_rows(.id = 'key') %>%
  mutate(key = gsub("n[:0-9:]", "nx", key),
         parameter = gsub('age\\.root\\.', 'TMRCA', parameter) %>%
           gsub('ucld\\.mean', 'evoRate', .)) %>%
  rename(`95lower` = lower,
         `95upper` = upper) %>%
  pivot_wider(names_from = parameter,
              values_from = c('95lower', 'med', '95upper'),
              names_glue =  "{parameter}_{.value}") %>%
  rename_with(.fn = ~gsub('_med', '', .x) , .cols = contains('med')) %>%
  
  mutate(cluster_profile = gsub('[[:alpha:]][[:alpha:]][[:digit:]]{0,1}_', '', key))


# left join most recent and oldest sequence for each dominant reassortant
sampling_dates <- read_csv('2024-09-09_meta.csv') %>%
  
  # Separate1_1_1_1_1_1_1_1 into A and B clades
  mutate(cluster_profile  = case_when(grepl('asia', collection_regionname) & cluster_profile == '1_1_1_1_1_1_1_1' ~ paste0(cluster_profile, 'A'), 
                                      !grepl('asia', collection_regionname) & cluster_profile == '1_1_1_1_1_1_1_1' ~ paste0(cluster_profile, 'B'), 
                                      .default = cluster_profile)) %>%
  group_by(cluster_profile) %>%
  slice(which(collection_date == max(collection_date, na.rm = T) | collection_date == min(collection_date, na.rm = T)), 
        .preserve = T) %>%
  group_by(collection_date, .add = TRUE) %>%
  slice_sample(n = 1) %>%
  ungroup(collection_date) %>%
  mutate(var = case_when(collection_date== max(collection_date, na.rm = T) ~ 'most_recent',
                         collection_date==min(collection_date, na.rm = T) ~'most_distant'),
         cluster_profile = gsub('_', '', cluster_profile)) %>%
  select(cluster_profile, collection_date, var) %>%
  ungroup() %>%
  filter(cluster_profile %in% c("32313212",
                                "21111111",
                                "11111111B",
                                "11111111A",
                                "21211111",
                                "11211111",
                                "16211111",
                                "11414114",
                                "26116111",
                                "21611411",
                                "71521312",
                                "43112113",
                                "51112113",
                                "54912111"))%>%
  pivot_wider(names_from = var,
              values_from = collection_date) %>%

  mutate(across(starts_with('most'), .fns = ~ decimal_date(.x)))

old_data <- read_csv('./2024Aug18/treedata_extractions/2024-09-20_combined_data.csv') %>%
  dplyr::select(key, starts_with('evoRate'), TMRCA) %>%
  filter(grepl("32313212$|21111111$|11111111$|11111111A$|21211111$|11211111$|16211111$|11414114$|26116111$|21611411$|71521312$|43112113$|51112113$|54912111$", key)) %>%
  mutate(key = case_when(grepl('11111111$', key) ~ gsub('11111111$', '11111111B', key), 
                         .default = key))

test <- expand_grid(segment = c('ha', 'mp', 'np', 'ns', 'nx', 'pa', 'pb1', 'pb2'),
                    cluster_profile = c("32313212",
                                        "21111111",
                                        "11111111B",
                                        "11111111A",
                                        "21211111",
                                        "11211111",
                                        "16211111",
                                        "11414114",
                                        "26116111",
                                        "21611411",
                                        "71521312",
                                        "43112113",
                                        "51112113",
                                        "54912111")) %>%
  left_join(sampling_dates) %>%
  left_join(imported_logs %>%
              mutate(segment = gsub('_.*', '', key)) %>%
              dplyr::select(-key)) %>%
  unite(key, segment, cluster_profile, remove = F) %>%
  rows_patch(old_data, by = 'key', unmatched = 'ignore')
  

test %>%
  ggplot() +
  
  # Geom objects
  geom_linerange(aes(x = cluster_profile,
                     ymin = most_distant, 
                     ymax = most_recent,
                     colour = segment), 
                 linewidth = 1.5, 
                 position = position_dodge(width = 0.75)) +
  
  geom_linerange(aes(x = cluster_profile,
                     ymin = TMRCA, 
                     ymax = most_recent,
                     colour = segment), 
                 linewidth = 2, 
                 alpha = 0.5,
                 position = position_dodge(width = 0.75)) +
  # Scales
  scale_x_discrete('Reassortant Profile')+
  scale_y_continuous()+
  scale_colour_brewer(palette = 'Dark2') +
  coord_flip(ylim = c(2015, 2025)) + 
 
  
  # Graphical
  #facet_grid(rows = vars(segment)) + 
  theme_minimal() 