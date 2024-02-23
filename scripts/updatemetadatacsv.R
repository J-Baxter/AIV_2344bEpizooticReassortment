metadatafiles <- list.files('./data/metadata/2024Jan24',
                            full.names = T,
                            pattern = '.tsv')

segnames <- c('HA', 'MP', 'N1', 'N2', 'N3', 'N6', 'N8', 'NP', 'NS', 'PA', 'PB1', 'PB2')

updated_reassortantdata <- read_csv('./data/metadata/h5_metadata_global_6280_update.csv') %>% select(c(isolate_id, clade))



Jan24_metadata <- lapply(metadatafiles, read_tsv) %>%
  setNames(., segnames) %>%
  lapply(., function(x) x %>% mutate(collection.tipdate = as.character(collection.tipdate))) %>%
  bind_rows(., .id = 'segment') %>%
  left_join(.,updated_reassortantdata , by = join_by(isolate.id == isolate_id)) %>%
  mutate(clade = coalesce(clade.y, clade.x)) %>%
  select(-clade.y) 

missing_latlong <- Jan24_metadata %>% filter(!is.na(collection.subdiv1.name)) %>% filter(is.na(collection.subdiv1.long))
missing_latlong %>% select(c(collection.country.name, collection.subdiv1.name)) %>% distinct() %>% arrange(collection.country.name) %>%
  mutate(collection.subdiv1.lat = NA) %>%
  mutate(collection.subdiv1.long = NA) %>%
  write.csv(file = 'missinglatlong.csv')

corrected_latong <- read_csv('./missinglatlong.csv') %>% select(-1)

reformatted <- Jan24_metadata %>%
  rows_update(corrected_latong, by = c('collection.country.name', 'collection.subdiv1.name')) %>%
  
  # segregate Russia europe/asian
  mutate(collection.region.name = case_when(collection.country.name == 'russia' & collection.subdiv1.long > 60 &  collection.subdiv1.long < 105 ~ 'central asia',
                                            collection.country.name == 'russia' & collection.subdiv1.long >= 105 ~ 'eastern asia',
                                            .default = collection.region.name)) %>%
  # replace wild NA with unknown
  mutate(host.class = case_when(is.na(host.order) ~ 'unknown',
                                grepl('environment', host.order) ~ 'unknown',
                                .default = host.class)) %>%
  # update best lat/lon
  select(-starts_with('best')) %>%
  mutate(tiplocation_lon = coalesce(collection.subdiv2.long , collection.subdiv1.long, collection.country.long)) %>%
  mutate(tiplocation_lat = coalesce(collection.subdiv2.lat , collection.subdiv1.lat, collection.country.lat)) %>%
  group_split(segment)

  
  # Update clade and cluster profile (per segment)
  
  

metadatafiles_tsv <- gsub('.csv$', '.tsv' ,metadatafiles)

Jan24_metadata_joint <- Jan24_metadata 
mapply(write_delim, 
       quote= 'needed',
       metadata_subsampled_beast, 
       metadatafiles_tsv)
