metadatafiles <- list.files('./data/metadata/2024Jan24',
                            full.names = T,
                            pattern = '.csv')

segnames <- c('HA', 'MP', 'N1', 'N2', 'N3', 'N6', 'N8', 'NP', 'NS', 'PA', 'PB1', 'PB2')

updated_reassortantdata <- read_csv('./data/metadata/h5_metadata_global_6280_update.csv') %>% select(c(isolate_id, clade))



Jan24_metadata <- lapply(metadatafiles, read_csv) %>%
  setNames(., segnames) %>%
  lapply(., function(x) x %>% mutate(collection.tipdate = as.character(collection.tipdate))) %>%
  bind_rows(., .id = 'segment') %>%
  left_join(.,updated_reassortantdata , by = join_by(isolate.id == isolate_id)) %>%
  mutate(clade = coalesce(clade.y, clade.x)) %>%
  select(-clade.y) 

missing_latlong <- Jan24_metadata %>% filter(!is.na(collection.subdiv1.name)) %>% filter(is.na(collection.subdiv1.long))

metadatafiles_tsv <- gsub('.csv$', '.tsv' ,metadatafiles)

Jan24_metadata_joint <- Jan24_metadata 
mapply(write_delim, 
       quote= 'needed',
       metadata_subsampled_beast, 
       metadatafiles_tsv)
