treefiles <- list.files(path = './data/alignments/subsampled_alignments/2024Jan10_strict/iqtrees',
                        recursive = FALSE,
                        include.dirs = FALSE, 
                        full.names = TRUE) %>%
  .[grep('treefile$', .)]

trees <- lapply(treefiles, read.newick) 


test_treedata <- tidytree::as_tibble(trees[[1]]) %>%
  left_join(., HA_strict, by = join_by(label == tipnames)) %>%
  as.treedata()



metadata_subsampled_strict_df <- metadata_subsampled_strict %>%
  bind_rows(., .id = 'virus.segment') %>%
  filter(virus.segment == 'HA') %>%
  ggplot() +
  geom_bar(aes(x = `collection.country.code`)) + 
  facet_wrap(host.~virus.segment)






####################################################################################################
# Subsample - strict

metadata_subsampled_strict_test <- metadata %>%
  bind_rows(., .id = 'virus.segment') %>%
  
  # join identical sequence groupings
  left_join(identical_seqs, 
            by = join_by(virus.segment, 
                         tipnames)) %>%
  # create groupings
  mutate(across(contains('collection'), .fns = ~ gsub('^NA$', NA, .x))) %>% 
  mutate(best_location_code = coalesce(collection.subdiv1.code,
                                       collection.country.code)) %>%
  mutate(group = as.factor(group)) %>%
  
  # Sample one virus in each segment that has unique combination of location, host family, 
  # reasssortment, and sequence identity (viruses in the same group are identical)
  
  group_by(virus.segment, 
           group,
           host.family, 
           cluster.genome) %>%
  mutate(todrop = case_when(collection.country.code %in% c('RU', 'JP', 'CN') && is.na(collection.subdiv1.code) && n()>1 ~ 'drop',
                            .default = 'keep')) %>%
  filter(todrop == 'keep') %>%
  ungroup() %>%
  group_by(virus.segment, 
           group,
           host.family, 
           cluster.genome,
           best_location_code) %>%
  slice_sample(n = 1) %>%
  ungroup() %>%
  
  # as list
  group_split(virus.segment, .keep = FALSE) %>%
  as.list() %>%
  setNames(segnames)


####################################################################################################
metadata_subsampled_relaxed_test <- metadata %>%
  bind_rows(., .id = 'virus.segment')%>%
  
  # join identical sequence groupings
  left_join(identical_seqs, 
            by = join_by(virus.segment, 
                         tipnames)) %>%
  # create groupings
  mutate(across(contains('collection'), .fns = ~ gsub('^NA$', NA, .x))) %>% 
  mutate(best_location_code = coalesce(collection.subdiv1.code,
                                       collection.country.code)) %>%
  mutate(group = as.factor(group)) %>%
  
  # Sample one virus in each segment that has unique combination of location, host family, 
  # reasssortment, and sequence identity (viruses in the same group are identical)
  
  group_by(group,
           host.family, 
           best_location_code, #Leave in if relaxed, leave out if strict
           cluster.genome) %>%
  slice_sample(n = 1) %>%
  ungroup() %>%
  
  # as list
  group_split(virus.segment, .keep = FALSE) %>%
  as.list() %>%
  setNames(segnames)
