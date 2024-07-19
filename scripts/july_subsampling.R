library(ape)
library(tidyverse)
####################################################################################################
# Import metadata and alignments
metadatafiles <- list.files(path = './2024Jul12/region_metadata',
                            pattern = 'csv',
                            full.names = T)



metadata <- lapply(metadatafiles, read_csv, col_types = cols(collection_tipdate = col_character()))
names(metadata) <- str_split_i(metadatafiles, '/', 4) %>%
  gsub('.csv', '', .)


aln_files <- list.files(path = './2024Jul12/region_alignments',
                        pattern = 'fasta',
                        full.names = T)

aln <- lapply(aln_files ,
              read.dna, 
              format = 'fasta') %>%
  lapply(., as.matrix) 

names(aln) <- str_split_i(aln_files, '/', 4) %>%
  gsub('.fasta', '', .) %>%
  gsub('h5_', '', .)


####################################################################################################
# Import list of sequences to remove due to lack of temporal signal
seqstoremove <- read_csv('./2024Jul12/region_iqtree/dropsequences.csv') %>%
  mutate(drop_temporal = 'drop') %>%
  mutate(isolate_id =  str_match(sequence_name, "EPI_ISL_\\d+[^.|]*")) %>%
  unite(label, segment, region) %>%
  select(-sequence_name) %>% 
  distinct()


####################################################################################################
# Identify duplicate sequences from alignment

identical_seqs <- lapply(aln,
                         FindIdenticalSeqs) %>%
  bind_rows(., .id = 'label')


####################################################################################################
# Subsample - strict

metadata_subsampled_strict_test <- metadata %>%
  bind_rows(., .id = 'label') %>%
  
  # join identical sequence groupings
  left_join(identical_seqs, 
            by = join_by(label, 
                         tipnames)) %>%
  #drop seqs
  left_join(., seqstoremove, 
            by = join_by(label == label, isolate_id)) %>% 
  filter(is.na(drop_temporal)) %>%
  select(-drop_temporal) %>%
  
  # only 2344b
  filter(clade == '2344b' | is.na(clade)) %>%
  
  # create groupings
  mutate(across(contains('collection'), .fns = ~ gsub('^NA$', NA, .x))) %>% 
  mutate(best_location_code = coalesce(collection_subdiv1code,
                                       collection_countrycode)) %>%
  mutate(group = as.factor(group)) %>%
  
  # Sample one virus in each segment that has unique combination of location, host family, 
  # reasssortment, and sequence identity (viruses in the same group are identical)
  
  group_by(label, 
           group,
           host_simplifiedhost, 
           cluster_profile) %>%
  mutate(todrop = case_when(collection_countrycode %in% c('RU', 'JP', 'CN') && is.na(collection_subdiv1code) && n()>1 ~ 'drop',
                            collection_countrycode == 'IRN' ~ 'drop',
                            .default = 'keep')) %>%
  filter(todrop == 'keep') %>%
  ungroup() %>%
  group_by(label, 
           group,
           host_simplifiedhost, 
           cluster_profile,
           best_location_code) %>%
  slice_sample(n = 1) %>%
  ungroup() %>%
  
  # as list
  group_split(label, .keep = FALSE) %>%
  as.list() %>%
  setNames(names(metadata))
