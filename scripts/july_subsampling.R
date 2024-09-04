library(ape)
library(tidyverse)
source('./scripts/FindIdenticalSeqs.R')

####################################################################################################
# Import metadata and alignments
metadata <- read_csv('2024-08-19_meta.csv')

aln_files <- list.files(path = './2024Aug18/region_alignments',
                        pattern = 'fasta',
                        full.names = T)

aln <- lapply(aln_files ,
              read.dna, 
              format = 'fasta',
              as.matrix = TRUE)

names(aln) <- str_split_i(aln_files, '/', 4) %>%
  gsub('h5_|.fasta', '', .) 


metadata_per_alignment <- lapply(aln, function(x) metadata %>% 
                                   filter(isolate_id %in% str_extract(rownames(x), "EPI_ISL_(china_){0,1}\\d+[^.|]*")))

names(metadata_per_alignment) <- str_split_i(aln_files, '/', 4) %>%
  gsub('h5_|.fasta', '', .) 

all(unlist(lapply(aln, nrow)) == unlist(lapply(metadata_per_alignment, nrow)))


####################################################################################################
# Import list of sequences to remove due to lack of temporal signal
seqstoremove <- read_csv('./2024Jul12/region_iqtree/dropsequences.csv') %>%
  mutate(drop_temporal = 'drop') %>%
  mutate(isolate_id =  str_extract(sequence_name, "EPI_ISL_(china_){0,1}\\d+[^.|]*")) %>%
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
metadata_prep <- metadata_per_alignment %>%
  bind_rows(., .id = 'label') %>%
  
  
  # join identical sequence groupings
  left_join(identical_seqs, 
            by = join_by(label, 
                         tipnames)) %>%
  #drop seqs
  left_join(., seqstoremove, 
            by = join_by(label == label, isolate_id)) %>% 
  filter(is.na(drop_temporal)) %>%
  filter(host_simplifiedhost != 'unknown') %>%
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
           #host_simplifiedhost, 
          cluster_profile
  ) %>%
  mutate(todrop = case_when(collection_countrycode %in% c('RU', 'JP', 'CN') && is.na(collection_subdiv1code) && n()>1 ~ 'drop',
                            collection_countrycode == 'IRN' ~ 'drop',
                            .default = 'keep')) %>%
  filter(todrop == 'keep') %>%
  ungroup()


# sample the earliest example of each cluster
metadata_subsampled_1 <- metadata_prep  %>%
  group_by(label, 
           #group,
           #host_simplifiedhost, 
           cluster_profile#,
           # best_location_code
  ) %>%
  slice_min(order_by = collection_date, n = 1) %>%
  #slice_sample(n = 1) %>%
  ungroup() %>%
  
  # as list
  group_split(label, .keep = FALSE) %>%
  as.list() %>%
  setNames(names(metadata_per_alignment))



# sample the most recent example of each cluster
metadata_subsampled_2 <- metadata_prep  %>%
  group_by(label, 
           #group,
           #host_simplifiedhost, 
           cluster_profile#,
           # best_location_code
  ) %>%
  slice_max(order_by = collection_date, n = 1) %>%
  #slice_sample(n = 1) %>%
  ungroup() %>%
  
  # as list
  group_split(label, .keep = FALSE) %>%
  as.list() %>%
  setNames(names(metadata_per_alignment))



# sample according to location
metadata_subsampled_3 <- metadata_prep %>%
  group_by(label, 
           #group,
           #host_order, 
           cluster_profile,
           best_location_code
  ) %>%
  slice_sample(n = 1) %>%
  ungroup() %>%
  
  # as list
  group_split(label, .keep = FALSE) %>%
  as.list() %>%
  setNames(names(metadata_per_alignment))



# sample according to location
metadata_subsampled_4 <- metadata_prep %>%
  group_by(label, 
           #group,
           host_simplifiedhost, 
           cluster_profile#,
           #best_location_code
  ) %>%
  slice_sample(n = 1) %>%
  ungroup() %>%
  
  # as list
  group_split(label, .keep = FALSE) %>%
  as.list() %>%
  setNames(names(metadata_per_alignment))

metadata_subsampled <- lapply(1:length(metadata_per_alignment), function(i) bind_rows(metadata_subsampled_1[[i]],
                                                                        metadata_subsampled_2[[i]],
                                                                        metadata_subsampled_3[[i]],
                                                                        metadata_subsampled_4[[i]]) %>%
                                distinct()) %>%
  setNames(names(metadata_per_alignment))

lapply(metadata_subsampled, nrow) %>% unlist()
