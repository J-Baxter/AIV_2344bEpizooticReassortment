library(ape)
library(tidyverse)
source('./scripts/FindIdenticalSeqs.R')

####################################################################################################
# Import metadata and alignments
metadata <- read_csv('2024-08-19_meta.csv') %>%
  filter(cluster_profile == '1_1_1_1_1_1_1_1') %>%
  mutate(cluster_profile  = case_when(grepl('asia', collection_regionname) ~ paste0(cluster_profile, 'A'), 
                                      .default = cluster_profile)) %>%
  filter(cluster_profile == '1_1_1_1_1_1_1_1A')

aln_files <- list.files(path = './2024Aug18/reassortant_alignments',
                        pattern = 'fasta',
                        full.names = T) %>%
  .[grepl('11111111',.)]

aln <- lapply(aln_files ,
              read.dna, 
              format = 'fasta',
              as.matrix = TRUE)

names(aln) <- str_split_i(aln_files, '/', 4) %>%
  gsub('h5_|.fasta', '', .) 

####################################################################################################
# subset alignment
aln_11111111a <- lapply(aln_files ,
                        read.dna, 
                        format = 'fasta',
                        as.matrix = TRUE) %>%
  lapply(., function(x) x[str_extract(rownames(x), "EPI_ISL_(china_){0,1}\\d+[^.|]*") %in% metadata$isolate_id,])

names(aln_11111111a) <- str_split_i(aln_files, '/', 4) %>%
  gsub('h5_|.fasta', '', .) %>%
  paste0(., 'A')

metadata_per_alignment <- lapply(aln, function(x) metadata %>% 
                                   filter(isolate_id %in% str_extract(rownames(x), "EPI_ISL_(china_){0,1}\\d+[^.|]*")))

names(metadata_per_alignment) <- str_split_i(aln_files, '/', 4) %>%
  gsub('h5_|.fasta', '', .) %>%
  paste0(., 'A')

all(unlist(lapply(aln_11111111a, nrow)) == unlist(lapply(metadata_per_alignment, nrow)))

####################################################################################################
# Import list of sequences to remove due to lack of temporal signal

seqstoremove <- read_csv('./2024Jul12/reassortant_iqtree/dropsequences.csv') %>%
  mutate(drop_temporal = 'drop') %>%
  mutate(isolate_id =  str_extract(sequence_name, "EPI_ISL_(china_){0,1}\\d+[^.|]*")) %>%
  mutate(reassortant = gsub('_', '', reassortant)) %>%
  unite(label, segment, reassortant) %>%
  select(-sequence_name) %>% 
  distinct()


####################################################################################################
# Identify duplicate sequences from alignment
identical_seqs <- lapply(aln,
                         FindIdenticalSeqs) %>%
  bind_rows(., .id = 'label') %>%
  rename(., sequence_group = group)


####################################################################################################
# Subsample - strict
metadata_subsampled <- metadata_per_alignment %>%
  
  bind_rows(., .id = 'label') %>%
  mutate(clade = gsub('\\.', '', clade)) %>%
  
  # join identical sequence groupings
  left_join(identical_seqs, 
            by = join_by(label, 
                         tipnames)) %>%
  #drop seqs
  left_join(., seqstoremove, 
            by = join_by(label== label, isolate_id)) %>% 
  filter(host_simplifiedhost != 'unknown') %>%
  filter(is.na(drop_temporal)) %>%
  select(-drop_temporal) %>%
  
  # only 2344b
  filter(clade == '2344b' | is.na(clade)) %>%
  
  # create groupings
  mutate(across(contains('collection'), .fns = ~ gsub('^NA$', NA, .x))) %>% 
  mutate(best_location_code = coalesce(collection_subdiv1code,
                                       collection_countrycode)) %>%
  mutate(sequence_group = as.factor(sequence_group)) %>%
  
  # Sample one virus in each segment that has unique combination of location, host family, 
  # reasssortment, and sequence identity (viruses in the same group are identical)
  
  group_by(label, 
           sequence_group,
           host_simplifiedhost
  ) %>%
  #mutate(todrop = case_when(collection_countrycode %in% c('RU', 'JP', 'CN') && is.na(collection_subdiv1code) && n()>1 ~ 'drop',
  #collection_countrycode == 'IRN' ~ 'drop',
  # .default = 'keep')) %>%
  #filter(todrop == 'keep') %>%
  ungroup() %>%
  group_by(label, 
           sequence_group,
           host_simplifiedhost, 
           best_location_code) %>%
  slice_sample(n = 1) %>%
  ungroup() %>%
  
  # as list
  group_split(label, .keep = FALSE) %>%
  as.list() %>%
  setNames(names(metadata_per_alignment))

####################################################################################################
# Subsample alignments
alignments_subsampled <- mapply(function(x,y) x[str_extract(rownames(x), "EPI_ISL_(china_){0,1}\\d+[^.|]*") %in% y$isolate_id,],
                                aln_11111111a,
                                metadata_subsampled,
                                SIMPLIFY = FALSE) 
####################################################################################################
#BEAST traits
metadata_subsampled_beast <- lapply(metadata_subsampled, 
                                    function(x) x %>% 
                                      select(c(tipnames,
                                               virus_subtype, 
                                               collection_regionname, 
                                               ends_with('long'),
                                               ends_with('lat'),
                                               cluster_profile,
                                               host_simplifiedhost)) %>%
                                      mutate(lat = coalesce(collection_subdiv1lat, 
                                                            collection_countrylat)) %>%
                                      mutate(long = coalesce(collection_subdiv1long, 
                                                             collection_countrylong)) %>%
                                      select(where(~n_distinct(.) > 1)) %>%
                                      select(-c(contains('date'),
                                                contains('subdiv'),
                                                collection_countrylat,
                                                collection_countrylong,
                                                #collection_original, 
                                                #collection_tipdate
                                      )))


####################################################################################################
# Write alignments and metadata to file



alignmentfiles_subsampled <- paste('./2024Aug18/11111111A',
                                   paste(names(aln_11111111a), 'subsampled.fasta', sep = '_'),
                                   sep = '/' )

mapply(ape::write.dna, 
       alignments_subsampled, 
       alignmentfiles_subsampled ,
       format = 'fasta')


metadatafiles_subsampled_beast <-paste('./2024Aug18/11111111A',
                                       paste(names(aln_11111111a), 'subsampled.txt',  sep = '_'),
                                       sep = '/' )

mapply(write_delim, 
       delim = '\t',
       quote= 'needed',
       metadata_subsampled_beast, 
       metadatafiles_subsampled_beast)  


####################################################################################################
# Write plain BEAUTI XML and metadata to file
# SRD06, relaxed lognormal and skygrid coalescent
# lognormal prior 


cmds <- paste0("./beastgen -date_order -1 -date_prefix \\| -date_precision -D '",
               'fileName=',
               gsub('.fasta$|.*beastsubsample/', '', alignmentfiles_subsampled), 
               '_relaxLn_constant', '_1',
               "' flu_constanttemplate ",
               gsub('.*beastsubsample/|./2024Aug18/reassortant_subsampled_alignments/', '', alignmentfiles_subsampled),
               ' ',
               gsub('.fasta|.*beastsubsample/|./2024Aug18/reassortant_subsampled_alignments/', '', alignmentfiles_subsampled),
               '_relaxLn_constant', '.xml')

write_lines(cmds,  paste('./2024Aug18/11111111A', 
                         'beastgen.sh',
                         sep = '/' ))


