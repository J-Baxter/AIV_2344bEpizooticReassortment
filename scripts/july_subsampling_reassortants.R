library(tidyverse)
library(ape)

####################################################################################################
# Import metadata and alignments
metadatafiles <- list.files(path = './2024Jul12/reassortant_metadata',
                            pattern = 'csv',
                            full.names = T)



metadata <- lapply(metadatafiles, read_csv, col_types = cols(collection_tipdate = col_character()))

names(metadata) <- str_split_i(metadatafiles, '/', 4) %>%
  gsub('h5_|.csv', '', .)


aln_files <- list.files(path = './2024Jul12/reassortant_alignments',
                        pattern = 'fasta',
                        full.names = T)

aln <- lapply(aln_files ,
              read.dna, 
              format = 'fasta') %>%
  lapply(., as.matrix) 

names(aln) <- str_split_i(aln_files, '/', 4) %>%
  gsub('h5_|.fasta', '', .) 


####################################################################################################
# Import list of sequences to remove due to lack of temporal signal

seqstoremove <- read_csv('./2024Jul12/reassortant_iqtree/dropsequences.csv') %>%
  mutate(drop_temporal = 'drop') %>%
  mutate(isolate_id =  str_match(sequence_name, "EPI_ISL_(china_){0,1}\\d+[^.|]*")[1]) %>%
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

metadata_subsampled <- metadata %>%

  bind_rows(., .id = 'label') %>%
  mutate(clade = gsub('\\.', '', clade)) %>%
  
  # join identical sequence groupings
  left_join(identical_seqs, 
            by = join_by(label, 
                         tipnames)) %>%
  #drop seqs
  left_join(., seqstoremove, 
            by = join_by(label== label, isolate_id)) %>% 
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
  setNames(names(metadata))

####################################################################################################
# Subsample alignments
alignments_subsampled <- mapply(function(x,y) x[rownames(x) %in% y$tipnames,],
                                aln,
                                metadata_subsampled,
                                SIMPLIFY = FALSE) 
####################################################################################################
#BEAST traits
metadata_subsampled_beast <- lapply(metadata_subsampled, 
                                    function(x) x %>% 
                                      select(c(tipnames,
                                               virus_subtype, 
                                               contains('collection'), 
                                               cluster_profile,
                                               host_simplifiedhost)) %>%
                                      mutate(lat = coalesce(collection_subdiv1lat, 
                                                            collection_countrylat)) %>%
                                      mutate(long = coalesce(collection_subdiv1long, 
                                                             collection_countrylong)) %>%
                                      select(where(~n_distinct(.) > 1)) %>%
                                      select(-c(collection_countrycode,
                                                contains('date'),
                                                contains('subdiv'),
                                                collection_countrylat,
                                                collection_countrylong
                                                #collection_original, 
                                                #collection_tipdate
                                      )))


####################################################################################################
# Write alignments and metadata to file



alignmentfiles_subsampled <- paste('./2024Jul12/reassortant_beastsubsample',
                                   paste(names(aln), 'subsampled.fasta', sep = '_'),
                                   sep = '/' )

mapply(ape::write.dna, 
       alignments_subsampled, 
       alignmentfiles_subsampled ,
       format = 'fasta')


metadatafiles_subsampled_beast <-paste('./2024Jul12/reassortant_beastsubsample',
                                       paste(names(aln), 'subsampled.txt',  sep = '_'),
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
               gsub('.*beastsubsample/', '', alignmentfiles_subsampled),
               ' ',
               gsub('.fasta|.*beastsubsample/', '', alignmentfiles_subsampled),
               '_relaxLn_constant', '_1', '.xml')

write_lines(cmds,  paste('./2024Jul12/reassortant_beastsubsample', 
                         'beastgen.txt',
                         sep = '/' ))


