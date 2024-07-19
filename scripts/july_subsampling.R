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
           #cluster_profile
           ) %>%
  mutate(todrop = case_when(collection_countrycode %in% c('RU', 'JP', 'CN') && is.na(collection_subdiv1code) && n()>1 ~ 'drop',
                            collection_countrycode == 'IRN' ~ 'drop',
                            .default = 'keep')) %>%
  filter(todrop == 'keep') %>%
  ungroup() %>%
  group_by(label, 
           group,
           host_simplifiedhost, 
           #cluster_profile,
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
                                alignments,
                                metadata_subsampled,
                                SIMPLIFY = FALSE) 
####################################################################################################
#BEAST traits
metadata_subsampled_beast <- lapply(metadata_subsampled, 
                                    function(x) x %>% 
                                      select(c(tipnames,
                                               virus_subtype, 
                                               contains('collection'), 
                                               cluster_number,
                                               host_order,
                                               host_simplifiedhost)) %>%
                                      mutate(tiplocation_lat = coalesce(collection_subdiv1lat, 
                                                                        collection_countrylat)) %>%
                                      mutate(tiplocation_long = coalesce(collection_subdiv1long, 
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



alignmentfiles_subsampled <- paste('.',
                                   ddmonthyy, 
                                   'alignments',
                                   paste(segnames[!segnames %in% problems], '2344b_subsampled.fasta', sep = '_'),
                                   sep = '/' )

mapply(ape::write.dna, 
       alignments_subsampled, 
       alignmentfiles_subsampled ,
       format = 'fasta')


metadatafiles_subsampled_beast <-paste('.',
                                       ddmonthyy, 
                                       'metadata',
                                       paste(segnames[!segnames %in% problems], '2344b_subsampled.txt',  sep = '_'),
                                       sep = '/' )

mapply(write_delim, 
       delim = '\t',
       quote= 'needed',
       metadata_subsampled_beast, 
       metadatafiles_subsampled_beast)  


metadatafiles_subsampled <-paste('.',
                                 ddmonthyy, 
                                 'metadata',
                                 paste(segnames[!segnames %in% problems], '2344b_subsampled.csv',  sep = '_'),
                                 sep = '/' )

mapply(write_csv, 
       quote= 'needed',
       metadata_subsampled, 
       metadatafiles_subsampled)  

####################################################################################################
# Write plain BEAUTI XML and metadata to file
# SRD06, relaxed lognormal and skygrid coalescent
# lognormal prior 

treeprior <- lapply(metadata_subsampled, function(x) x %>% filter(!is.na(collection_datedecimal)) %>% summarise(prior = (as.numeric(max(collection_datedecimal)) - as.numeric(min(collection_datedecimal))))) %>% unlist() %>% ceiling()

cmds <- paste0("./beastgen -date_order -1 -date_prefix . -date_precision  -D ",
               "'skygrid_PopSize=",
               treeprior*4,
               ",skygrid_numGridPoints=",
               format(treeprior*4-1, nsmall = 1),
               ",skygrid_cutOff=",
               format(treeprior, nsmall = 1),
               ',fileName=',
               gsub('.prior', '', names(treeprior)), 
               '_relaxLn_Skygrid', treeprior, '-', treeprior*4, '_1',
               "' skygridtemplate ",
               alignmentfiles_subsampled,
               ' ',
               gsub('.prior', '', names(treeprior)), 
               '_relaxLn_Skygrid', treeprior, '-', treeprior*4, '_1', '.xml')

write_lines(cmds,  paste('.',
                         ddmonthyy, 
                         'beastgen.txt',
                         sep = '/' ))
