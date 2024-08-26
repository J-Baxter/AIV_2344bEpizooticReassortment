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
####################################################################################################
# Subsample - No dominants
summary_data <- read_csv('2024-06-05_reassortant_summary.csv') %>%
  select(cluster_profile, group) %>%
  mutate(cluster_profile = gsub('_', '', cluster_profile))

# sample the earliest example of each cluster
metadata_subsampled_nodominant_1 <- metadata_prep  %>%
  left_join(summary_data, by = join_by(cluster_profile)) %>%
  filter(cluster_profile != 'dominant') %>%
  
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
  setNames(names(metadata))



# sample the most recent example of each cluster
metadata_subsampled_nodominant_2 <- metadata_prep  %>%
  left_join(summary_data, by = join_by(cluster_profile)) %>%
  filter(cluster_profile != 'dominant') %>%
  
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
  setNames(names(metadata))



# sample according to location
metadata_subsampled_nodominant_3 <- metadata_prep %>%
  left_join(summary_data, by = join_by(cluster_profile)) %>%
  filter(cluster_profile != 'dominant') %>%
  
  group_by(label, 
           #group,
           host_simplifiedhost, 
           cluster_profile,
           best_location_code
  ) %>%
  slice_sample(n = 1) %>%
  ungroup() %>%
  
  # as list
  group_split(label, .keep = FALSE) %>%
  as.list() %>%
  setNames(names(metadata))


# sample according to host
metadata_subsampled_nodominant_4 <- metadata_prep %>%
  left_join(summary_data, by = join_by(cluster_profile)) %>%
  filter(cluster_profile != 'dominant') %>%
  
  group_by(label, 
           #group,
           host_order, 
           cluster_profile#,
           #best_location_code
  ) %>%
  slice_sample(n = 1) %>%
  ungroup() %>%
  
  # as list
  group_split(label, .keep = FALSE) %>%
  as.list() %>%
  setNames(names(metadata))

metadata_subsampled_nodominant <- lapply(1:length(metadata), function(i) bind_rows(metadata_subsampled_nodominant_1[[i]],
                                                                        metadata_subsampled_nodominant_2[[i]],
                                                                        metadata_subsampled_nodominant_3[[i]],
                                                                        metadata_subsampled_nodominant_4[[i]]) %>%
                                distinct()) 

lapply(metadata_subsampled_nodominant, nrow) %>% unlist()

####################################################################################################
# Subsample alignments
alignments_subsampled <- mapply(function(x,y) x[str_extract(rownames(x), "EPI_ISL_(china_){0,1}\\d+[^.|]*") %in% y$isolate_id,],
                                aln,
                                metadata_subsampled,
                                SIMPLIFY = FALSE) 


alignments_subsampled_nodominant <- mapply(function(x,y) x[rownames(x) %in% y$tipnames,],
                                aln,
                                metadata_subsampled_nodominant,
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
                                               cluster_number,
                                               host_simplifiedhost)) %>%
                                      mutate(lat = coalesce(collection_subdiv1lat, 
                                                                        collection_countrylat)) %>%
                                      mutate(long = coalesce(collection_subdiv1long, 
                                                                         collection_countrylong)) %>%
                                      mutate(across(c(lat, long), .fns = ~ as.numeric(.x))) %>%
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
alignmentfiles_subsampled <- paste('./2024Aug18/region_subsampled_alignments',
                                   paste(names(aln), 'subsampled.fasta', sep = '_'),
                                   sep = '/' )

mapply(ape::write.dna, 
       alignments_subsampled, 
       alignmentfiles_subsampled ,
       format = 'fasta')


metadatafiles_subsampled_beast <-paste('./2024Jul12/region_beastsubsample',
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

treeprior <- lapply(metadata_subsampled, function(x) x %>% 
                      filter(!is.na(collection_date)) %>%
                      mutate(date_dec = decimal_date(ymd(collection_date))) %>%
                      summarise(prior = (as.numeric(max(date_dec)) - as.numeric(min(date_dec))))) %>% unlist() %>% ceiling()

cmds_1 <- mapply(function(x,y) paste0("./beastgen -date_order -1 -date_prefix \\| -date_precision -D ",
                                      "'skygrid_PopSize=",
                                      min(40, x*4),
                                      ",skygrid_numGridPoints=",
                                      str_trim(format(min(40, x*4)-1, nsmall = 1)),
                                      ",skygrid_cutOff=",
                                      str_trim(format(x, nsmall = 1)),
                                      ',fileName=',
                                      gsub('.fasta$|./2024Aug18/region_subsampled_alignments/', '', y), 
                                      '_relaxLn_Skygrid', x, '-', min(40, x*4), '_1',
                                      "' flu_skygridtemplate ",
                                      gsub('./2024Aug18/region_subsampled_alignments/', '', y),
                                      ' ',
                                      gsub('.fasta|./2024Aug18/region_subsampled_alignments/', '', y),
                                      '_relaxLn_Skygrid', x, '-', min(40, x*4), '_1', '.xml'),
                 treeprior,
                 alignmentfiles_subsampled,
                 SIMPLIFY = F)
  


 
cmds_2 <- mapply(function(x,y) paste0("./beastgen -date_order -1 -date_prefix \\| -date_precision -D ",
                                      "'skygrid_PopSize=",
                                      min(40, x*4),
                                      ",skygrid_numGridPoints=",
                                      str_trim(format(min(40, x*4)-1, nsmall = 1)),
                                      ",skygrid_cutOff=",
                                      str_trim(format(x, nsmall = 1)),
                                      ',fileName=',
                                      gsub('.fasta$|./2024Aug18/region_subsampled_alignments/', '', y), 
                                      '_relaxLn_Skygrid', x, '-', min(40, x*4), '_2',
                                      "' flu_skygridtemplate ",
                                      gsub('./2024Aug18/region_subsampled_alignments/', '', y),
                                      ' ',
                                      gsub('.fasta|./2024Aug18/region_subsampled_alignments/', '', y),
                                      '_relaxLn_Skygrid', x, '-', min(40, x*4), '_2', '.xml'),
                 treeprior,
                 alignmentfiles_subsampled,
                 SIMPLIFY = F)


write_lines(c(cmds_1, cmds_2) ,  paste('./2024Aug18/region_subsampled_alignments', 
                         'beastgen.sh',
                         sep = '/' ))
