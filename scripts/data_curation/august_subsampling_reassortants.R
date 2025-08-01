library(ape)
library(tidyverse)
source('./scripts/FindIdenticalSeqs.R')

####################################################################################################
# Import metadata 
metadata <- read_csv('2024-09-09_meta.csv') %>%
  
  # Separate1_1_1_1_1_1_1_1 into A and B clades
  mutate(cluster_profile  = case_when(grepl('asia', collection_regionname) & cluster_profile == '1_1_1_1_1_1_1_1' ~ paste0(cluster_profile, 'A'), 
                                      !grepl('asia', collection_regionname) & cluster_profile == '1_1_1_1_1_1_1_1' ~ paste0(cluster_profile, 'B'), 
                                      .default = cluster_profile))

# Import alignments
aln_files <- list.files(path = './2024Aug18/reassortant_alignments',
                        pattern = 'fasta',
                        full.names = T) 

aln <- lapply(aln_files ,
              read.dna, 
              format = 'fasta',
              as.matrix = TRUE)

names(aln) <- str_split_i(aln_files, '/', 4) %>%
  gsub('h5_|.fasta', '', .) 


# Split 1_1_1_1_1_1_1_1 into two separate alignments
temp_11111111 <- lapply(aln[grepl( '_11111111$', names(aln))],
                        function(x) list(x[str_extract(rownames(x), "EPI_ISL_(china_){0,1}\\d+[^.|]*") %in% (metadata %>%
                                                                                                               filter(cluster_profile == '1_1_1_1_1_1_1_1A') %>%
                                                                                                               pull(isolate_id)),], 
                                         x[str_extract(rownames(x), "EPI_ISL_(china_){0,1}\\d+[^.|]*") %in% (metadata %>% 
                                                                                                               filter(cluster_profile == '1_1_1_1_1_1_1_1B') %>%
                                                                                                               pull(isolate_id)),])) %>%
  lapply(., setNames, c('A', 'B')) %>%
  list_flatten(., name_spec = "{outer}{inner}")

# join split alignments and remove old
aln %<>% 
  .[!grepl( '_11111111$', names(aln))] %>%
  c(., temp_11111111) %>%
  .[order(names(.))]

# check all sequences in metadata
aln %<>%
  lapply(., function(x) x[str_extract(rownames(x), "EPI_ISL_(china_){0,1}\\d+[^.|]*") %in% metadata$isolate_id,])

# split metadata by alignment
metadata_per_alignment <- lapply(aln, function(x) metadata %>% 
                                   filter(isolate_id %in% str_extract(rownames(x), "EPI_ISL_(china_){0,1}\\d+[^.|]*")))


# sanity check 
all(unlist(lapply(aln, nrow)) == unlist(lapply(metadata_per_alignment, nrow)))

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
groups <- lapply(aln, GroupSequences, snp_threshold = 2)%>%
  bind_rows(., .id = 'label') %>%
  mutate(isolate_id = str_extract(tipnames, "EPI_ISL_(china_){0,1}\\d+[^.|]*"))



# 40 largest countries by land mass
large_countries <- c('russia',
                     'canada',
                     'china',
                     'united states',
                     'brazil',
                     'australia',
                     'india',
                     'argentina',
                     'kazakhstan',
                     'algeria',
                     'saudi arabia',
                     'méxico',
                     'indonesia',
                     'sudan',
                     'libya',
                     'iran',
                     'mongolia',
                     'peru',
                     'chad',
                     'niger',
                     'angola',
                     'mali',
                     'south africa',
                     'colombia',
                     'Ethiopia',
                     'Bolivia' ,
                     'Mauritania',
                     'Egypt',
                     'Tanzania' ,
                     'Nigeria' ,
                     'Venezuela' ,
                     'Pakistan' ,
                     'Namibia' ,
                     'Mozambique' ,
                     'Turkey',
                     'Chile' ,
                     'Zambia' ,
                     'Myanmar')
       
                     
 ####################################################################################################
# Subsample - strict
metadata_prep <- metadata_per_alignment %>%
  bind_rows(., .id = 'label') %>%
  
  #drop seqs
  left_join(., seqstoremove, 
            by = join_by(label == label, isolate_id)) %>% 
  filter(is.na(drop_temporal)) %>%
  filter(host_simplifiedhost != 'unknown') %>%
  select(-drop_temporal) %>%
  
  # join identical sequence groupings
  left_join(groups %>% select(-tipnames), 
            by = join_by(label, 
                         isolate_id)) %>%
  
  # only 2344b
  filter(clade == '2344b' | is.na(clade)) %>%
  
  # create groupings
  mutate(across(contains('collection'), .fns = ~ gsub('^NA$', NA, .x))) %>% 
  mutate(best_location_code = case_when(collection_countryname %in% large_countries ~ coalesce(collection_subdiv1code,
                                                                                                collection_countrycode),
         .default  = collection_countrycode)) %>%
  
  # Sample one virus in each segment that has unique combination of location, host family, 
  # reasssortment, and sequence identity (viruses in the same group are identical)
  
  group_by(label, 
           cluster_profile
  ) %>%
  mutate(todrop = case_when(collection_countrycode == 'CN' & is.na(collection_subdiv1code) & sum(collection_countrycode == 'CN') > 1~ 'drop',
                            collection_countrycode == 'RU' & is.na(collection_subdiv1code) & sum(collection_countrycode == 'RU') > 1 ~ 'drop',
                            collection_countrycode == 'JP' & is.na(collection_subdiv1code) & sum(collection_countrycode == 'JP') > 1~ 'drop',
                            collection_countrycode == 'IRN' ~ 'drop',
                            .default = 'keep')) %>%
  filter(todrop == 'keep') %>%
  ungroup() %>%
  select(-c(6,34,35))


####################################################################################################
# Sample first and last sequence from reassortant
metadata_subsampled_1 <- metadata_prep  %>%
  group_by(label, 
           cluster_profile
  ) %>%
  
  # Sample first and last sequence from reassortant
  slice(which(collection_date == max(collection_date, na.rm = T) | collection_date == min(collection_date, na.rm = T)), 
        .preserve = T) %>%
  group_by(collection_date, .add = TRUE) %>%
  slice_sample(n = 1) %>%
  ungroup() %>%
  
  # as list
  group_split(label, .keep =  FALSE) %>%
  as.list() %>%
  setNames(names(metadata_per_alignment))


####################################################################################################
# sample according to location
metadata_subsampled_2 <- metadata_prep %>%
  group_by(label, 
           cluster_profile,
           sequence_group,
           best_location_code
  ) %>%
  slice_sample(n = 1) %>%
  ungroup() %>%
  
  # as list
  group_split(label, .keep = FALSE) %>%
  as.list() %>%
  setNames(names(metadata_per_alignment))


####################################################################################################
# sample according to host
metadata_subsampled_3 <- metadata_prep %>%
  group_by(label, 
           cluster_profile,
           sequence_group) %>% 
  slice(which(host_class == 'mammalia'), 
        .preserve = T) %>%
  slice_sample(n = 1) %>%
  ungroup() %>%
  
  full_join(tibble(label = names(metadata_per_alignment))) %>%
  
  # as list
  group_split(label, .keep = FALSE) %>%
  as.list() %>%
  setNames(names(metadata_per_alignment))

# Combined all subsamples together
metadata_subsampled <- lapply(1:length(metadata_per_alignment), function(i) bind_rows(metadata_subsampled_1[[i]],
                                                                                      metadata_subsampled_2[[i]],
                                                                                      metadata_subsampled_3[[i]]) %>%
                                distinct()) %>%
  setNames(names(metadata_per_alignment))

lapply(metadata_subsampled, nrow) %>% unlist()

####################################################################################################
# Subsample alignments
alignments_subsampled <- mapply(function(x,y) x[str_extract(rownames(x), "EPI_ISL_(china_){0,1}\\d+[^.|]*") %in% y$isolate_id,],
                                aln,
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



alignmentfiles_subsampled <- paste('./2024Nov18/reassortant_subsampled',
                                   paste(names(aln), 'subsampled.fasta', sep = '_'),
                                   sep = '/' )

mapply(ape::write.FASTA, 
       alignments_subsampled, 
       alignmentfiles_subsampled )


metadatafiles_subsampled_beast <-paste('./2024Nov18/reassortant_subsampled',
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
               gsub('.fasta$|.*beastsubsample/|./2024Nov18/reassortant_subsampled/', '', alignmentfiles_subsampled), 
               '_relaxLn_constant', '_2',
               "' flu_constanttemplate ",
               gsub('.*beastsubsample/|2024Nov18/', '', alignmentfiles_subsampled),
               ' ',
               gsub('.fasta|.*beastsubsample/|2024Nov18/', '', alignmentfiles_subsampled),
               '_relaxLn_constant_2', '.xml')

write_lines(cmds,  paste('./2024Nov18', 
                         'beastgen_reassortants.sh',
                         sep = '/' ))


