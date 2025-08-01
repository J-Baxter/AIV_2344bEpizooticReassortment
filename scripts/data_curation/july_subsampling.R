library(ape)
library(tidyverse)
library(phangorn)
#source('./scripts/FindIdenticalSeqs.R')

# Calculates pairwise hamming distance, then groups at a threshold difference
GroupSequences <- function(aln, snp_threshold = 0){
  
  require(phangorn)
  require(igraph)
  
  # Ensure alignment is correctly formatted
  if(class(aln) != 'PhyDat'){
    aln_formatted <- as.phyDat(aln) 
  }else{
    aln_formatted <- aln
  }
  
  # Calculate hamming distance
  hd_normalised <- dist.hamming(aln_formatted) %>%
    as.matrix()
  hd_raw <- hd_normalised * ncol(aln)
  
  # Obtain groups of sequences for which HD < SNP threshold
  
  groups <- which(hd_raw <= snp_threshold, 
                  arr.ind = TRUE) %>%
    dplyr::as_tibble(rownames = 'tipnames') %>% 
    filter(row != col) %>%
    dplyr::select(-tipnames) %>%
    
    # Infer network from HDs
    igraph::graph_from_data_frame(.,
                                  directed = F) %>%
    
    components() %>%
    getElement('membership') %>%
    stack() %>%
    as_tibble() %>%
    mutate(ind = as.numeric(as.character(ind))) %>%
    mutate(tipnames = map_chr(ind, ~ rownames(aln)[.x])) %>%
    dplyr::select(c(tipnames, values)) %>%
    dplyr::distinct() %>%
    dplyr::rename(sequence_group = values)
  
  out <- tibble(tipnames = rownames(aln)) %>%
    left_join(groups) %>%
    mutate(sequence_group = 
             ifelse(is.na(sequence_group), 
                    max(sequence_group, na.rm = T) + row_number() + 1, 
                    sequence_group))
  
  
  return(out)
}
####################################################################################################
# Import metadata and alignments
metadata <- read_csv('2024-09-09_meta.csv')

aln_files <- list.files(path = './2024Aug18/region_alignments',
                        pattern = 'fasta',
                        full.names = T) #%>%
  #.[grep('southamerica', .)]

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
groups <- lapply(aln, GroupSequences, snp_threshold = 5)%>%
  bind_rows(., .id = 'label') %>%
  mutate(isolate_id = str_extract(tipnames, "EPI_ISL_(china_){0,1}\\d+[^.|]*"))


#ml_trees <- list.files('./2023Dec02/ml_phylo', 
                      # pattern = 'treefile$',
                       #full.names = T) %>%
  #lapply(., read.tree) %>%
 # lapply(., function(x) {x$tip.label <- gsub('\\.', '|', x$tip.label) 
  #return(x)})

#ExtractTerminalBranchLengths <- function(tree, aln_length, snp_threshold = 0){
  
  #create tidy tree object
  #tidytree <- as_tibble(tree)
  
  #sequence_groups <- tidytree %>%
    #mutate(branch.length = ifelse(branch.length<0, 0, branch.length)) %>%
    #mutate(snps = floor(branch.length * aln_length)) %>% 
    #filter(snps <= snp_threshold & grepl('\\|', label)) %>% 
    #group_by(parent) %>%
    #mutate(n = n()) %>% 
    #filter(n > 1) %>%
    #mutate(sequence_group = cur_group_id()) %>%
    #ungroup() %>%
    #select(c(label, sequence_group)) %>%
    #rename(tipnames = label)
  
  
 # return(sequence_groups)
#}


#groups <- mapply(ExtractTerminalBranchLengths,
                 #ml_trees,
                 #lapply(aln, ncol),
                 #snp_threshold = 0,
                 #SIMPLIFY = F) %>%
  #setNames(names(aln)) %>%
  #bind_rows(., .id = 'label') %>%
  #mutate(isolate_id = str_extract(tipnames, "EPI_ISL_(china_){0,1}\\d+[^.|]*"))

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
                     'ethiopia',
                     'bolivia' ,
                     'mauritania',
                     'egypt',
                     'tanzania' ,
                     'nigeria' ,
                     'venezuela' ,
                     'pakistan' ,
                     'namibia' ,
                     'mozambique' ,
                     'turkey',
                     'chile' ,
                     'zambia' ,
                     'myanmar')


####################################################################################################
# Subsample - strict
metadata_prep <- metadata_per_alignment %>%
  bind_rows(., .id = 'label') %>%
  
  
  # join identical sequence groupings
  left_join(groups %>% select(-tipnames), 
            by = join_by(label, 
                         isolate_id)) %>%
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
  mutate(best_location_code = case_when(collection_countryname %in% large_countries ~ coalesce(collection_subdiv1code,
                                                                                               collection_countrycode),
                                        .default  = collection_countrycode)) %>%
  
  mutate(sequence_group = as.factor(sequence_group)) %>%
  
  # Sample one virus in each segment that has unique combination of location, host family, 
  # reasssortment, and sequence identity (viruses in the same group are identical)
  
  group_by(label, 
           sequence_group,
           #host_simplifiedhost, 
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
  group_split(label, .keep = FALSE) %>%
  as.list() %>%
  setNames(names(metadata_per_alignment))


####################################################################################################
# sample according to location 
metadata_subsampled_2 <- metadata_prep %>%
  group_by(label, 
           cluster_profile,
           sequence_group,
           best_location_code) %>%
  slice_sample(n = 1) %>%
  ungroup() %>%
  
  # as list
  group_split(label, .keep = FALSE) %>%
  as.list() %>%
  setNames(names(metadata_per_alignment))


####################################################################################################
# sample according to host (mammals including human)
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
                                                                                      metadata_subsampled_3[[i]]
                                                                                      ) %>%
                                distinct()) %>%
  setNames(names(metadata_per_alignment))

lapply(metadata_subsampled, nrow)



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



alignmentfiles_subsampled <- paste('./2024Nov18/region_subsampled',
                                   paste(names(aln), 'subsampled.fasta', sep = '_'),
                                   sep = '/' )

mapply(ape::write.FASTA, 
       alignments_subsampled, 
       alignmentfiles_subsampled)


metadatafiles_subsampled_beast <-paste('./2024Nov18/region_subsampled',
                                       paste(names(aln), 'subsampled.txt',  sep = '_'),
                                       sep = '/' )

mapply(write_delim, 
       delim = '\t',
       quote= 'needed',
       metadata_subsampled_beast, 
       metadatafiles_subsampled_beast)  


####################################################################################################
# Write plain BEAUTI XML and metadata to file

treeprior <- lapply(metadata_subsampled, function(x) x %>% filter(!is.na(collection_date)) %>% summarise(prior = (max(decimal_date(ymd(collection_date)), na.rm = T) - min(decimal_date(ymd(collection_date)), na.rm = T)))) %>% unlist() %>% ceiling()

cmds <- paste0("./beastgen -date_order -1 -date_prefix \\| -date_precision  -D ",
               "'skygrid_PopSize=",
               treeprior*4,
               ",skygrid_numGridPoints=",
               str_trim(format(treeprior*4-1, nsmall = 1)),
               ",skygrid_cutOff=",
               str_trim(format(treeprior, nsmall = 1)),
               ',fileName=',
               gsub('.prior', '', names(treeprior)), 
               '_relaxLn_Skygrid', treeprior, '-', treeprior*4, '_2',
               "' flu_skygridtemplate ",
               gsub('2024Nov18\\/' ,'' , alignmentfiles_subsampled),
               ' ',
               ' ./region_subsampled/',
               gsub('.prior', '', names(treeprior)), 
               '_relaxLn_Skygrid', treeprior, '-', treeprior*4, '_2', '.xml')

write_lines(cmds,  paste('.',
                         '2024Nov18', 
                         'beastgen_regions.sh',
                         sep = '/' ))
