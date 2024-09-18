#dependencies
library(ape)
library(tidyverse)
library(zoo)

# create beast txt files from subsampled alignments
ImputeCladeandCluster <- function(metadata, alignment, ordered = FALSE){
  
  seg <- metadata %>% 
    pull(segment) %>%
    unique()
  
  print(seg)
  
  if(ordered == FALSE){
    tree = nj(dist.dna(alignment))
    seqnames = tree$tip.label
    
  }else{
    seqnames = rownames(alignment)
  }
  
  if(!any(grepl('cluster_[a-z]{2}\\d{0,1}$', colnames(metadata)))){
    metadata <-  metadata %>% 
      separate_wider_delim(cluster_profile, '_', 
                           names = paste('cluster',
                                         tolower(c("PB2", "PB1", "PA", "HA", "NP", "NA", "MP", "NS")), 
                                         sep = '_'),
                           cols_remove = FALSE)
  }
  
  out <- metadata %>%
    mutate(tipnames = ordered(tipnames,
                              levels = seqnames)) %>%
    arrange(tipnames) %>%
    
    # Impute clade
    mutate(clade = case_when(
      is.na(clade) & na.locf0(clade, fromLast = TRUE) == na.locf0(clade, fromLast = FALSE) ~ na.locf0(clade, fromLast = TRUE), 
      .default = clade)) %>%
    
    # Impute cluster (column is dependent on segment)
    rename(profile = cluster_profile) %>%
    #rename(genome = cluster_genome) %>%
    pivot_longer(matches('cluster_[a-z]{2}\\d{0,1}$'), 
                 values_to = 'cluster_number',
                 names_to = 'cluster_segment') %>%
    mutate(cluster_segment = gsub('.*_', '',
                                  cluster_segment)) %>%
    mutate(cluster_segment = case_when(
      grepl('^N[:0-9:]{0,1}$', segment, ignore.case = T) & cluster_segment == 'na' ~ tolower(segment),
      .default = cluster_segment)) %>% 
    filter(tolower(segment) == cluster_segment) %>%
    mutate(cluster_number = case_when(
      is.na(cluster_number) & na.locf0(cluster_number, fromLast = TRUE) == na.locf0(cluster_number, fromLast = FALSE) ~ na.locf0(cluster_number, fromLast = TRUE), 
      .default = cluster_number))
  
  
  return(out)
}


# import alignments
aln_files <- list.files('./2024Sept16/reassortant_subsampled_alignments',
                        pattern = '.fasta',
                        full.names = T)

aln <- lapply(aln_files, read.dna, format = 'fasta', as.matrix = T) 

names(aln) <- str_split_i(aln_files, '/', 4) %>%
  gsub('h5_|.fasta|_subsampled', '', .) 


# import metadata
data <- read_csv('2024-09-09_meta.csv')


# create subsampled dataframes
metadata_per_alignment <- lapply(aln, function(x) data %>% 
                                   filter(isolate_id %in% str_extract(rownames(x), "EPI_ISL_(china_){0,1}\\d+[^.|]*"))) 


# sanity check
stopifnot(all(unlist(lapply(aln, nrow)) == unlist(lapply(metadata_per_alignment, nrow))))

names(metadata_per_alignment) <- str_split_i(aln_files, '/', 4) %>%
  gsub('h5_|.fasta|_subsampled', '', .) 

# infer cluster number --------
test <-  metadata_per_alignment %>%
  mapply(function(x, y) x %>% mutate(segment = gsub('_.*', '', y)) ,x= .,  y= as.list(names(aln)), SIMPLIFY = F) %>%
  mapply(ImputeCladeandCluster,
               .,
               aln,
         ordered = FALSE,
         SIMPLIFY = F)

# sanity checks
stopifnot(all(unlist(lapply(aln, nrow)) == unlist(lapply(test, nrow))))

mapply(function(alignment, data) all(rownames(alignment) %in% data$tipnames),
       aln,
       test,
       SIMPLIFY = F)


# output to txt file
metadata_subsampled_beast <- lapply(metadata_per_alignment, 
                                    function(x) x %>% 
                                      select(c(tipnames,
                                               virus_subtype, 
                                               collection_regionname, 
                                               collection_countryname,
                                               ends_with('long'),
                                               ends_with('lat'),
                                               #cluster_number,
                                               host_simplifiedhost)) %>%
                                      mutate(lat = coalesce(collection_subdiv1lat, 
                                                            collection_countrylat)) %>%
                                      mutate(long = coalesce(collection_subdiv1long, 
                                                             collection_countrylong)) %>%
                                      mutate(across(c(lat, long), .fns = ~ as.numeric(.x))) %>%
                                      #mutate(cluster_number = str_pad(cluster_number, 4, pad = "0") %>% paste0('profile',.)) %>%
                                      select(where(~n_distinct(.) > 1)) %>%
                                      select(-c(contains('date'),
                                                contains('subdiv'),
                                                collection_countrylat,
                                                collection_countrylong,
                                                #collection_original, 
                                                #collection_tipdate
                                      ))) 


# write to file

metadatafiles_subsampled_beast <-paste('./2024Sept16/reassortant_subsampled_traits',
                                       paste(names(aln), 'subsampled.txt',  sep = '_'),
                                       sep = '/' )

mapply(write_delim, 
       delim = '\t',
       quote= 'needed',
       metadata_subsampled_beast, 
       metadatafiles_subsampled_beast)  
