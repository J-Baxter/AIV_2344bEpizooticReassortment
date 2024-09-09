# Split regional alignments by dominant subtype
# NB useful list manipulation using purr (lines 62 - 66)

# Dependencies
library(tidyverse)
library(ape)


SplitAlignment <- function(alignment, data){
  
  subset <- alignment[str_extract(rownames(alignment), "EPI_ISL_(china_){0,1}\\d+[^.|]*")  %in% data$isolate_id, ] 
  
  return(as.list(subset))}

# Import metadata

new_clusters <- c('7_1_5_2_1_3_1_2',
                  '5_4_9_1_2_1_1_1',
                  '1_6_2_1_1_1_1_1')

metadata <- read_csv('2024-08-19_meta.csv') %>%
  filter(cluster_profile %in% new_clusters) 


aln_files <- list.files(path = './2024Aug18/region_alignments',
                        pattern = 'fasta',
                        full.names = T) 

aln <- lapply(aln_files ,
              read.dna, 
              format = 'fasta',
              as.matrix = TRUE)

names(aln) <- str_split_i(aln_files, '/', 4) %>%
  gsub('h5_|.fasta', '', .) 

####################################################################################################
# subset alignment

# import alignments and filter by required clusters (ie those filtered above)
aln_new <- lapply(aln_files ,
                   read.dna, 
                   format = 'fasta',
                   as.matrix = TRUE) %>%
  lapply(., function(x) x[str_extract(rownames(x), "EPI_ISL_(china_){0,1}\\d+[^.|]*") %in% metadata$isolate_id,])

names(aln_new) <- str_split_i(aln_files, '/', 4) %>%
  gsub('h5_|.fasta', '', .)


# Split alignments by reassortant
SplitAlignment <- function(alignment, data){
  
  subset <- alignment[str_extract(rownames(alignment), "EPI_ISL_(china_){0,1}\\d+[^.|]*")  %in% data$isolate_id, ] 
  
  return(as.list(subset))
}


aln_split <- aln_new %>%
  split(.,  str_split_i(names(aln_new), '_', 1)) %>% #split list of alignments by segment
  map(., compact) %>% #remove empty matrices
  compact() %>% #remove empty lists
  lapply(., function(x) lapply(x, as.list.DNAbin)) %>%
  map(., ~do.call(c, .x)) %>%
  lapply(., function(x) split(x, str_extract(names(x), '[0-9]_[0-9]_[0-9]_[0-9]_[0-9]_[0-9]_[0-9]_[0-9]'))) %>%
  map(., compact) %>% #remove empty matrices
  list_flatten() %>%
  lapply(., function(x) {names(x) <- gsub('.*\\.', '', names(x))
  return(x)}) %>%
  lapply(., function(x) x[unique(names(x))]) 


# write to file
filenames <- names(aln_split) %>%
  str_remove_all(., '(?<=_\\d)_') %>%
  
  paste0('./2024Aug18/reassortant_alignments/h5_',., '.fasta' ) 

mapply(write.dna,
       aln_split,
       filenames,
       format = 'fasta')

