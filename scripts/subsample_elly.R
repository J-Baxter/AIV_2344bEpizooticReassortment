library(ape)
library(tidyverse)
source('./scripts/FindIdenticalSeqs.R')

####################################################################################################
# Import metadata and alignments
elly_reassortants <- c('4_5_1_1_2_1_1_3',
                       '8_5_1_1_2_1_1_3',
                       '1_1_4_1_4_1_1_4',
                       '1_1_2_1_1_1_1_1') 

clusters_to_downsample <- c('1_1_4_1_4_1_1_4', '1_1_2_1_1_1_1_1')

metadata <- read_csv('2024-08-19_meta.csv') %>%
  filter(cluster_profile %in% elly_reassortants) %>% #filter by required reassortants
  split(~cluster_profile)  %>%  # split the dataframe by the group column
  map_at(clusters_to_downsample, ~ group_by(.x,
                                            collection_regionname,
                                            collection_datemonth) %>%
           slice_sample(n = 2) %>%
           ungroup()) %>%  # downsample only the specified groups
  bind_rows() %>%
  distinct()


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
aln_elly <- lapply(aln_files ,
                        read.dna, 
                        format = 'fasta',
                        as.matrix = TRUE) %>%
  lapply(., function(x) x[str_extract(rownames(x), "EPI_ISL_(china_){0,1}\\d+[^.|]*") %in% metadata$isolate_id,])

names(aln_elly) <- str_split_i(aln_files, '/', 4) %>%
  gsub('h5_|.fasta', '', .)


# Split alignments by reassortant
SplitAlignment <- function(alignment, data){
  
  subset <- alignment[str_extract(rownames(alignment), "EPI_ISL_(china_){0,1}\\d+[^.|]*")  %in% data$isolate_id, ] 
  
  return(as.list(subset))
}


aln_split <- aln_elly %>%
  split(.,  str_split_i(names(aln_elly), '_', 1)) %>% #split list of alignments by segment
  map(., compact) %>% #remove empty matrices
  compact() %>% #remove empty lists
  map(., ~do.call(rbind.DNAbin, .x))%>%
  lapply(., function(x) x[unique(rownames(x)),]) #bind together within segments

elly_fileneames <- paste0('./2024Aug18/elly_subsample/', names(aln_split), '.fasta') 



mapply(ape::write.dna, 
       aln_split, 
       elly_fileneames ,
       format = 'fasta')