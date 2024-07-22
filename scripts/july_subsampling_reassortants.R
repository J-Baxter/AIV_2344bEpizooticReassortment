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