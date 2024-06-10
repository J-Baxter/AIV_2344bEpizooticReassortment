# dependencies 
library(ape)
library(tidyverse)

# load old alignments and metadata
old_aln_files <- list.files(path = './2023Dec02/alignments',
                            full.names = TRUE,
                            pattern = '.fasta')
nom <- str_split_i(old_aln_files, '/', 4) %>%
  gsub('_clade2344b.*', '', .)

old_aln <- lapply(old_aln_files, read.dna, format = 'fasta', as.matrix = TRUE)

old_meta_files <- list.files(path = './2023Dec02/metadata',
                            full.names = TRUE,
                            pattern = '.csv')

old_meta <- lapply(old_meta_files, read_csv, col_types = cols(collection_tipdate = col_character())) 


# filter metadata by missing cluster profile (because this cannot be inferred)
old_meta_blastonly <- lapply(old_meta, function(x) x %>% filter(is.na(profile)))

# subsample alignments, leaving a list of alignments stratified by region and segment
old_aln_split <- mapply(SplitAlignment, old_aln, old_meta_blastonly, SIMPLIFY = FALSE) %>%
  setNames(nom)


SplitAlignment <- function(alignment, data){
  seqnames <- str_split_i(rownames(alignment), '\\.', 3)
  subset <- alignment[seqnames  %in% data$isolate_id, ] 
  
  return(subset)
}

# import new alignments
aln_files <- list.files(path = './2024Jun01',
                        pattern = '[:a-z:].fasta',
                        include.dirs = FALSE,
                        full.names = TRUE) 

aln_all <- lapply(aln_files, 
                  read.dna, 
                  as.matrix = T, 
                  format = 'fasta')


names(old_aln_split) <- nom
names(aln_all) <- str_split_i(aln_files, '/', 3) %>%
  gsub('*h5_|\\..*', '', .)


# append blast sequences to new alignments
matched_aln <- list()
newnames <- c()
for(i in 1:length(aln_all)){
  if(names(aln_all)[[i]] %in% names(old_aln_split)){
    matching <- which(names(aln_all)[[i]] == names(old_aln_split))
    matched_aln[[i]] <- c(as.list(aln_all[[i]]), as.list(old_aln_split[[matching]]))
    newnames[i] <- names(aln_all)[[i]]
  }
}


# export
filenames <- paste0('./2024Jun01/region_alignments_withBLAST/h5_',newnames[-which(!lengths(matched_aln))], '.fasta' ) 

mapply(write.dna,
       s,
       filenames,
       format = 'fasta')






