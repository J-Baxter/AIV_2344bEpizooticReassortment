library(tidyverse)
library(ape)


SplitAlignment <- function(alignment, data){
  seqnames <- gsub('\\|.*', '', rownames(alignment))
  subset <- alignment[seqnames  %in% data$isolate_id, ] 
  
  return(subset)
}



# Import metadata
load("./2024Jun01/h5nx_2344b_clusters_20240513.Rda")
meta <- meta %>% 
  as_tibble() %>%
  mutate(date_year = as.double(date_year),
         location_1 = case_when(location_1 == 'Antarctica' ~ 'South America', 
                                .default = location_1))


# Import alignments
aln_files <- list.files(path = './2024Jun01/master',
                        pattern = '.fasta',
                        include.dirs = FALSE,
                        full.names = TRUE)

aln_all <- lapply(aln_files, 
                  read.dna, 
                  as.matrix = T, 
                  format = 'fasta')


# set region names
regions <- meta %>% 
  pull(location_1) %>%
  unique %>%
  sort() %>%
  tolower() %>%
  gsub(' ', '', .)

# set segment names
segments <- aln_files %>%
  str_split_i(., '/', 4) %>%
  gsub('aln_|.fasta', '', .) 


# Split dataframe by region
split_by_region <- meta %>%
  group_split(location_1) %>%
  setNames(regions)

# Split alignments by regions
aln_split <- lapply(aln_all, function(x) lapply(split_by_region, SplitAlignment, alignment = x)) %>%
  setNames(segments) %>%
  flatten()

# write to file
filenames <- apply(expand.grid(segments, regions), 1, paste, collapse="_") %>%
  paste0('./2024Jun01/h5_',., '.fasta' ) %>%
  sort()

mapply(write.dna,
       aln_split,
       filenames,
       format = 'fasta')


# separate tree including only reassortants that originated in east asia.
