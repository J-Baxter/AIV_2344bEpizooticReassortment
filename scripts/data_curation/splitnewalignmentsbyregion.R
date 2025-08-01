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


# Import alignments (all except na)
aln_files <- list.files(path = './2024Jun01/master',
                        pattern = '.fasta',
                        include.dirs = FALSE,
                        full.names = TRUE) %>%
  .[!grepl('.*na.fasta',.)]

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


############ Seperate protocol for splitting neuraminidase ########
aln_na <- read.dna('./2024Jun01/master/aln_na.fasta',  as.matrix = T, 
                   format = 'fasta') 


# Split dataframe by region
split_by_region_na <- meta %>%
  mutate(location_1 = gsub(' ', '', tolower(location_1)),
         subtype = gsub('h[:0-9:]','', tolower(subtype))) %>%
  unite(., na_split, subtype, location_1 ) %>%
  split(., .$na_split)

# Split alignments by regions
aln_region_split_na <- lapply(split_by_region_na, SplitAlignment, alignment = aln_na) %>%
  setNames(names(split_by_region_na)) 



# write to file
filenames <- paste0('./2024Jun01/h5_',names(split_by_region_na), '.fasta' ) 

mapply(write.dna,
       aln_region_split_na,
       filenames,
       format = 'fasta')




# separate tree including only reassortants that originated in east asia.
