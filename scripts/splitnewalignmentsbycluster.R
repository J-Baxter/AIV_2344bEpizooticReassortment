# Split regional alignments by dominant subtype
# NB useful list manipulation using purr (lines 62 - 66)

# Dependencies
library(tidyverse)
library(ape)


SplitAlignment <- function(alignment, data){
  
  subset <- alignment[str_extract(rownames(alignment), "EPI_ISL_(china_){0,1}\\d+[^.|]*")  %in% data$isolate_id, ] 
  
  return(as.list(subset))
}

# Import metadata
metadatafiles <- list.files(path = './2024Jul12/region_metadata',
                            full.names = TRUE)

summary_data <- read_csv( '2024-06-05_reassortant_summary.csv')

new_clusters <- meta %>% select(c(isolate_id, cluster_profile))

metadata_dominant <- lapply(metadatafiles, read_csv, col_types = cols(collection_tipdate = col_character())) %>% 
  bind_rows() %>%
  select(-c(cluster_profile, clade, cluster_number, date_frac)) %>%
  left_join(new_clusters) %>%
  left_join(summary_data, by = join_by(cluster_profile)) %>%
  filter(group == 'dominant') %>%
  distinct()

check_meta <- meta %>% 
  left_join(summary_data, by = join_by(cluster_profile)) %>%
  filter(group == 'dominant') %>%
  summarise(n = n(), .by = cluster_profile)

all(unique(meta$isolate_id) %in% metadata_dominant$isolate_id)  # Check all isolate ID are present in metadata


# Import alignments
aln_files <- list.files(path = './2024Jul12/region_alignments',
                        pattern = '.fasta',
                        include.dirs = FALSE,
                        full.names = TRUE)

aln_all <- lapply(aln_files, 
                  read.dna, 
                  as.matrix = T, 
                  format = 'fasta')

aln_isolates <- lapply(aln_all, function(x) str_extract(rownames(x), "EPI_ISL_(china_){0,1}\\d+[^.|]*"))
all(unique(meta$isolate_id) %in% unlist(unique(aln_isolates)))  # Check all isolate ID are present in alignments

# set cluster names
clusters <- summary_data %>% 
  select(c(cluster_profile, group)) %>%
  filter(group == 'dominant') %>%
  pull(cluster_profile) %>%
  unique() %>%
  sort() %>%
  tolower() %>%
  gsub('_', '', .)


# set segment names
segments <- aln_files %>%
  str_split_i(., '/', 4) %>% gsub('h5_|.fasta', '',.)


# Split dataframe by region
split_by_cluster <- metadata_dominant %>%
  group_split(cluster_profile) %>%
  setNames(clusters)


# Split alignments by reassortant
aln_split <- lapply(aln_all, function(x) lapply(split_by_cluster, SplitAlignment, alignment = x)) %>%
  split(.,  str_split_i(segments, '_', 1)) %>%
  lapply(., function(x)  pmap(.l = x, .f = c)) %>%
  list_flatten() %>%
  compact() %>%
  lapply(., function(x) x[unique(names(x))])


# write to file
filenames <- names(aln_split) %>%
  paste0('./2024Jul12/reassortant_alignments/h5_',., '.fasta' ) %>%
  sort()

mapply(write.dna,
       aln_split,
       filenames,
       format = 'fasta')


# metadata
reassortant_metadata <- lapply(aln_split, function(x) metadata_dominant %>% filter(isolate_id %in% str_extract(names(x), "EPI_ISL_(china_){0,1}\\d+[^.|]*")))

filenames <- names(aln_split) %>%
  paste0('./2024Jul12/reassortant_metadata/h5_',., '.csv' ) %>%
  sort()

mapply(write_csv,
       reassortant_metadata,
       filenames)