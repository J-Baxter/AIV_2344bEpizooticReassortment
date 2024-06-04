# 
library(tidyverse)
load("/Users/s1506888/Downloads/H5_predictive_paper/h5nx_2344b_clusters_20240513.Rda") 
meta <- meta %>% as_tibble() %>%
  mutate(date_year = as.double(date_year))


split_by_region <- meta %>% group_split(location_1)

regions <- sapply(split_by_region, `[`, 'location_1') %>%
  sapply(., head, 1) %>%
  as.vector() %>%
  tolower() %>%
  gsub(' ', '', .)
names(split_by_region) = regions
# read HA alignment
library(ape)
# separate tree including only reassortants that originated in east asia.

ha_aln <- read.dna('~/Downloads/aln_ha.fasta',
                   as.matrix = T,
                   format = 'fasta')

pb2_aln <- read.dna('~/Downloads/aln_pb2.fasta',
                   as.matrix = T,
                   format = 'fasta')

ha_by_region <- lapply(split_by_region, function(x) ha_aln[gsub('\\|.*', '', rownames(ha_aln)) %in% x$isolate_id, ] )
names(ha_by_region) <- regions
pb2_by_region <- lapply(split_by_region, function(x) pb2_aln[gsub('\\|.*', '', rownames(pb2_aln)) %in% x$isolate_id, ] )
names(pb2_by_region) <- regions

ha_filenames <- paste0('~/Downloads/ha_', regions, '.fasta' )
mapply(write.dna,
       ha_by_region,
       ha_filenames,
       format = 'fasta')


pb2_filenames <- paste0('~/Downloads/pb2_', regions, '.fasta' )
mapply(write.dna,
       pb2_by_region,
       pb2_filenames,
       format = 'fasta')
