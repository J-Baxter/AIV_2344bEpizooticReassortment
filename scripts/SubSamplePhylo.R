####################################################################################################
# Basic filtering of maximum likelihood phylogenies

# Requires a phylogeny (branch lengths inferred from genetic distance only) and accompanying 
# alignment

####################################################################################################
# Package dependencies
library(ape)
library(tidyverse)
library(treeio)
library(TreeTools)
library(Rcpp)
library(igraph)

source('./scripts/FindIdenticalSeqs.R')
sourceCpp("./scripts/getTaxonomyForName.cpp") # ~8mins
sourceCpp('./scripts/getLocation.cpp') #~3mins




MakeTipNames <- function(data){
  out <- data %>%
    unite(tipnames, 
          virus.subtype,
          clade,
          isolate.id,
          host.order,
          collection.country.name,
          cluster.genome,
          collection.tipdate,
          sep = '|',
          remove = FALSE) %>%
    mutate(tipnames = gsub(' ', '_', tipnames))
  
  return(out)
}


ReNamePhylo <- function(trees, metadata, x){
  
  trees[[x]][['tip.label']] <- metadata[[x]][['tipnames']]
  
  return(trees[[x]])
}


ReNameAlignment <- function(alignment, data){
  isolates <- str_extract(rownames(alignment), "([^|]*)\\|")%>%
    str_replace_all("\\|", "")
  
  new_seqnames <- sapply(isolates, function(x) data$tipnames[data$isolate.id %in% x]) %>% 
    as.vector() 
  
  rownames(alignment) <-  new_seqnames
  
  return(alignment)
}


SubsampleAlignment <- function(alignment, data, removepipe = TRUE){
  out <- alignment[rownames(alignment) %in% data$tipnames,]
  
  if(removepipe == TRUE){
    rownames(out) <- gsub('\\|', '.', rownames(out)) 
  }else{
    rownames(out) <- rownames(out)
  }
  
  return(out)
}


####################################################################################################
# Import alignments
alignmentfiles <- list.files('./2023Dec02/alignments',
                             full.names = T, 
                             pattern="asia")

alignments <- lapply(alignmentfiles,
                     read.dna, 
                     format = 'fasta',
                     as.matrix = T)

#import metadata
metadatafiles <-  list.files('./2023Dec02/metadata',
                             full.names = T, 
                             pattern="asia")


metadata <- lapply(metadatafiles,
                   read_csv,
                   col_types = cols(collection.date = col_date(format = "%Y-%m-%d"),
                                    collection.tipdate = col_character(),
                                    `...26` = col_skip())) 


# Segment names
segnames <- str_split(treefiles,  '_') %>% 
  lapply(., tail, n = 3) %>% 
  lapply(., `[`, 1:2) %>%
  lapply(., paste, collapse = '_') %>%
  unlist() %>%
  tolower()

names(alignments) <- segnames
names(metadata) <- segnames


# Import list of sequences to remove due to lack of temporal signal
seqstoremove <- read_csv('./2023Dec02/unclocklikeseqs.csv') %>%
  select(c(1,2,3)) %>%
  mutate(drop.temporal = 'drop') %>%
  mutate(seqname = gsub('*2.3.4.4.b*', '2344b', seqname)) %>%
  setNames(c('collection.region', 'virus.segment', 'tipnames', 'drop.temporal')) %>%
  unite(id, collection.region, virus.segment)


####################################################################################################

####################################################################################################
# Identify duplicate sequences from alignment

identical_seqs <- lapply(alignments,
                         FindIdenticalSeqs) %>%
  bind_rows(., .id = 'id')


####################################################################################################
# Subsample - strict
metadata_subsampled <- metadata %>%
  bind_rows(., .id = 'id') %>%
  
  #drop seqs
  left_join(., seqstoremove) %>% 
  filter(is.na(drop.temporal)) %>%
  select(-drop.temporal) %>%
  
  mutate(collection.datemonth = ym(collection.datemonth)) %>%
  
  # join identical sequence groupings
  left_join(identical_seqs, 
            by = join_by(id, 
                         tipnames)) %>%
  # create groupings
  mutate(across(contains('collection'), .fns = ~ gsub('^NA$', NA, .x))) %>% 
  mutate(best_location_code = coalesce(collection.subdiv1.code,
                                       collection.country.code)) %>%
  mutate(group = as.factor(group)) %>%
  
  # Remove seqs in RU/JP/CN that do not have subdivision data (unless it is a unique seq)
  #group_by(virus.segment, 
  # group) %>%
  #mutate(todrop = case_when(collection.country.code %in% c('RU', 'JP', 'CN') && is.na(collection.subdiv1.code) && n()>1 ~ 'drop',
  # .default = 'keep')) %>%
  #filter(todrop == 'keep') %>%
  #ungroup() %>%
  
  group_by(host.order) %>%
  mutate(todrop = case_when(n()<10 ~ 'drop',
                            .default = 'keep')) %>%
  filter(todrop == 'keep') %>%
  ungroup() %>%
  
  # Sample one virus in each segment that has unique combination of location, host family, 
  # reasssortment, and sequence identity (viruses in the same group are identical)
  group_by(id, 
           group,
           host.family, 
           genome,
           best_location_code) %>%
  slice_sample(n = 1)   %>%
  ungroup() %>%
  
  # as list
  group_split(id, .keep = FALSE) %>%
  as.list() %>%
  setNames(segnames)


####################################################################################################
# Subsample alignments
alignments_subsampled <- alignments %>%
  mapply(SubsampleAlignment,
         .,
         metadata_subsampled,
         SIMPLIFY = FALSE) 
####################################################################################################
#BEAST traits
metadata_subsampled_beast <- metadata_subsampled %>%
  lapply(., function(x) x %>% 
           select(c(tipnames, virus.subtype, contains('collection'), cluster.number, host.order)) %>%
           #mutate(tiplocation_discrete = coalesce(collection.subdiv2.code, collection.subdiv1.code, collection.country.code)) %>%
           mutate(tiplocation_lon = coalesce(collection.subdiv1.long, collection.country.long)) %>%
           mutate(tiplocation_lat = coalesce(collection.subdiv1.lat, collection.country.lat)))

####################################################################################################
# Write alignments and metadata to file

ddmonthyy <- format(Sys.Date(), '%Y%b%d')
check_dirs <- paste('.', ddmonthyy,c('metadata', 'alignments'),  sep = '/')
dirs <- list.dirs()

for (check_dir in check_dirs){
  
  if (!(check_dir %in% dirs)){
    dir.create(check_dir, recursive = T)
  }
  
}


alignmentfiles_subsampled <- paste('.',
                                   ddmonthyy, 
                                   'alignments',
                                   paste(segnames, '2344b_subsampled.fasta', sep = '_'),
                                   sep = '/' )

mapply(write.dna, 
       alignments_subsampled, 
       alignmentfiles_subsampled ,
       format = 'fasta')


metadatafiles_subsampled <-paste('.',
                                 ddmonthyy, 
                                 'metadata',
                                 paste(segnames, '2344b_subsampled.txt',  sep = '_'),
                                 sep = '/' )

mapply(write_delim, 
       delim = '\t',
       quote= 'needed',
       metadata_subsampled_beast, 
       metadatafiles_subsampled)  
####################################################################################################
# END #
####################################################################################################