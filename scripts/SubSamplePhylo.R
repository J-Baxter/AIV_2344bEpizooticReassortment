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
# Import tree data
treefiles <- list.files(path = './data/phylo_ml/original_mlphylos',
                        recursive = FALSE,
                        include.dirs = FALSE, 
                        full.names = TRUE)

trees <- lapply(treefiles, read.newick) 


# Subset alignments
alignmentfiles <- list.files('./data/alignments/init_alignments',
                             full.names = T, 
                             pattern="\\.")

alignments <- lapply(alignmentfiles,
                     read.dna, 
                     format = 'fasta',
                     as.matrix = T)


# Segment names
segnames <- str_split(treefiles,  '_') %>% 
  lapply(., tail, n = 1) %>% 
  unlist() %>%
  toupper() %>%
  gsub('NA', 'N', .)

names(alignments) <- segnames
names(trees) <- segnames

# Import existing metadata
reassortant_metadata <- read_csv('./data/metadata/h5_metadata.csv') %>%
  mutate(date = dmy(date) %>%
           as.Date())


# Import list of sequences to remove due to lack of temporal signal
seqstoremove <- read_csv('./data/erroneous_seqs.csv') %>%
  mutate(drop.temporal = 'drop') %>%
  setNames(c('virus.segment', 'tipnames', 'drop.temporal'))


####################################################################################################
# Extract metadata and format tipnames
tipnames <- lapply(trees, TipLabels) %>% 
  setNames(segnames) # Only required because we are working backwards from tipnames


# Extract metadata from tipnames
metadata_unformatted <- lapply(tipnames, ExtractMetadata)

# Format metadata
reassortant_metadata_formatted <- FormatMetadata(reassortant_metadata) %>%
  select(-matches('location.[cn]')) %>%
  select(-c(domestic.status, 
            domestic.wild)) %>% 
  unite('collection.original', 
        starts_with('location'), 
        sep = ' / ')

metadata_formatted <- lapply(metadata_unformatted, FormatMetadata) 

# Join phylo metadata with reassortant data
MyFunc <- function(data, newdata){
  df <- data  %>% 
    left_join(newdata, by = join_by(isolate.id)) %>%
    mutate(across(ends_with(".x"), ~coalesce(., get(sub("\\.x$", ".y", cur_column()))), .names = "{.col}")) %>%
    select(-ends_with(".y")) %>%
    rename_with(~gsub('.x', '', .x)) %>%
    select(-c(id.unsure, week.date, is.problem.bird, virus.species)) %>%
    mutate(collection.tipdate = case_when(is.na(collection.date) ~ collection.datemonth,
                                          .default = as.character(collection.date))) 
  
  return(df)
}

metadata_joined <- lapply(metadata_formatted, MyFunc, reassortant_metadata_formatted ) %>%
  setNames(segnames) %>%
  lapply(., MakeTipNames)
  

metadata_joined_df <- metadata_joined %>% 
  bind_rows(., .id = 'virus.segment')

####################################################################################################
# Rename phylogenies and alignments with new-format tipnames
renamed_phylos <- lapply(segnames,
                         ReNamePhylo, 
                         trees = trees, 
                         metadata = metadata_joined)

mapply(write.tree,
       renamed_phylos,
       file = paste0(segnames,  '2.tree'))



renamed_alignments <- alignments %>%
  mapply(ReNameAlignment, 
         .,
         metadata_joined , 
         SIMPLIFY = FALSE)

mapply()
####################################################################################################
# Identify duplicate sequences from alignment

identical_seqs <- lapply(renamed_alignments,
                         FindIdenticalSeqs) %>%
  bind_rows(., .id = 'virus.segment')


####################################################################################################
# Subsample - strict
metadata_subsampled <- metadata_joined_df %>%
  
  #drop seqs
  left_join(., seqstoremove) %>% 
  filter(is.na(drop.temporal)) %>%
  select(-drop.temporal) %>%
  
  # join identical sequence groupings
  left_join(identical_seqs, 
            by = join_by(virus.segment, 
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
  group_by(virus.segment, 
           group,
           host.family, 
           cluster.genome,
           best_location_code) %>%
  slice_sample(n = 1)   %>%
  ungroup() %>%
  
  # as list
  group_split(virus.segment, .keep = FALSE) %>%
  as.list() %>%
  setNames(segnames)


####################################################################################################
# Subsample alignments
alignments_subsampled <- renamed_alignments %>%
  mapply(SubsampleAlignment,
         .,
         metadata_subsampled,
         SIMPLIFY = FALSE) 
####################################################################################################
#BEAST traits
metadata_subsampled_beast <- metadata_subsampled %>%
  lapply(., function(x) x %>% 
           select(c(tipnames, virus.subtype, contains('collection'), cluster.profile, host.order)) %>%
           mutate(tiplocation_discrete = coalesce(collection.subdiv2.code, collection.subdiv1.code, collection.country.code)) %>%
           mutate(tiplocation_lon = coalesce(collection.subdiv2.long , collection.subdiv1.long, collection.country.long)) %>%
           mutate(tiplocation_lat = coalesce(collection.subdiv2.lat , collection.subdiv1.lat, collection.country.lat)))

####################################################################################################
# Write alignments and metadata to file
alignmentfiles_subsampled <- alignmentfiles  %>%
  gsub('Re_H5_', '',. ) %>%
  gsub('_Asia_blast_', '_asia_',. ) %>%
  gsub('.trim.fasta$', '_subsampled.fasta',. ) %>%
  gsub('/init_alignments/', '/subsampled_alignments/2024Jan24/',. )

mapply(write.dna, 
       alignments_subsampled, 
       alignmentfiles_subsampled ,
       format = 'fasta')


metadatafiles_subsampled <- paste('./data/metadata/2024Jan24/', 
                                  segnames, 
                                  '_subsampled.tsv',  
                                  sep = '')
mapply(write_tsv, 
       metadata_subsampled_beast, 
       metadatafiles_subsampled)

  
####################################################################################################
# END #
####################################################################################################