################################################################################
## Script Name:        <INSERT_SCRIPT_NAME_HERE>
## Purpose:            Reformat initial alignments and phylogenies. Requires a 
##                     phylogeny (branch lengths inferred from genetic distance 
##                     only) and accompanying alignment
## Author:             James Baxter
## Date Created:       
################################################################################

############################### SYSTEM OPTIONS #################################
options(
  scipen = 6,     # Avoid scientific notation
  digits = 7      # Set precision for numerical display
)
memory.limit(30000000)

############################### DEPENDENCIES ###################################
# Load required libraries
library(tidyverse)
library(magrittr)
library(ape)
library(treeio)
library(TreeTools)
library(igraph)
library(zoo)
source('./scripts/FindIdenticalSeqs.R')
source("./scripts/Subsamplefunctions.R")
source('./scripts/FormatBirds.R')
source('./scripts/FormatMammals.R')


MakeTipNames <- function(data){
  out <- data %>%
    unite(tipnames, 
          virus_subtype,
          clade,
          isolate_id,
          host_order,
          collection_countryname,
          cluster_profile,
          collection_tipdate,
          sep = '|',                        
          remove = FALSE) %>%
    mutate(tipnames = gsub(' ', '_', tipnames))
  
  return(out)
}


isDateError <- function(dataframe){
  out <- dataframe %>%
    mutate(collection_dateerror = case_when(
      is.na(collection_tipdate) ~ TRUE,
      as.numeric(gsub('-.*', '', collection_tipdate))>2024 ~ TRUE,
      as.numeric(gsub('-.*', '', collection_tipdate))<1990 ~ TRUE,
      .default = FALSE))
  
  return(out)
}

ReNamePhylo <- function(trees, metadata, x){
  
  trees[[x]][['tip.label']] <- metadata[[x]][['tipnames']]
  
  return(trees[[x]])
}


ReNameAlignment <- function(alignment, data){
  z <- rownames(alignment)
  
  isolates <-  regmatches(z, gregexpr("EPI_ISL_(china_){0,1}\\d+[^.|]*", z)) %>% unlist()
  
  new_seqnames <- sapply(isolates, function(x) data$tipnames[data$isolate_id %in% x]) %>% 
    as.vector() 
  
  rownames(alignment) <-  new_seqnames
  
  return(alignment)
}


################################### DATA #######################################
# Read and inspect data
geodata <- read_csv('./annotated_geodata.csv')
birds <- read_csv('bird_taxonomy.csv')
mammals <- read_csv('mammal_taxonomy.csv')
cluster_classes <- read_csv('clusterprofile_summary.csv')


alignmentfiles <- list.files('./2023Dec01/alignments',
                             full.names = T, 
                             pattern="\\.")

alignments <- lapply(alignmentfiles,
                     read.dna, 
                     format = 'fasta',
                     as.matrix = T)


# Import existing metadata
reassortant_metadata <- read_csv('./2023Dec01/metadata/h5_metadata_global_6280_update.csv') %>%
  mutate(date = dmy(date) %>%
           as.Date())


################################### MAIN #######################################
# Main analysis or transformation steps

# Segment names
segnames <- str_split(alignmentfiles,  '_') %>% 
  lapply(., tail, n = 3) %>% 
  lapply(., `[`, c(1,3)) %>%
  lapply(., function(x) gsub("\\..*$", "",x)) %>%
  lapply(., paste0 , collapse = '_') %>%
  unlist() %>%
  tolower() %>%
  gsub('na', 'n', .)

names(alignments) <- segnames

# Extract metadata and format tipnames
seqnames <- lapply(alignments, rownames) %>% 
  setNames(segnames) # Only required because we are working backwards from tipnames


# Extract metadata from tipnames
metadata_unformatted <- lapply(seqnames, ExtractMetadata)

# Format metadata
reassortant_metadata_formatted <- FormatMetadata(reassortant_metadata) %>%
  select(-matches('location.[cn]')) %>%
  select(-c(domestic_status, 
            domestic_wild)) %>% 
  unite('collection_original', 
        starts_with('location'), 
        sep = ' / ') %>%
  mutate(clade =  gsub('[[:punct:]]+','', clade)) %>%
  left_join(cluster_classes, by = join_by(cluster_profile)) %>%
  mutate(cluster_profile = case_when(
    collection_regionname == 'eastern asia' & cluster_profile == '1_1_1_1_1_1_1_1' ~ '1_1_1_1_1_1_1_1A',
    collection_regionname != 'eastern asia' & cluster_profile == '1_1_1_1_1_1_1_1' ~ '1_1_1_1_1_1_1_1B',
    .default = cluster_profile))


metadata_formatted <- lapply(metadata_unformatted, FormatMetadata) 


metadata_joined <- lapply(metadata_formatted, 
                          MergeReassortantData,
                          reassortant_metadata_formatted) %>%
  setNames(segnames) %>%
  lapply(., MakeTipNames) %>%
  mapply(function(x, y) x %>% 
           mutate(segment = gsub('.*_', '', y)) 
         ,x= ., 
         y= as.list(segnames), 
         SIMPLIFY = F) %>%
  lapply(., isDateError)


# impute clade and cluster (from NJ tree)
temp_alignments <- alignments %>%
  mapply(ReNameAlignment, 
         .,
         metadata_joined , 
         SIMPLIFY = FALSE)

metadata_joined_imputed <- mapply(ImputeCladeandCluster,  
                                  metadata_joined, 
                                  temp_alignments, 
                                  SIMPLIFY = F) %>%
  lapply(., function(x) x %>% 
           mutate(tipnames = gsub( '\\|', '\\.',  tipnames))) %>%
  lapply(., function(x) x %>% 
           mutate(cluster_number = paste0('profile', str_pad(cluster_number, 3, pad = "0"))))


# Rename alignments with new-format tipnames
renamed_alignments <- alignments %>%
  mapply(ReNameAlignment, 
         .,
         metadata_joined_imputed, 
         SIMPLIFY = FALSE)


# Check year and remove if problematic
metadata_joined_imputed_dateschecked <- metadata_joined_imputed %>%
  lapply(., function(x) x %>% filter(!collection_dateerror))

prob_seqs <- lapply(metadata_joined_imputed, function(x) x %>% 
                      filter(collection_dateerror) %>% 
                      pull(tipnames))

renamed_alignments_dateschecked <- renamed_alignments %>%
  mapply(function(alignment,seqnames_to_drop) alignment[!rownames(alignment) %in% seqnames_to_drop,], 
         .,
         prob_seqs, 
         SIMPLIFY = FALSE)

################################### OUTPUT #####################################
# Save output files, plots, or results
alignmentfiles_newnames <- alignmentfiles  %>%
  gsub('.*Re_H5_', '',. ) %>%
  gsub('_.*', '',. ) %>%
  paste0(gsub('.*_','',segnames), '_', gsub('_.*','',segnames),'_', ., '.fasta') %>%
  paste0('./2023Dec02/alignments/', .)

mapply(write.dna, 
       renamed_alignments_dateschecked, 
       alignmentfiles_newnames,
       format = 'fasta')


metadatafiles <- paste('./2023Dec02/metadata/', 
                       gsub('.*_','',segnames), '_', gsub('_.*','',segnames),
                       '.csv',  
                       sep = '')

mapply(write_csv, 
       quote= 'needed',
       metadata_joined_imputed_dateschecked, 
       metadatafiles)  
#################################### END #######################################
################################################################################