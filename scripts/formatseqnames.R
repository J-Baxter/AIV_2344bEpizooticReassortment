####################################################################################################
# Reformat initial alignments and phylogenies


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
source("./scripts/Subsamplefunctions.R")
#sourceCpp('./scripts/getLocation.cpp') #~3mins


# Join alignment metadata with reassortant data
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
# Import tree data
#treefiles <- list.files(path = './2023/phylo_ml/original_mlphylos',
                      #  recursive = FALSE,
                       # include.dirs = FALSE, 
                        #full.names = TRUE)

#trees <- lapply(treefiles, read.newick) 


# Subset alignments
alignmentfiles <- list.files('./2023Dec01/alignments',
                             full.names = T, 
                             pattern="\\.")

alignments <- lapply(alignmentfiles,
                     read.dna, 
                     format = 'fasta',
                     as.matrix = T)


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
#names(trees) <- segnames

# Import existing metadata
reassortant_metadata <- read_csv('./2023Dec01/metadata/h5_metadata_global_6280_update.csv') %>%
  mutate(date = dmy(date) %>%
           as.Date())


####################################################################################################
# Extract metadata and format tipnames
seqnames <- lapply(alignments, rownames) %>% 
  setNames(segnames) # Only required because we are working backwards from tipnames


# Extract metadata from tipnames
metadata_unformatted <- lapply(seqnames, ExtractMetadata)

# Format metadata
reassortant_metadata_formatted <- FormatMetadata(reassortant_metadata) %>%
  select(-matches('location.[cn]')) %>%
  select(-c(domestic.status, 
            domestic.wild)) %>% 
  unite('collection.original', 
        starts_with('location'), 
        sep = ' / ')

metadata_formatted <- lapply(metadata_unformatted, FormatMetadata) 



metadata_joined <- lapply(metadata_formatted, MyFunc, reassortant_metadata_formatted ) %>%
  setNames(segnames) %>%
  lapply(., MakeTipNames)


metadata_joined_df <- metadata_joined %>% 
  bind_rows(., .id = 'virus.segment') 

######### Clade and cluster guess here ##########

####################################################################################################
# Rename alignments with new-format tipnames
renamed_alignments <- alignments %>%
  mapply(ReNameAlignment, 
         .,
         metadata_joined , 
         SIMPLIFY = FALSE)

stop()