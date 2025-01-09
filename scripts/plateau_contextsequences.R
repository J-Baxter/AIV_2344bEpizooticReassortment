####################################################################################################
####################################################################################################
## Script name: Create initial context sequences for Qinghai data
##
## Purpose of script:
##
## Date created: 2025-01-08
##
##
########################################## SYSTEM OPTIONS ##########################################
options(scipen = 6, digits = 7) 
memory.limit(30000000) 
source('./scripts/FindIdenticalSeqs.R')
  
########################################## DEPENDENCIES ############################################
# Packages
library(tidyverse)
library(magrittr)
library(ape)
library(igraph)
library(phytools)
library(phangorn)

# User functions

# From ./scripts/july_subsampling.R
GroupSequences <- function(aln, snp_threshold = 0){
  require(igraph)
  
  # Ensure alignment is correctly formatted
  if(class(aln) != 'PhyDat'){
    aln_formatted <- as.phyDat(aln) 
  }else{
    aln_formatted <- aln
  }
  
  # Calculate hamming distance
  hd_normalised <- dist.hamming(aln_formatted) %>%
    as.matrix()
  hd_raw <- hd_normalised * ncol(aln)
  
  # Obtain groups of sequences for which HD < SNP threshold
  
  if( any(hd_raw[lower.tri(hd_raw, diag = FALSE)] <= snp_threshold)){
    groups <- which(hd_raw <= snp_threshold, 
                    arr.ind = TRUE) %>%
      dplyr::as_tibble(rownames = 'tipnames') %>% 
      filter(row !=col) %>%
      dplyr::select(-tipnames) %>%
      
      # Infer network from HDs
      igraph::graph_from_data_frame(.,
                                    directed = F) %>%
      
      components() %>%
      getElement('membership') %>%
      stack() %>%
      as_tibble() %>%
      mutate(ind = as.numeric(as.character(ind))) %>%
      mutate(tipnames = map_chr(ind, ~ rownames(aln)[.x])) %>%
      dplyr::select(c(tipnames, values)) %>%
      dplyr::distinct() %>%
      dplyr::rename(sequence_group = values)
  }else{
    print('all sequences above threshold.')
    groups <- tibble(tipnames = rownames(aln)) %>%
      rowid_to_column(., var = 'sequence_group')
  }
  
  
  out <- tibble(tipnames = rownames(aln)) %>%
    left_join(groups) %>%
    mutate(sequence_group = 
             ifelse(is.na(sequence_group), 
                    max(sequence_group, na.rm = T) + row_number() + 1, 
                    sequence_group))
  
  
  return(out)
}



ReNameAlignment2 <- function(alignment, data){
  z <- names(alignment)
  
  isolates <-  regmatches(z, gregexpr("EPI_ISL_(china_){0,1}\\d+[^.|]*", z)) %>% unlist()
  
  new_seqnames <- sapply(isolates, function(x) data %>% filter(isolate_id %in% x) %>% pull(new_tipnames)) %>% 
    as.vector() 
  
  names(alignment) <-  new_seqnames
  
  return(alignment)
}

#data$new_tipnames[data$isolate_id %in% x]


############################################## DATA ################################################

# Import plateau alignments
qinghai_ha <- read.dna('./new_data/2024Dec23_CAG_h5nX_ha_formatted_aligned.fasta',
                       format = 'fasta',
                       as.matrix = T)

qinghai_pb2 <- read.dna('./new_data/2024Dec23_CAG_h5nX_pb2_formatted_aligned.fasta',
                       format = 'fasta',
                       as.matrix = T)



# Import other alignments
# HA
ha_alignmentfiles <- c(list.files('./2024Aug18/region_alignments',
                                pattern = 'ha_',
                                full.names = T), 
                     list.files('./2024Aug18/reassortant_alignments',
                                pattern = 'ha_',
                                full.names = T))


ha_alignments <- lapply(ha_alignmentfiles, read.dna, format = 'fasta', as.matrix = F)
names(ha_alignments) <-  gsub('.*\\/h5_|.fasta$', '', ha_alignmentfiles)


# PB2
pb2_alignmentfiles <- c(list.files('./2024Aug18/region_alignments',
                                    pattern = 'pb2_',
                                    full.names = T), 
                         list.files('./2024Aug18/reassortant_alignments',
                                    pattern = 'pb2_',
                                    full.names = T))


pb2_alignments <- lapply(pb2_alignmentfiles, read.dna, format = 'fasta', as.matrix = F)
names(pb2_alignments) <-  gsub('.*\\/h5_|.fasta$', '', pb2_alignmentfiles)


# concatenate alignments (within segment)
ha_temp <-  do.call(c, ha_alignments) 
names(ha_temp) <- gsub('.*\\.', '', names(ha_temp))
ha_alignment <- ha_temp[!duplicated(names(ha_temp))]


pb2_temp <-  do.call(c, pb2_alignments) 
names(pb2_temp) <- gsub('.*\\.', '', names(pb2_temp))
pb2_alignment  <- pb2_temp[!duplicated(names(pb2_temp))]


# Import all meta
all_meta <- read_csv('./2024-09-09_meta.csv')
qinghai_meta <- read_csv('./new_data/cag_data_2025Jan06.csv')


all_meta %<>% 
  #drop_na(collection_date) %>% 
  #filter(clade  == '2344b'| is.na(clade) ) %>%
  unite(new_tipnames, 
        virus_subtype,
        clade,
        isolate_id,
        host_order,
        collection_countryname,
        cluster_profile,
        collection_tipdate,
        sep = '|',                        
        remove = FALSE) %>%
  mutate(new_tipnames = gsub(' ', '_', new_tipnames))   %>%
  mutate(lat = coalesce(collection_subdiv1lat, 
                        collection_countrylat)) %>%
  mutate(long = coalesce(collection_subdiv1long, 
                         collection_countrylong)) %>%
  mutate(across(c(lat, long), .fns = ~ as.numeric(.x)))

# rename  
ha_alignment <- ReNameAlignment2(ha_alignment, all_meta) %>% as.matrix()
pb2_alignment <- ReNameAlignment2(pb2_alignment, all_meta) %>% as.matrix()


############################################## MAIN ################################################
ha_groups <-GroupSequences(ha_alignment, snp_threshold = 3) %>%
  mutate(isolate_id = str_extract(tipnames, "EPI_ISL_(china_){0,1}\\d+[^.|]*"))

pb2_groups <-  GroupSequences(pb2_alignment, snp_threshold = 8) %>%
  mutate(isolate_id = str_extract(tipnames, "EPI_ISL_(china_){0,1}\\d+[^.|]*"))


pb2_generalcontext <- pb2_groups %>% 
  mutate(sequence_group = as.factor(sequence_group)) %>%
  group_by(sequence_group) %>%
  slice_sample(n = 1) %>%
  ungroup()

ha_generalcontext <-ha_groups %>% 
  mutate(sequence_group = as.factor(sequence_group)) %>%
  group_by(sequence_group) %>%
  slice_sample(n = 1) %>%
  ungroup()

pb2_generalcontext_aln <- pb2_alignment[rownames(pb2_alignment) %in% pb2_generalcontext$tipnames,]
ha_generalcontext_aln <- ha_alignment[rownames(ha_alignment) %in% ha_generalcontext$tipnames,]


pb2_generalcontextwithqinghai_aln <- rbind(pb2_generalcontext_aln, qinghai_pb2) 
ha_generalcontextwithqinghai_aln <- rbind(ha_generalcontext_aln, qinghai_ha)
############################################## WRITE ###############################################

write.FASTA(pb2_generalcontextwithqinghai_aln, 
            './new_data/2025Jan09_qinghaigeneral_pb2.fasta')

write.FASTA(ha_generalcontextwithqinghai_aln, 
            './new_data/2025Jan09_qinghaigeneral_ha.fasta')

############################################## END #################################################
####################################################################################################
####################################################################################################