####################################################################################################
####################################################################################################
## Script name: Create alignments and metadata for Global subsample
##
## Purpose of script: Create alignments and metadata for Global subsample. HA and PB2 Only (Matched).
## n ~ 800
##
##
## Date created: 2024-12-23
##
##
########################################## SYSTEM OPTIONS ##########################################
#options(scipen = 6, digits = 7) 
memory.limit(30000000) 


########################################## DEPENDENCIES ############################################
# Packages
library(tidyverse)
library(magrittr)


# User functions
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


ReNameAlignment2 <- function(alignment, data){
  z <- names(alignment)
  
  isolates <-  regmatches(z, gregexpr("EPI_ISL_(china_){0,1}\\d+[^.|]*", z)) %>% unlist()
  
  new_seqnames <- sapply(isolates, function(x) data$new_tipnames[data$isolate_id %in% x]) %>% 
    as.vector() 
  
  names(alignment) <-  new_seqnames
  
  return(alignment)
}


############################################## DATA ################################################
all_meta <- read_csv('./2024-09-09_meta.csv')

alignment_files <- c(list.files('./2024Aug18/region_alignments',
                                pattern = 'ha_',
                                full.names = T), 
                     list.files('./2024Aug18/reassortant_alignments',
                                pattern = 'ha_',
                                full.names = T))


alignments <- lapply(alignment_files, read.dna, format = 'fasta', as.matrix = F)
names(alignments) <-  gsub('.*\\/h5_|.fasta$', '', alignment_files)


pb2_alignment_files <- c(list.files('./2024Aug18/region_alignments',
                                pattern = 'pb2_',
                                full.names = T), 
                     list.files('./2024Aug18/reassortant_alignments',
                                pattern = 'pb2_',
                                full.names = T))


pb2_alignments <- lapply(pb2_alignment_files, read.dna, format = 'fasta', as.matrix = F)
names(pb2_alignments) <-  gsub('.*\\/h5_|.fasta$', '', pb2_alignment_files)


# concatenate alignments (within segment)
temp <-  do.call(c, alignments) 
names(temp) <- gsub('.*\\.', '', names(temp))
temp <- temp[!duplicated(names(temp))]


pb2_temp <-  do.call(c, pb2_alignments) 
names(pb2_temp ) <- gsub('.*\\.', '', names(pb2_temp ))
pb2_temp  <- pb2_temp [!duplicated(names(pb2_temp ))]

############################################## MAIN ################################################

ha_global_subsample <- all_meta %>% 
  drop_na(collection_date) %>% 
  filter(grepl('H5', virus_subtype)) %>% 
  filter(isolate_id %in% str_extract(names(temp), "EPI_ISL_(china_){0,1}\\d+[^.|]*")) %>%
  filter(clade  == '2344b') %>%
  drop_na(cluster_profile) %>%
  group_by(collection_datemonth, collection_regionname, cluster_profile) %>% 
  slice_sample(n = 1) 

ha_global_subsample %<>% 
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

# select only genomes that are included 
ha_selected_genomes <- temp[str_extract(names(temp), "EPI_ISL_(china_){0,1}\\d+[^.|]*") %in% ha_global_subsample$isolate_id]
ha_selected_genomes <- ha_selected_genomes[!duplicated(names(ha_selected_genomes))]

# rename  
renamed_ha_alignment <- ReNameAlignment2(ha_selected_genomes, ha_global_subsample)



# PB2
pb2_global_subsample <- all_meta %>% 
  drop_na(collection_date) %>% 
  #filter(grepl('H5', virus_subtype)) %>% 
  filter(isolate_id %in% str_extract(names(pb2_temp), "EPI_ISL_(china_){0,1}\\d+[^.|]*")) %>%
  filter(clade  == '2344b') %>%
  group_by(collection_datemonth,  collection_regionname, cluster_profile) %>% 
  slice_sample(n = 1)


pb2_global_subsample %<>% 
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
  mutate(new_tipnames = gsub(' ', '_', new_tipnames))  %>%
  mutate(lat = coalesce(collection_subdiv1lat, 
                        collection_countrylat)) %>%
  mutate(long = coalesce(collection_subdiv1long, 
                         collection_countrylong)) %>%
  mutate(across(c(lat, long), .fns = ~ as.numeric(.x)))

# select only genomes that are included 
pb2_selected_genomes <- pb2_temp[str_extract(names(pb2_temp), "EPI_ISL_(china_){0,1}\\d+[^.|]*") %in% pb2_global_subsample$isolate_id]
pb2_selected_genomes <- pb2_selected_genomes[!duplicated(names(pb2_selected_genomes))]

# rename  
renamed_pb2_alignment <- ReNameAlignment2(pb2_selected_genomes, pb2_global_subsample)





############################################## WRITE ################################################
write.FASTA(renamed_ha_alignment, './2025Feb26/globalsubsample/ha_global_subsample.fasta')
write.FASTA(renamed_pb2_alignment, './2025Feb26/globalsubsample/pb2_global_subsample.fasta')


write_delim(pb2_global_subsample %>% ungroup()%>% select(new_tipnames, lat, long) %>% rename(tipnames = new_tipnames),
            delim = '\t',
            quote= 'needed',
            './2025Feb26/globalsubsample/pb2_global_subsample_traits.txt')

write_delim(ha_global_subsample %>% ungroup() %>% select(new_tipnames, lat, long) %>% rename(tipnames = new_tipnames),
            delim = '\t',
            quote= 'needed',
            './2025Feb26/globalsubsample/ha_global_subsample_traits.txt')
############################################## END ################################################
####################################################################################################
####################################################################################################