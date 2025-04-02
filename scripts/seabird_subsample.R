####################################################################################################
####################################################################################################
## Script name:
##
## Purpose of script:
##
## Date created: 2025-04-02
##
##
########################################## SYSTEM OPTIONS ##########################################
options(scipen = 6, digits = 7) 
memory.limit(30000000) 

  
########################################## DEPENDENCIES ############################################
# Packages
library(tidyverse)
library(magrittr)
library(ape)

# User functions

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


seabirds <- read_csv('~/Downloads/seabird_H5Nx_meta.csv') %>% 
  mutate(is_seabird = 1) %>%
  select(tipnames, is_seabird) 

alignment_files <- c(list.files('./2024Aug18/region_alignments',
                                pattern = 'ha_',
                                full.names = T), 
                     list.files('./2024Aug18/reassortant_alignments',
                                pattern = 'ha_',
                                full.names = T))


alignments <- lapply(alignment_files, read.dna, format = 'fasta', as.matrix = F)

############################################## MAIN ################################################
temp <-  do.call(c, alignments) 
names(temp) <- gsub('.*\\.', '', names(temp))
temp <- temp[!duplicated(names(temp))]


groups <- GroupSequences(as.matrix(temp), snp_threshold = 4)

groups_with_seabirds <- groups %>% left_join(seabirds) %>%
  group_by(sequence_group) %>%
  replace_na(list(is_seabird = 0)) %>%
  filter(sum(is_seabird)>0) %>%
  ungroup()


ha_seabird_subsample <- all_meta %>% 
  drop_na(collection_date) %>% 
  filter(grepl('H5', virus_subtype)) %>% 
  filter(isolate_id %in% str_extract(names(temp), "EPI_ISL_(china_){0,1}\\d+[^.|]*")) %>%
  filter(clade  == '2344b') %>%
  left_join(groups_with_seabirds) %>%
  drop_na(sequence_group) %>%
  drop_na(cluster_profile) %>%
  
  mutate(new_host = case_when(is_seabird == 1 ~ 'seabird',
                              .default = host_simplifiedhost) ) %>%
  
  mutate(collection_regionname = case_when(grepl('europe', collection_regionname) ~ 'europe',
                                           grepl('africa', collection_regionname) ~ 'africa',
                                           grepl('asia', collection_regionname) ~ 'asia',
                                           grepl('(central|northern) america|caribbean', collection_regionname) ~ 'central & northern america',
                                           grepl('south america|southern ocean', collection_regionname) ~ 'south america',
                                           grepl('australia|melanesia', collection_regionname) ~ 'australasia',
                                           .default = collection_regionname)) %>%
  
  group_by(sequence_group, collection_regionname, collection_datemonth, host_simplifiedhost) %>% 
  slice_sample(n = 1) 

ha_seabird_subsample %<>% 
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
ha_selected_genomes <- temp[str_extract(names(temp), "EPI_ISL_(china_){0,1}\\d+[^.|]*") %in% ha_seabird_subsample$isolate_id]
ha_selected_genomes <- ha_selected_genomes[!duplicated(names(ha_selected_genomes))]

# rename  
renamed_ha_alignment <- ReNameAlignment2(ha_selected_genomes, ha_seabird_subsample)





############################################## WRITE ###############################################
write.FASTA(renamed_ha_alignment, './2025Apr2/ha_seabird_subsample.fasta')


write_delim(ha_seabird_subsample %>% ungroup()%>% select(new_tipnames, lat, long, collection_regionname, new_host, cluster_profile) %>% rename(tipnames = new_tipnames) %>% mutate(collection_regionname = gsub(' \\& | ', '_', collection_regionname)),
            delim = '\t',
            quote= 'needed',
            './2025Apr2/ha_seabird_subsample_traits.txt')




############################################## END #################################################
####################################################################################################
####################################################################################################