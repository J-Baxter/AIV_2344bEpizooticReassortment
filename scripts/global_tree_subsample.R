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

############################################## DATA ################################################
all_meta <- read_csv('./2024-09-09_meta.csv')

alignments <- c(list.files('./2024Aug18/reassortant_alignments',
                         pattern = 'ha'),
                grep(list.files('./2024Aug18/region_alignments',
                                pattern = 'ha_'), pattern='southamerica', invert=TRUE, value=TRUE),
                list.files('./2024Sept16/south_america',
                           pattern = 'ha_.*fasta'))

############################################## MAIN ################################################

ha_global_subsample <- all_meta %>% 
  drop_na(collection_date) %>% 
  filter(grepl('H5', virus_subtype)) %>% 
  group_by(collection_datemonth, collection_regionname) %>% 
  slice_sample(n = 2)

pb2_global_subsample <- all_meta %>% 
  drop_na(collection_date) %>% 
  #filter(grepl('H5', virus_subtype)) %>% 
  group_by(collection_datemonth, collection_regionname) %>% 
  slice_sample(n = 1)

# HA subsample: n = 948
# PB2 subsample: n = 839
############################################## WRITE ################################################


############################################## END ################################################
####################################################################################################
####################################################################################################