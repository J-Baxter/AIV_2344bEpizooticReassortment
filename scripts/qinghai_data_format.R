####################################################################################################
####################################################################################################
## Script name: Create alignments and metadata for Qinghai (altitude) subsample
##
## Purpose of script: Create alignments and metadata for Qinghai (altitude) subsample, combining
## new data with existing HA and PB2 Only (Matched). NB Temporary script which will need to be added
## to new project directory.
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
library(ggmap)
library(elevatr)
library(ape)
library(lubridate)
library(auk)

# User functions

############################################## DATA ################################################

cag_data <- read_csv('./new_data/cag_data_2023Dec23.csv')%>%
  mutate(collection_year= year(collection_date) %>% as.character())

cag_alns <- list.files('./new_data/',
                       pattern = '.fas',
                       full.names = T) %>%
  lapply(., read.FASTA)
  
############################################## MAIN ################################################

# Format Metadata and Name sequences
# Vars to include: date, lat-long + region, altitude (continuous + binned), host species
sequence_name_data <- cag_alns %>%
  lapply(., names) %>%
  unlist() %>%
  str_replace_all(., '_$|\\)|_(?=Lake)', '') %>% # remove trailing _
  str_replace_all(., '/|__| |\\(', '_') %>%
  as_tibble_col(column_name = 'temp') %>%
  distinct() %>%
  separate_wider_delim(temp, delim = '_', too_few = 'align_end', names_sep = '_',cols_remove = F) %>%
  mutate(across(starts_with('temp'), .fn = ~ gsub('^A$', NA_character_, .x))) %>%
  unite(host, temp_1, temp_2, temp_3, temp_4, na.rm = T, sep = ' ') %>%
  mutate(temp_7 = gsub('-.*', '', temp_7),
         temp_5 = gsub('QinghaiLake', 'Qinghai Lake', temp_5)) %>%
  rename(sequence_name = temp_temp,
         collection_region =temp_5,
         sequence_id = temp_6,
         collection_year = temp_7,
         virus_subtype = temp_8)


temp <- cag_data %>% 
  full_join(sequence_name_data,
            by = join_by(sequence_id,
                         virus_subtype, 
                         collection_year)) %>%
  
  # join collection location to formulate search query]
  unite(location_query, 
        collection_location, 
        collection_region,
        sep = '  ', 
        na.rm = TRUE, 
        remove = FALSE) %>%
  rowwise() %>%
  mutate(location_query = c(str_split(location_query,  ' ')) %>%
           unlist() %>%
           str_unique() %>%
           str_c(collapse= ' ') %>%
           str_replace_all('[:punct:]', '') %>%
           str_to_lower() %>%
           str_replace_all('  ', ' ') %>%
           str_c(' china')) %>%
  as_tibble() %>%

    # geocode to obtain lat-lon coordinates
  mutate(location = geocode(location_query)) %>%  
  unnest(location) %>%
  st_as_sf(coords = c('lon', 'lat'), crs = 4326) %>%
  
  
  
  # obtain elevations for inferred coordinates
  rowwise() %>%
  get_elev_point(prj = 4326, src = "aws") %>%
  select(-elev_units) %>%
  as_tibble()

  # group bird species into orders
temp_3 <- temp_2 %>%
  
  
  # format dates 
  




# Filter those for which we have all segments available (ie full genomes)


# 

############################################## WRITE ################################################


############################################## END ################################################
####################################################################################################
####################################################################################################