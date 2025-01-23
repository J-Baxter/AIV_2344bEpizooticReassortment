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
library(sf)

# User functions
source('./scripts/FormatBirds.R')

############################################## DATA ################################################

cag_data <- read_csv('./new_data/cag_data_2023Dec23.csv')%>%
  mutate(collection_year= year(collection_date) %>% as.character())

cag_alns <- list.files('./new_data/',
                       pattern = '.fas',
                       full.names = T) %>%
  lapply(., read.FASTA)


# external data (shapefile and bird taxonomy)
china_level2 <- read_sf('./new_data/gadm41_CHN_shp/gadm41_CHN_2.shp')
birds <- read_csv('bird_taxonomy.csv')


############################################## MAIN ################################################

# Format Metadata and Name sequences
# Vars to include: date, lat-long + region, altitude (continuous + binned), host species
sequence_name_key <- cag_alns %>%
  lapply(., names) %>%
  unlist() %>%
  as_tibble_col(column_name = 'original_name') %>%
  mutate(temp =   str_replace_all(original_name, '_$|\\)|_(?=Lake)', '') %>% # remove trailing _
           str_replace_all(., '/|__| |\\(', '_'))

sequence_name_data <- sequence_name_key %>%
  select(temp) %>%
  distinct() %>%
  separate_wider_delim(temp, delim = '_', too_few = 'align_end', names_sep = '_',cols_remove = F) %>%
  mutate(across(starts_with('temp'), .fn = ~ gsub('^A$', NA_character_, .x))) %>%
  unite(host, temp_1, temp_2, temp_3, temp_4, na.rm = T, sep = ' ') %>%
  mutate(temp_7 = gsub('-.*', '', temp_7),
         temp_5 = gsub('QinghaiLake', 'Qinghai Lake', temp_5) %>%
           gsub('[tT]umuji|dongpaozi north', 'Tumuji Nature Reserve', .) %>% # Note this is inferred
           gsub('[bB]ird [iI]sland', 'Niaodao Scenic Spot',.)%>% # Note this is inferred
           gsub('[Hh]ada [Bb]each', 'Shadao',.)) %>% # Note this is inferred
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
  
  mutate(location_query =   gsub('*[tT]umuji.*|*dongpaozi north.*', 'Tumuji Nature Reserve, China', location_query) %>% # Note this is inferred
           gsub('[bB]ird [iI]sland', 'Niaodao Scenic Spot',.)%>% # Note this is inferred
           gsub('Sanshui Factory', '',.)%>% # Note this is an approximation. I (JB) cannot find the exact location
           gsub(p,.)) %>%  # Note this is inferred
  as_tibble() %>%

    # geocode to obtain lat-lon coordinates
  mutate(location = geocode(location_query)) %>%  
  unnest(location) %>%
  st_as_sf(coords = c('lon', 'lat'), crs = 4326) %>%
  
  
  
  # obtain elevations for inferred coordinates
  rowwise() %>%
  get_elev_point(prj = 4326, src = "aws") %>%
  #select(-elev_units) %>%
  #as_tibble()
  #st_as_sf(coords = c('lon', 'lat'), crs = 4326) %>% 
  st_join(.,
        china_level2 , 
        join = st_within, 
        # left = FALSE,
        largest = TRUE) %>%
  
  # Set deprecated levels as NA (ie where we only have the province level location to start)
  mutate(across(ends_with('_2'), .fns = ~ case_when(collection_region == NAME_1 & is.na(collection_location) ~ NA_character_, 
                                                    .default = .x))) %>% 
  as_tibble()
  
  



  # group bird species into orders
temp_3 <- temp %>%
  rowwise() %>%
  mutate(primary_com_name = FormatBirds(tolower(host))) %>%
  left_join(birds) %>%
  
  # cut elevation into classes
  mutate(elevation_class = cut(elevation,
                             breaks = c(0, 500, 2000, 3000, 5500, Inf),
                             labels = c('sea-level', 'low', 'medium', 'high', 'extreme')))
  
  
  
# drop deprecated columns and tidy
out <- temp_3 %>%
  dplyr::select(sequence_isolate,
                sequence_id,
                virus_subtype,
                collection_date,
                collection_year,
                primary_com_name,
                sci_name,
                order,
                family,
                class,
                COUNTRY,
                NAME_1,
                NL_NAME_1,
                NAME_2,
                NL_NAME_2,
                geometry,
                elevation,
                elevation_class,
                sequence_name) %>%
  rename(
         host_commonname = primary_com_name,
         host_sciname = sci_name,
         host_order = order,
         host_family = family,
         host_class =  class,
         collection_country = COUNTRY,
         collection_subdiv1 = NAME_1,
         collection_subdiv1_hanzi = NL_NAME_1,
         collection_subdiv2 = NAME_2,
         collection_subdiv2_hanzi = NL_NAME_2) %>%
  
  mutate(unique_id = coalesce(sequence_isolate, sequence_id),
         tipdate = coalesce(as.character(collection_date), collection_year)) %>%
  
  mutate(tipnames = paste(virus_subtype, 
                           '2344b',
                          unique_id,
                           host_order,
                           tolower(collection_country),
                           NA,
                          tipdate,
                           sep = '|')) %>%
  select(-c(unique_id, tipdate)) %>%
  rename(old_tipnames = sequence_name) %>%
  relocate(tipnames)
  
# make new tipnames

data_with_old_names <- sequence_name_key %>%
  left_join(out, by = join_by(temp == old_tipnames)) %>%
  distinct()

named_alignments <- lapply(cag_alns,  function(x) setNames(x, data_with_old_names %>% 
                                                             filter(original_name %in% names(x)) %>%
                                                             pull(tipnames)))


############################################## WRITE ################################################
cag_data <- read_csv('./new_data/cag_data_2023Dec23.csv')%>%
  mutate(collection_year= year(collection_date) %>% as.character())

out_filenames <- list.files('./new_data',
                       pattern = '.fas',
                       full.names = T) %>%
  gsub('\\.fas$', '_formatted.fasta',.)

mapply(write.FASTA,
       named_alignments,
       out_filenames )

write_csv(out, './new_data/cag_data_2025Jan06.csv')
############################################## END ################################################
####################################################################################################
####################################################################################################