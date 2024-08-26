# Dependencies
library(tidyverse)
library(ape)
source('./scripts/FormatBirds.R')
source('./scripts/FormatMammals.R')

geodata <- read_csv('./annotated_geodata.csv')
birds <- read_csv('bird_taxonomy.csv')
mammals <- read_csv('mammal_taxonomy.csv')

load("./2024Jun01/h5nx_2344b_clusters_20240513.Rda")
meta <- meta %>% 
  as_tibble() %>%
  select(-date_frac) %>%
  mutate(date_year = as.double(date_year),
         location_1 = case_when(location_1 == 'Antarctica' ~ 'South America', 
                                .default = location_1))




all_meta <- read_csv('./2024Jul12/all_metadata.csv')

meta_key <- all_meta %>%
  select(-c(
    label,
    tipnames,
    cluster_number,
    collection_dateerror
  )) %>%
  distinct()

duplicated <- meta_key %>%
  group_by(isolate_id) %>%
  mutate(n = n()) %>%

  mutate(across(starts_with('host'), .fns = ~ case_when(n >1 & grepl(' ', host_class) ~ NA,
                                                        .default = .x))) %>%
  summarise(across(everything(), ~ na.omit(.x)[1])) %>%
  select(-n) %>%
  ungroup() 


length(duplicated$isolate_id) == length(unique(duplicated$isolate_id))

test <- duplicated %>%
  select(isolate_id,
             host_sciname) %>%
  rename(source = host_sciname)
  

  
test3 <- meta %>% 
  filter(host_class == 'Mammal') %>%
  select(c(starts_with('isolate'),
           starts_with('host'))) %>%
  mutate(source = case_when(str_count(isolate_name, '/') < 4 ~ NA, 
                            .default = str_extract(isolate_name, "(?<=/)[^/]+"))) %>%
  mutate(source =  str_to_lower(source) %>%
           str_replace_all("[[:punct:]]", ' ') %>%
           str_trim(.)) %>%
  rows_patch(., test, by = 'isolate_id', unmatched = 'ignore') %>%
  rowwise() %>%
  mutate(primary_com_name = FormatMammal(source)) %>%
  as_tibble() %>%
  left_join(mammals,  by = join_by(primary_com_name)) %>%
  select(-c(host_order,
            host_class,
            source)) %>%
  rename(host_class = class,
         host_order = order,
         host_family = family,
         host_sciname = sci_name,
         host_commonname = primary_com_name)


test4 <- duplicated %>%
  filter(host_class %in% c('unknown', 'mammal') | is.na(host_class)) %>%
  rows_update(., test3, by = 'isolate_id', unmatched = 'ignore') %>%
  filter(is.na(host_sciname)) %>%
  select(-starts_with('host')) %>%
  mutate(source = case_when(!grepl('unknown', isolate_name) ~ 'homo sapiens',
                            .default = 'unknown')) %>% 
  rowwise() %>%
  mutate(primary_com_name = FormatMammal(source)) %>%
  as_tibble() %>%
  left_join(mammals,  by = join_by(primary_com_name)) %>%
  rename(host_class = class,
         host_order = order,
         host_family = family,
         host_sciname = sci_name,
         host_commonname = primary_com_name) %>%
  select(-source)

test5 <- duplicated %>%
  rows_update(., test3, by = 'isolate_id') %>%
  rows_update(., test4, by = 'isolate_id')
  

test6 <- test5 %>% 
  filter(host_class %in% c("anseriformes domestic",
                           "avian other",
                           "anseriformes wild", 
                           "charadriiformes", 
                           "galliformes")) %>%
  select(c(starts_with('isolate'),
           starts_with('host'),
           collection_regionname)) %>%
  mutate(source = case_when(str_count(isolate_name, '/') < 4 ~ NA, 
                            .default = str_extract(isolate_name, "(?<=/)[^/]+"))) %>%
  mutate(source =  str_to_lower(source) %>%
           str_replace_all("[[:punct:]]", ' ') %>%
           str_trim(.)) %>%
  rowwise() %>%
  mutate(primary_com_name = FormatBirds(source)) %>%
  as_tibble() %>%
  mutate(primary_com_name = case_when(
    
    # location-specific birds
    grepl('sparrowhawk', 
          primary_com_name) & grepl('europe', collection_regionname) ~ 'eurasian sparrowhawk',
    grepl("^magpie$|^common magpie",
          primary_com_name) & grepl('europe', collection_regionname) ~ "eurasian magpie",
    grepl("^magpie$|^common magpie",
          primary_com_name) & grepl('eastern asia', collection_regionname) ~ "oriental/eurasian magpie",
    grepl("^magpie$|^common magpie", 
          primary_com_name) & grepl('northern america', collection_regionname) ~ "black-billed magpie",
    grepl("^kestrel$",
          primary_com_name) & grepl('northern america', collection_regionname) ~ "american kestrel",
    grepl("^kestrel$", 
          primary_com_name) & grepl('europe', collection_regionname) ~ "lesser/eurasian kestrel",
    grepl("^guinea {0,1}fowl", 
          primary_com_name) & grepl('africa', collection_regionname) ~ "crested guineafowl sp.",
    grepl('^guinea {0,1}fowl', 
          primary_com_name) & !grepl('africa', collection_regionname) ~ "helmeted guineafowl (domestic type)",
    grepl('^fulmar$', 
          primary_com_name) & !grepl('south america', collection_regionname) ~ "northern fulmar",
    grepl('^fulmar$',
          primary_com_name) & grepl('south america', collection_regionname) ~ "southern fulmar",
    grepl('^gannet$', 
          primary_com_name) & grepl('europe', collection_regionname) ~ "northern gannet",
    grepl('jungle crow', 
          primary_com_name) ~ "crow/raven sp.",
    grepl('oystercatcher', 
          primary_com_name) & grepl('europe', collection_regionname) ~ 'eurasian oystercatcher',
    grepl("turkey|^pavo$", 
          primary_com_name) & grepl('central america', collection_regionname) ~ "wild turkey", 
    grepl("turkey|^pavo$",
          primary_com_name) & !grepl('central america', collection_regionname) ~ "wild turkey (domestic type)", 
    grepl("turnstone", 
          primary_com_name) & grepl('europe', collection_regionname) ~ "ruddy turnstone",
    primary_com_name == 'sea eagle' & grepl('europe', collection_regionname) ~"white-tailed eagle",
    grepl("quail", 
          primary_com_name) & !grepl('america', collection_regionname) ~ "old world quail sp.", 
    grepl("quail", 
          primary_com_name) & grepl('america', collection_regionname) ~ "new world quail sp.", 
    grepl('^vulture$', 
          primary_com_name) &  grepl('america', collection_regionname) ~ "new world vulture sp.",
    grepl('^vulture$', 
          primary_com_name) & !grepl('america', collection_regionname) ~ "old world vulture sp.",
    
    # unknown
    grepl('^an$|unknown', 
          primary_com_name) ~ 'unknown',
    is.na(primary_com_name) ~ 'unknown',
    
    
    .default = primary_com_name
  )) %>%
  left_join(birds,  by = join_by(primary_com_name)) %>%
  select(-c(host_order,
            host_class,
            host_family,
            host_sciname,
            host_class,
            host_commonname,
            host_isbird,
            host_isdomestic,
            host_simplifiedhost)) %>%
  rename(host_class = class,
         host_order = order,
         host_family = family,
         host_sciname = sci_name,
         host_commonname = primary_com_name) %>%
  # environment
  mutate(across(c(host_order, host_class, host_family),
                ~ case_when(host_commonname == 'environment' ~ 'environment', TRUE ~ .)),
         across(c(host_order, host_class, host_family),
                ~ case_when(host_commonname == 'unknown' ~ 'unknown', TRUE ~ .)))  %>%
  select(-source) 
  
  


test_final <- duplicated %>%
  rows_update(., test3, by = 'isolate_id') %>%
  rows_update(., test4, by = 'isolate_id') %>%
  rows_update(., test6, by = 'isolate_id') %>%
  # binary host information
  mutate(host_isdomestic = case_when(grepl("domestic|homo sapiens", host_sciname) ~ "domestic",
                                     .default = "wild"
  )) %>%
  #rowwise() %>%
  mutate(host_isbird = case_when(host_class == 'aves' ~ TRUE,
                                 .default = FALSE))%>% 
 # as_tibble() %>%

  mutate(host_class = str_to_lower(host_class)) %>%
  mutate(host_simplifiedhost = case_when(
    host_order %in% c('anseriformes', 'galliformes', 'charadriiformes') ~ paste(host_order, 
                                                                                host_isdomestic,
                                                                                sep = '-'),
    host_class == 'mammalia' & host_sciname != 'homo sapiens' ~ 'mammal',
    host_class == 'mammalia' & host_sciname == 'homo sapiens' ~ 'human',
    host_order == 'environment' | host_commonname == 'environment' ~ 'environment',
    host_class == 'aves' & !(host_order %in% c('anseriformes', 'galliformes', 'charadriiformes')) ~ 'other-bird',
    .default = 'unknown'
  ))  %>%
  
  #make tipnames
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
  
  
write_csv(test_final, '2024-08-19_meta.csv')


# Rename all alignments
ReNameAlignment <- function(alignment, data){
  z <- rownames(alignment)
  
  isolates <-  str_extract(rownames(alignment), "EPI_ISL_(china_){0,1}\\d+[^.|]*")
  
  new_seqnames <- sapply(isolates, function(x) data$tipnames[data$isolate_id %in% x]) %>% 
    as.vector() 
  
  stopifnot(length(new_seqnames) == nrow(alignment))
  
  rownames(alignment) <-  new_seqnames
  
  return(alignment)
}

all_meta <- read_csv('2024-08-19_meta.csv')

# Reassortants
reassortant_files <- list.files('./2024Aug18/reassortant_alignments',
                                full.names = T)

reassortant_aln <- lapply(reassortant_files,
                          read.dna, 
                          format = 'fasta',
                          as.matrix = TRUE) %>%
  lapply(., function(x) x[!duplicated(rownames(x)),])

reassortant_aln_renamed <- reassortant_aln %>%
  lapply(., ReNameAlignment, all_meta)
  
  

mapply(write.dna, 
       reassortant_aln_renamed, 
       reassortant_files,
       format = 'fasta')


# Regions
region_files <- list.files('./2024Aug18/region_alignments',
                           full.names = T)

region_aln <- lapply(region_files,
                          read.dna, 
                          format = 'fasta',
                          as.matrix = TRUE) %>%
  lapply(., function(x) x[!duplicated(rownames(x)),])

region_aln_renamed <- region_aln %>%
  lapply(., ReNameAlignment, all_meta)

mapply(write.dna, 
       region_aln_renamed, 
       region_files,
       format = 'fasta')