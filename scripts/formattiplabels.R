library(ape)
library(tidyverse)
library(treeio)
library(TreeTools)

pb2_tree <- read.newick('./data/phylo_ml/RAxML_bipartitions.H5_large_pb2')
ggtree(pb2_tree)

metadata <- read_csv('./data/metadata/h5_metadata.csv')
will_metadata <- read_csv('./data/metadata/Re_H5_clade2344b_Asia_blast_PB2_meta.will.2.csv')

t <- left_join(metadata, will_metadata )
tip_labels <- TipLabels(pb2_tree) 


joined_data <- tip_labels %>%
  sapply(., str_split_1, pattern = '\\|') %>%
  sapply(., `[`, 4) %>%
  bind_cols() %>%
  setNames('isolate_id') %>%
  left_join(metadata)

prob_seqs <- joined_data %>% filter(is.na(isolate_name))

tree <- pb2_tree
n<-length(tree$tip.label)
ee<-setNames(tree$edge.length[sapply(4:n,function(x,y)   which(y==x),y=tree$edge[,2])],tree$tip.label)


data_from_tipnames <- tip_labels %>% 
  gsub('\\|$', '', .) %>%
  sapply(., str_split_1, pattern = '\\|') %>%
  sapply(., function(x) ifelse(!grepl('\\.', tail(x, n = 4)), return(x), NA)) %>%
  .[!is.na(.)] %>%
  
  do.call(rbind.data.frame,. )


#Isolate id = starts with EPI_ISL
#subtype = HxNx
#isolate name = 
#segment
#date
birds <- read_csv('ebird_taxonomy_v2023.csv') %>% 
  rename_with( .fn = ~tolower(.x)) %>%
  dplyr::select(-c(contains('taxon'), species_code, report_as)) %>%
  filter(category %in% c('species', 'slash', 'domestic', 'spuh')) %>%
  mutate(across(contains('name'), .fns = ~ tolower(.x))) %>%
  mutate(primary_com_name = gsub('gray', 'grey', primary_com_name))


remainder <- 
MyFunc4 <- function(x){
  
  isolate.id = x[grep('EPI_ISL_*' ,x)]
  subtype = x[grep('H[[:digit:]]{1,}N[[:digit:]]*' ,x)]
  segment = x[grep('HA|PB4|PB2|PA|NP|NA|NS' ,x)]
  date = x[grep('^[[:digit:]][[:digit:]][[:digit:]][[:digit:]]-[[:digit:]][[:digit:]]-[[:digit:]][[:digit:]]' ,x)]
  isolate.name = x[grep('/[^/]*/[^/]*/' ,x)]
  
  if(identical(isolate.id, character(0))){
    isolate.id <- NA
  }
  
  if(identical(subtype, character(0))){
    subtype <- NA
    print(subtype)
  }
  if(identical(segment, character(0))){
    segment <- NA
  } 
  if(identical(date, character(0))){
    date <- NA
  } 
  if(identical(isolate.name, character(0))){
    isolate.name <- NA
  }
  
  out <- cbind.data.frame('isolate.id' = isolate.id,
                    'subtype' = subtype,
                    'segment' = segment,
                    'date' = date,
                    'isolate.name' = isolate.name) #%>%
  # as_tibble()
  
  
  return(out)
}

df <- lapply(data_from_tipnames, MyFunc4) %>%
  bind_rows() %>%
  mutate(subtype = str_extract(subtype, 'H[[:digit:]]{1,}N[[:digit:]]*')) %>%
  mutate(date = parsedate::parse_date(date) %>% 
           as.Date()) %>%
  mutate(decimal.date = decimal_date(date)) %>%
  mutate(isolate.name = gsub("_(\\d{4})|(\\d{4})_", "\\1",  perl = TRUE, isolate.name)) %>%
  
  mutate(isolate.name = case_when(!grepl( "/(\\d{4})$", isolate.name) & !grepl('/$', isolate.name) ~ paste0(isolate.name, '/NA'),
                                  !grepl( "/(\\d{4})$", isolate.name) & grepl('/$', isolate.name) ~ paste0(isolate.name, 'NA'),
                                  .default = isolate.name)) %>%
 
  separate_wider_delim(isolate.name,
                       delim = '/', 
                       names = c('virus_species', 'source', 'location', 'id_unsure', 'year'),
                       cols_remove = F, 
                       too_few = 'align_end',
                       too_many = 'debug') %>%
  dplyr::select(-c(year, contains('isolate.name_'))) %>%
  mutate(virus_species = case_when(str_length(source) == 1  ~ source,
                                   .default = virus_species)) %>%
  mutate(source = case_when(str_length(source) == 1 ~ NA,
                                   .default = source)) %>%

  
  
  # Bird names/order
  mutate(source =  gsub('_', ' ',tolower(source))) %>%
  mutate(primary_com_name = case_when(grepl('common|eurasian teal|green-winged-teal', source) ~ 'green-winged teal',
                                      grepl('spot-billed duck', source) ~ 'eastern spot-billed duck',
                                      grepl('crested grebe', source) ~ 'great crested grebe',
                                      grepl('eastern curlew', source) ~ 'far eastern curlew',
                                      grepl('bean goose', source) ~ 'taiga/tundra bean-goose',
                                      grepl('surf scooter', source) ~ 'surf scoter',
                                      grepl('northern goshawk', source)  ~ 'eurasian goshawk', # add & for country
                                      grepl('rosss goose', source) ~ "ross's goose",
                                      grepl('turnstone', source) & location == 'Netherlands' ~ 'ruddy turnstone',
                                      grepl('madarin duck', source) ~ 'mandarin duck',
                                      grepl('bar headed goose', source) ~ 'bar-headed goose',
                                      grepl('shoveler', source) ~ 'northern shoveler',
                                      grepl('white-fronted goose', source) ~ 'greater/lesser white-fronted goose',
                                      grepl('mallard duck',source) ~ 'mallard',
                                      grepl('chicken|layer hen|laying hen', source) ~ 'red junglefowl (domestic type)',
                                      grepl('domestic goose', source) ~ 'domestic goose sp. (domestic type)', 
                                      grepl('domestic duck', source) ~ 'mallard (domestic type)', 
                                      grepl('turkey', source) & location != 'Mexico' ~ 'wild turkey',
                                      
                                      # latin names
                                      grepl('anas crecca', source) ~ 'green-winged teal',
                                      grepl('anas platyrhynchos', source) ~ 'mallard',
                                      grepl('anser albifrons', source) ~ 'greater white-fronted goose',
                                      grepl('anser brachyrhynchus', source) ~ 'pink-footed goose',
                                      grepl('anser anser', source) ~ 'greylag goose',
                                      grepl('anser fabalis', source) ~ 'taiga/tundra bean-goose',
                                      grepl('chlidonias hybrida', source) ~ 'whiskered tern',
                                      grepl('cygnus columbianus', source) ~ 'tundra swan',
                                      grepl('falco peregrinus', source) ~ 'peregrine falcon',
                                      grepl('anser fabalis', source) ~ 'taiga/tundra bean-goose',
                                      
                                      #generics
                                      grepl('\\bswan\\b', source) ~ 'swan sp.',
                                      grepl('\\begret\\b', source) ~ 'white egret sp.',
                                      grepl('\\bcrane\\b', source) ~ 'crane sp.',
                                      grepl('\\bgoose\\b|wild geese|wild goose', source) ~ 'goose sp.',
                                      grepl('\\begret\\b', source) ~ 'egret sp.',
                                      grepl('waterfowl', source) ~ 'waterfowl sp.',
                                      grepl('peacock', source) ~ 'green peafowl',
                                      grepl('quail', source) & location != 'America' ~ 'old world quail sp.', # requires better location
                                      grepl('shorebird', source) ~ 'shorebird sp.',
                                      grepl('black-backed gull|\\bgull\\b', source) ~ 'gull sp.',
                                      grepl('\\bcrow\\b', source) ~ 'crow/raven sp.',
                                      grepl('\\bduck\\b', source) ~ 'duck sp.',
                                      grepl('\\bgrebe\\b', source) ~ 'grebe sp.',
                                      grepl('\\bowl\\b', source) ~ 'owl sp.',
                                      grepl('\\bpelican\\b', source) ~ 'pelican sp.',
                                      grepl('\\bpheasant\\b', source) ~ 'pheasant sp.',
                                      grepl('\\bpigeon\\b', source) ~ 'pigeon/dove sp.',
                                      grepl('\\bteal\\b', source) ~ 'teal sp.',
                                      
                                      
                                      
                                      
                                      

                                      
                                      .default = source)) %>%
  #rename(primary_com_name = source) %>%
  #mutate(sci_name = primary_com_name) %>%
  left_join(birds, join_by(primary_com_name)) %>%

# family
  mutate(family = case_when(
    primary_com_name %in% c('en', 'env', 'environment', 'environmental', 'water', 'wild bird feces') ~ 'environment',
    primary_com_name %in% c('skunk') ~ 'mephitidae',
    primary_com_name %in% c('tanuki', 'ezo red fox') ~ 'canidae',
    primary_com_name %in% c('rattus norvegicus') ~ 'muridae',
    primary_com_name %in% c('feline') ~ 'felidae',
    
  .default = family)) %>%
  
  # order
  mutate(order = case_when(
    primary_com_name %in% c('en', 'env', 'environment', 'environmental', 'water', 'wild bird feces') ~ 'environment',
    primary_com_name %in% c('skunk', 'tanuki', 'ezo red fox', 'feline') ~ 'carnivora',
    primary_com_name %in% c('rattus norvegicus') ~ ' Rodentia',
    .default = order)) %>%
  
  # species group
  mutate(species_group = case_when(
    primary_com_name %in% c('en', 'env', 'environment', 'environmental', 'water', 'wild bird feces') ~ 'Environment',
    primary_com_name %in% c('skunk', 'tanuki', 'ezo red fox', 'feline') ~ 'Carnivores',
    primary_com_name %in% c('rattus norvegicus') ~ 'Rodents',
    primary_com_name %in% c('aquatic bird', 'avian', 'bird', 'wild bird', 'wild bird ') ~ 'Unknown bird',
    .default = species_group)) %>%
  
  
  # Location (iso name and subdivision)

  

  