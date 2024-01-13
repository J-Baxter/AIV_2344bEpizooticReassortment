####################################################################################################
# Basic filtering of maximum likelihood phylogenies

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
sourceCpp('./scripts/getLocation.cpp') #~3mins

# Functions
FormatTipData <- function(x) {
  isolate.id <- x[grep("EPI_ISL_*", x)]
  
  subtype <- x[grep("H[[:digit:]]{1,}N[[:digit:]]*", x)] %>%
    str_extract(., "H[[:digit:]]{1,}N[[:digit:]]*")
  
  date <- x[grep("^[[:digit:]][[:digit:]][[:digit:]][[:digit:]]-[[:digit:]][[:digit:]]-[[:digit:]][[:digit:]]", x)]
  
  isolate.name <- x[grep("/[^/]*/[^/]*/", x)]
  
  if (identical(isolate.id, character(0))) {
    isolate.id <- NA
  }
  
  if (identical(subtype, character(0))) {
    subtype <- NA
  }
  
  if (identical(date, character(0))) {
    date <- NA
  }
  
  if (identical(isolate.name, character(0))) {
    isolate.name <- NA
  }
  
  out <- cbind.data.frame(
    "isolate.id" = isolate.id,
    "subtype" = subtype,
    "date" = date,
    "isolate.name" = isolate.name
  ) %>%
    distinct()
  
  
  return(out)
}


ExtractMetadata <- function(tiplabels){
  out <- tiplabels %>%
    gsub("\\|$", "", .) %>%
    sapply(., str_split_1, pattern = "\\|") %>%
    lapply(., FormatTipData) %>%
    bind_rows() %>%
    as_tibble() %>%
    
    # Format Dates
    mutate(date = parsedate::parse_date(date) %>%
             as.Date()) %>%
    mutate(decimal.date = format(round(decimal_date(date), 2), 
                                 nsmall = 2) ) %>%
    mutate(decimal.date = suppressWarnings(as.double(decimal.date))) %>%
    mutate(week.date = format(date, "%Y-%V")) %>%
    mutate(month.date = format(date, "%Y-%m")) %>%
    
    # Format isolate name
    mutate(isolate.name = gsub("_(\\d{4})|(\\d{4})_", "\\1", perl = TRUE, isolate.name)) %>%
    mutate(isolate.name = case_when(!grepl("/(\\d{4})$", isolate.name) & !grepl("/$", isolate.name) ~ 
                                      paste0(isolate.name, "/NA"),
                                    !grepl("/(\\d{4})$", isolate.name) & grepl("/$", isolate.name) ~ paste0(isolate.name, "NA"),
                                    .default = isolate.name
    )) %>%
    
    # Extract data from isolate name
    separate_wider_delim(isolate.name,
                         delim = "/",
                         names = c("virus_species", "source", "location", "id_unsure", "year"),
                         cols_remove = F,
                         too_few = "align_end",
                         too_many = "debug"
    ) %>%
    
    # Drop extraneous cols and format
    dplyr::select(-c(year, contains("isolate.name_"))) %>%
    mutate(virus_species = case_when(str_length(source) == 1 ~ source,
                                     .default = virus_species
    )) %>%
    mutate(source = case_when(str_length(source) == 1 ~ NA,
                              .default = source
    )) %>%
    
    # Format source before taxa allocation
    mutate(source = gsub("_", " ", tolower(source))) %>%
    mutate(source = str_trim(source)) %>%
    mutate(source= gsub('gray', 'grey', source)) %>%
    mutate(primary_com_name = case_when(
      grepl("common|eurasian teal|green-winged-teal", source) ~ "green-winged teal",
      grepl("spot-billed duck", source) ~ "eastern spot-billed duck",
      grepl("crested grebe", source) ~ "great crested grebe",
      grepl("eastern curlew", source) ~ "far eastern curlew",
      grepl("bean goose", source) ~ "taiga/tundra bean-goose",
      grepl("surf scooter", source) ~ "surf scoter",
      grepl("northern goshawk|^goshawk$", source) ~ "eurasian goshawk", # add & for country
      grepl('sparrowhawk', source) & grepl('England', location) ~ 'eurasian sparrowhawk',
      grepl("coopers s hawk", source) ~ "cooper's hawk",
      grepl("rosss goose", source) ~ "ross's goose",
      grepl("turnstone", source) & location == "Netherlands" ~ "ruddy turnstone",
      grepl("madarin duck", source) ~ "mandarin duck",
      grepl("bar headed goose", source) ~ "bar-headed goose",
      grepl('grey goose|graylag goose', source) ~ 'greylag goose',
      grepl("shoveler", source) ~ "northern shoveler",
      grepl("white-fronted goose", source) ~ "greater/lesser white-fronted goose",
      grepl("mallard duck", source) ~ "mallard",
      grepl("chicken|layer hen|laying hen|broiler|poultry|chichen", source) ~ "red junglefowl (domestic type)",
      grepl("domestic goose|pomeranian goose", source) ~ "domestic goose sp. (domestic type)",
      grepl("domestic duck|pekin duck|mule duck", source) ~ "mallard (domestic type)",
      grepl('^buzzard$', source) ~ 'common buzzard',
      grepl('knot wader', source) ~ 'red knot',
      grepl("turkey", source) & location != "Mexico" ~ "wild turkey",
      grepl('^kestrel$', source) & grepl('Germany|England|Italy', location) ~ 'eurasian kestrel',
      grepl('brent goose', source) ~ 'brant',
      grepl('eagle owl', source) ~'eurasian eagle-owl',
      grepl('jungle crow', source) & location == 'Japan' ~ 'large-billed crow',
      grepl('^magpie$' ,source) & location == 'Idaho' ~ 'black-billed magpie',
      grepl('oystercatcher', source) & grepl('Germany', location) ~ 'eurasian oystercatcher',
      grepl('canade goose', source) ~ 'canada goose',
      grepl('european herring gull', source) ~ 'herring gull',
      grepl('towny owel', source) ~ 'tawny owl',
      grepl('gadwall duck', source) ~ 'gadwall',
      grepl('lesser snow goose', source) ~ 'snow goose',
      grepl("pink footed goose", source) ~ "pink-footed goose",
      grepl("european wigeon", source) ~ "eurasian wigeon",
      source == 'sea eagle' & grepl('Norway', location)~"white-tailed eagle",
      grepl('guinea fowl', source) & grepl('Italy|Estonia', location) ~ "helmeted guineafowl (domestic type)",
      grepl('western jackdaw', source) ~ 'eurasian jackdaw',
      grepl('^gannet$', source) & grepl('Scotland', location) ~ "northern gannet",
      grepl('spotbill duck', source) ~ "indian spot-billed duck",
      grepl("whistling duck", source) ~ "whistling-duck sp.",
      grepl('pink footed goose', source) ~ 'pink-footed goose',
      grepl('falcated teal', source) ~ 'falcated duck',
      grepl('^pintail$', source) ~ 'northern pintail',
      grepl('grey plover', source)~ "black-bellied plover",
      grepl('greenwing duck', source) ~ "green-winged teal",
      grepl('franklins gull', source) ~ "franklin's gull",
      
      
      
      # latin names
      grepl("anas crecca", source) ~ "green-winged teal",
      grepl('tyto alba', source) ~ 'barn owl',
      grepl('gallus gallus', source) ~ 'red junglefowl (domestic type)',
      grepl("anas platyrhynchos", source) ~ "mallard",
      grepl("anser albifrons", source) ~ "greater white-fronted goose",
      grepl("anser brachyrhynchus", source) ~ "pink-footed goose",
      grepl("anser anser", source) ~ "greylag goose",
      grepl("anser fabalis", source) ~ "taiga/tundra bean-goose",
      grepl("chlidonias hybrida", source) ~ "whiskered tern",
      grepl("cygnus columbianus", source) ~ "tundra swan",
      grepl("falco peregrinus", source) ~ "peregrine falcon",
      grepl("anser fabalis", source) ~ "taiga/tundra bean-goose",
      grepl('larus argentatus', source) ~ 'herring gull',
      grepl('pica pica', source) ~ 'eurasian magpie',
      grepl("ciconia ciconia", source) ~ 'white stork',
      grepl("buteo buteo", source) ~ 'common buzzard',
      grepl('podiceps cristatus', source) ~ 'great crested grebe',
      grepl('anas poecilorhyncha|anas poecilohyncha', source) ~ "indian spot-billed duck",
      grepl('garrulus glandarius', source) ~ 'eurasian jay',
      grepl('columba palumbus', source) ~ "common wood-pigeon",
      grepl('branta canadensis', source) ~ 'canada goose',
      grepl('cygnus olor', source) ~"mute swan",
      grepl('spatula clypeata', source) ~ "northern shoveler",
      grepl('tadorna tadorna', source) ~ "common shelduck",
      grepl('anseriformes', source) ~ 'waterfowl sp.',
      grepl('sterna hirundo', source) ~ "common tern",
      
      grepl('^seal$', source) ~ 'seal sp.', 
      grepl('^mink$', source) ~ 'mink sp.',
      grepl('^otter$', source) & grepl('Scotland', location) ~ 'eurasian otter',
      grepl('badger', source) ~ 'european badger',
      
      
      # generics
      grepl("^swan$", source) ~ "swan sp.",
      grepl("^cormorant$", source) ~ 'cormorant sp.',
      grepl("^egret$", source) ~ "white egret sp.",
      grepl("^crane$|black-throated crane", source) ~ "crane sp.",
      grepl("^goose$|wild geese|wild goose", source) ~ "goose sp.",
      grepl('^stork$', source) ~ 'stork sp.',
      grepl('^ostrich$' ,source) ~ "common ostrich",
      grepl("^egret$", source) ~ "egret sp.",
      grepl('lapwing', source) ~ 'lapwing sp.',
      grepl('^pochard$', source) ~ "aythya sp.",
      grepl('coot', source) ~ 'coot sp.',
      grepl("waterfowl", source) ~ "waterfowl sp.",
      grepl("peacock|peafowl", source) ~ "indian peafowl (domestic type)",
      grepl("quail", source) & location != "America" ~ "old world quail sp.", # requires better location
      grepl("shorebird|sandpiper", source) ~ "shorebird sp.",
      grepl("black-backed gull|^gull$|seagull", source) ~ "gull sp.",
      grepl("^crow$", source) ~ "crow/raven sp.",
      grepl("^duck$|wild duck|migratory duck", source) ~ "duck sp.",
      grepl("^grebe$", source) ~ "grebe sp.",
      grepl("^owl$", source) ~ "owl sp.",
      grepl('^hawk$', source) ~ "hawk sp.",
      grepl('^eagle$', source) ~ 'eagle sp.',
      grepl('steamer duck', source) ~ "steamer-duck sp.",
      grepl('^falcon$', source) ~ "falcon sp.",
      grepl("^pelican$", source) ~ "pelican sp.",
      grepl("^pheasant$", source) ~ "pheasant sp.",
      grepl("^pigeon$", source) ~ "pigeon/dove sp.",
      grepl("^teal$", source) ~ "teal sp.",
      grepl('^ibis$', source) ~ 'ibis sp.',
      source %in% c("aquatic bird", "avian", "bird", "wild bird", "wild bird", "backyard bird") ~ "bird sp.",
      source %in% c("en", "env", "environment", "enviroment", "environmental", 'environment sample',
                    "water", "wild bird feces") ~ "environment",
      grepl("red fox|^fox$", source) ~ "red fox",
      grepl("attus norvegicus", source) ~ "brown rat",
      .default = source
    )) %>%
    
    # Get host taxonomy
    rowwise() %>%
    mutate(taxonomy = getTaxonomyForName(primary_com_name)) %>%
    as_tibble() %>%
    unnest(taxonomy) %>%
    rename_with(.fn = ~ tolower(.x)) %>%
    
    # binary host information
    mutate(is.domestic = case_when(grepl("domestic type", sci.name) ~ "domestic",
                                   .default = "wild"
    )) %>%
    mutate(is.bird = case_when(grepl("ormes$", order) ~ "bird",
                               .default = "other"
    )) %>%
    
    # Location (iso name and subdivision)
    #mutate(location = tolower(location)) %>%
    mutate(location = gsub("_", " ", tolower(location))) %>%
    
    mutate(location = case_when(
      
      # Korea
      grepl("korea", location) ~ "south korea",
      
      # PRC
      grepl("yunnan.*", location) ~ "yunnan",
      grepl('rongcheng', location) ~ 'shandong',
      grepl("inner mongolia|tumuji", location) ~ "nei mongol",
      grepl("hong kong|hongkong", location) ~ "xianggang (hong-kong)",
      grepl("xizang$", location) ~ "xizang",
      grepl("ningxia$", location) ~ "ningxia",
      grepl("xinjiang$", location) ~ "xinjiang",
      grepl("guangxi$|guilin|hechi", location) ~ "guangxi",
      grepl("tibet", location) ~ "xizang",
      grepl("hangzhou|nanji$", location) ~ "zhejiang",  # confirm nanji
      grepl("sanmenxia|dongting", location) ~ "henan",
      grepl("heinan", location) ~ "hainan", # confirm 
      grepl("wuhan", location) ~ "hebei",
      grepl('^fujain$', location) ~ 'fujian',
      grepl('suzhou|xuzhou', location) ~ 'jiangsu',
      grepl("changsha", location) ~ "hunan",
      grepl('tianjing', location) ~ 'tianjin',
      grepl('qingyuan|foshan', location) ~ 'guangdong',
      grepl("northern china", location) ~ "beijing", # Interpretation from paper - check with Lu
      grepl("southwestern china", location) ~ "shanxi", # Interpretation from paper - check with Lu
      grepl('eastern china', location) ~ 'zhejiang', #Inferred - if we need to be 100% sure then change to country level
      
      
      # Canada
      grepl("ab", location) ~ "alberta",
      grepl("bc", location) ~ "british columbia",
      grepl('mb', location) ~'manitoba',
      
      # Russian Federation
      grepl('russian federation', location) ~ 'russia',
      grepl("astrakhan", location) ~ "astrakhanskaya oblast'",
      grepl("novosibirsk region|chany lake|^chany$", location) ~ "novosibirskaya oblast'",
      grepl("omsk.*", location) ~ "omskaya oblast'",
      grepl("rostov-on-don", location) ~ "rostovskaya oblast'",
      grepl("russia primorje", location) ~ "primorskiy kray",
      grepl("sakhalin", location) ~ "sakhalinskaya oblast'",
      grepl("yakutia", location) ~ "sakha, respublika [yakutiya]",
      grepl("kostroma", location) ~ "kostromskaya oblast'",
      grepl("amur region", location) ~ "amurskaya oblast'",
      grepl("buryatia", location) ~ "buryatiya, respublika",
      grepl('north ossetia-alania', location) ~ 'severnaya osetiya-alaniya, respublika',
      grepl('dagestan', location) ~ 'dagestan, respublika',
      grepl('stavropol', location) ~ "stavropol'skiy kray",
      grepl('krasnodar', location) ~ "krasnodarskiy kray", 
      grepl('tyumen', location) ~ "tyumenskaya oblast'",
      grepl("magadan", location) ~ "magadanskaya oblast'",
      grepl('saratov', location) ~ "saratovskaya oblast'",
      grepl('kurgan', location) ~ "kurganskaya oblast'",
      grepl('central russia' ,location) ~ 'russia',
      
      # USA
      grepl("north dakota", location) ~ "north dakota",
      
      # UK
      grepl("england", location) ~ "england and wales",
      # Japan
      grepl("tsukuba", location) ~ "ibaraki",
      
      # Spain
      grepl("castillalamancha", location) ~ "castilla-la mancha",
      
      # Germany 
      location == 'germany-mv' ~ 'mecklenburg-vorpommern',
      location == 'germany-nw' ~ 'nordrhein-westfalen',
      location == 'germany-be' ~ 'berlin',
      location == 'germany-by' ~ 'bayern',
      location =='germany-st' ~ 'sachsen-anhalt',
      location == 'germany-bb' ~ 'brandenburg',
      location == 'germany-ni' ~ 'niedersachsen',
      location == 'germany-sh' ~ 'schleswig-holstein',
      location == 'germany-sn' ~ 'sachsen',
      location == 'germany-he' ~ 'hessen',
      location == 'germany-hh' ~ 'hamburg',
      location == 'germany-bw' ~ 'baden-württemberg',
      location == 'germany-th' ~ 'thüringen',
      location == 'germany-rp' ~ 'rheinland-pfalz',
      
      
      grepl('bosnia and herzegovina', location) ~ 'bosnia-herzegovina',
      
      grepl('sva[:0-9:]', location) ~ 'NA',
      
      # Indonesia
      grepl("hulu sungai utara|banjarbaru", location) ~ "kalimantan selatan",
      grepl("republic of georgia", location) ~ "georgia",
      grepl("czech republic|czechia", location) ~ "czech republic",
      # grepl('laos', location) ~ "lao people's democratic republic",
      grepl("vietnam", location) ~ "viet nam",
      grepl("quang ninh" , location) ~  "quảng ninh" ,
      .default = location
    )) %>%
    rowwise() %>%
    mutate(loc = getLocation(location)) %>%
    as_tibble() %>%
    unnest(loc) %>%
    
    # Format column names
    dplyr::select(-c(virus_species, id_unsure)) %>%
    dplyr::rename(
      virus.subtype = subtype,
      collection.date = date,
      collection.datedecimal = decimal.date,
      collection.dateweek = week.date,
      collection.region.name = region,
      collection.country.name = name_country,
      collection.country.code = code_country,
      collection.country.lat = lat_country,
      collection.country.long = long_country,
      collection.subdiv1.name = name_subdiv1,
      collection.subdiv1.code = code_subdiv1,
      collection.subdiv1.lat = lat_subdiv1,
      collection.subdiv1.long = long_subdiv1,
      collection.subdiv2.name = name_subdiv2,
      collection.subdiv2.code = code_subdiv2,
      collection.subdiv2.lat = lat_subdiv2,
      collection.subdiv2.long = long_subdiv2,
      host.order = order,
      host.family = family,
      host.sciname = sci.name,
      host.commonname = primary_com_name,
      host.isbird = is.bird,
      host.isdomestic = is.domestic
    ) %>%
    dplyr::relocate(
      virus.subtype,
      isolate.id,
      isolate.name,
      collection.date,
      collection.datedecimal,
      collection.dateweek,
      host.order,
      host.family,
      host.sciname,
      host.commonname,
      host.isbird,
      host.isdomestic,
      collection.region.name,
      collection.country.name,
      collection.country.code,
      collection.country.lat,
      collection.country.long,
      collection.subdiv1.name,
      collection.subdiv1.code,
      collection.subdiv1.lat,
      collection.subdiv1.long,
      collection.subdiv2.name,
      collection.subdiv2.code,
      collection.subdiv2.lat,
      collection.subdiv2.long
    ) %>%
    # Select columns
    select(-c(source, location)) 
  
  
  return(out)
}


MergeReassortantData <- function(data, newdata = newdata){
  
  stopifnot(any('isolate_id' %in% colnames(newdata)))
  
  newdata <- newdata %>% 
    select(c(isolate_id, contains('cluster'), date_year_month)) %>%
    rename_with(.fn = ~gsub('_', '.', .x))
  
  # Get reassortant data from metadata
  out <- left_join(data, 
              newdata,
              by = join_by(isolate.id)) %>%
    
    # Get missing dates from metadata csv
    mutate(collection.dateweek = dplyr::coalesce(collection.dateweek, 
                                          date.year.month)) %>%
    
    mutate(collection.tipdate = case_when(is.na(collection.date) ~ collection.dateweek,
                                          .default = as.character(collection.date))) %>%
    select(-c(date.month.year, month.date))
  
  return(out)
}


MakeTipNames <- function(data){
  out <- data %>%
    unite(tipnames, 
          virus.subtype,
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
treefiles <- list.files(path = './data/phylo_ml/original_mlphylos',
                        recursive = FALSE,
                        include.dirs = FALSE, 
                        full.names = TRUE)

trees <- lapply(treefiles, read.newick) 


# Subset alignments
alignmentfiles <- list.files('./data/alignments/init_alignments',
                             full.names = T, 
                             pattern="\\.")

alignments <- lapply(init_filenames,
                     read.dna, 
                     format = 'fasta',
                     as.matrix = T)


# Segment names
segnames <- str_split(treefiles,  '_') %>% 
  lapply(., tail, n = 1) %>% 
  unlist() %>%
  toupper()

names(alignments) <- segnames
names(trees) <- segnames


# Import existing metadata
reassortant_metadata <- read_csv('./data/metadata/h5_metadata.csv')


# Import list of sequences to remove due to lack of temporal signal
seqstoremove <- read_csv('./data/erroneous_seqs.csv')


####################################################################################################
# Extract metadata and format tipnames
tipnames <- lapply(trees, TipLabels) %>% 
  setNames(segnames) # Only required because we are working backwards from tipnames


# Extract metadata from tipnames
metadata <- lapply(tipnames, ExtractMetadata) %>% 
  lapply(., MergeReassortantData, 
         newdata = reassortant_metadata) %>%
  lapply(., MakeTipNames) %>%
  setNames(segnames) 


# Rename phylogenies and alignments with new-format tipnames
renamed_phylos <- lapply(segnames,
                         ReNamePhylo, 
                         trees = trees, 
                         metadata = metadata)


renamed_alignments <- alignments %>%
  mapply(ReNameAlignment, 
         .,
         metadata, 
         SIMPLIFY = FALSE)


####################################################################################################
# Identify duplicate sequences from alignment

identical_seqs <- lapply(renamed_alignments,
                         FindIdenticalSeqs) %>%
  bind_rows(., .id = 'virus.segment')


####################################################################################################
# Subsample - strict

metadata_subsampled_strict <- metadata %>%
  bind_rows(., .id = 'virus.segment')%>%
  
  # join identical sequence groupings
  left_join(identical_seqs, 
            by = join_by(virus.segment, 
                         tipnames)) %>%
  # create groupings
  mutate(across(contains('collection'), .fns = ~ gsub('^NA$', NA, .x))) %>% 
  mutate(best_location_code = coalesce(collection.subdiv1.code,
                                       collection.country.code)) %>%
  mutate(group = as.factor(group)) %>%
  
  # Sample one virus in each segment that has unique combination of location, host family, 
  # reasssortment, and sequence identity (viruses in the same group are identical)
  
  group_by(group,
           host.family, 
           #best_location_code, #Leave in if relaxed, leave out if strict
           cluster.genome) %>%
  slice_sample(n = 1) %>%
  ungroup() %>%
  
  # as list
  group_split(virus.segment, .keep = FALSE) %>%
  as.list() %>%
  setNames(segnames)


####################################################################################################
metadata_subsampled_relaxed <- metadata %>%
  bind_rows(., .id = 'virus.segment')%>%
  
  # join identical sequence groupings
  left_join(identical_seqs, 
            by = join_by(virus.segment, 
                         tipnames)) %>%
  # create groupings
  mutate(across(contains('collection'), .fns = ~ gsub('^NA$', NA, .x))) %>% 
  mutate(best_location_code = coalesce(collection.subdiv1.code,
                                       collection.country.code)) %>%
  mutate(group = as.factor(group)) %>%
  
  # Sample one virus in each segment that has unique combination of location, host family, 
  # reasssortment, and sequence identity (viruses in the same group are identical)
  
  group_by(group,
           host.family, 
           best_location_code, #Leave in if relaxed, leave out if strict
           cluster.genome) %>%
  slice_sample(n = 1) %>%
  ungroup() %>%
  
  # as list
  group_split(virus.segment, .keep = FALSE) %>%
  as.list() %>%
  setNames(segnames)

####################################################################################################
# Subsample alignments

alignments_subsampled <- renamed_alignments %>%
  mapply(SubsampleAlignment,
         .,
         metadata_subsampled,
         SIMPLIFY = FALSE)


####################################################################################################
# Write alignments and metadata to file
alignmentfiles_subsampled <- alignmentfiles  %>%
  gsub('Re_H5_', '',. ) %>%
  gsub('_Asia_blast_', '_asia_',. ) %>%
  gsub('.trim.fasta$', '_subsampled.fasta',. ) %>%
  gsub('/init_alignments/', '/subsampled_alignments/',. )

mapply(write.dna, 
       alignments_subsampled, 
       alignmentfiles_subsampled ,
       format = 'fasta')


metadatafiles_subsampled <- paste('./data/metadata/clade2344b_asia', 
                                  segnames, 
                                  'subsampled.csv',  
                                  sep = '_')
mapply(write.csv, 
       metadata_subsampled, 
       metadatafiles_subsampled)

  
####################################################################################################
# END #
####################################################################################################