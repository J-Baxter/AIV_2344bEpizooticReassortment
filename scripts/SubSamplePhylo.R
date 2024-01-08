# Basic filtering of maximum likelihood phylogenies
# A function to create an initial subsample of nodes based on the following broad criteria:
# 1. Genetic distance - select a single example of trait/location/date from a clade in which 
#     all descending branches are less than 0
# 2. Remove outliers (>2sd from the root-to-tip regression)


# Filter phylogenetic tree according to branch length
FilterBL <- function(phylo,  bl_threshold = 0){
  require(ape)
  
  return(tiplist)
}

# filter phylogenetic tree accoding to root-to-tip regression
FilterR2T <- function(phylo,  sd_threshold = 0){
  require(ape)
  
  return(tiplist)
}

# root to tip regression plots
library(ape)
library(tidyverse)
library(treeio)
library(TreeTools)
library(Rcpp)

sourceCpp("./scripts/getTaxonomyForName.cpp") # ~8mins
sourceCpp('./scripts/getLocation.cpp') #~3mins



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

Main <- function(tiplabels){
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
    
    # Format columna names
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
    select(-c(source, location)) %>%
    
    # Get reassortant data from metadata
    left_join(reassortants,
              by = join_by(isolate.id)) %>%
    
    # Get missing dates from metadata csv 
    left_join(metadata %>% 
                select(isolate_id, date_year_month), 
              join_by(isolate.id == isolate_id)) %>%
    mutate(collection.dateweek = coalesce(collection.dateweek, date_year_month)) %>%
    mutate(collection.tipdate = case_when(is.na(collection.date) ~ collection.dateweek,
                                          .default = as.character(collection.date))) %>%
    
    # Create tipnames
    # subtype|isolatename|hostorder|country|reassortant|decimaldate
    #H5N1_Avian_AB188824_A/chicken/Kyoto/3/2004|2004.500
    unite(tipnames, 
          virus.subtype,
          isolate.id,
          host.order,
          collection.country.name,
          cluster.genome,
          collection.tipdate,
          sep = '|',
          remove = FALSE)%>%
    mutate(tipnames = gsub(' ', '_', tipnames))
  
  return(out)
}

# Load data
treefiles <- list.files(path = './data/phylo_ml/',
                        recursive = FALSE,
                        include.dirs = FALSE, 
                        full.names = TRUE)

segnames <- str_split(treefiles,  '_') %>% 
  lapply(., tail, n = 1) %>% 
  unlist() %>%
  toupper()

trees <- lapply(treefiles, read.newick) %>% 
  setNames(segnames)

tipnames <- lapply(trees, TipLabels) %>% 
  setNames(segnames)

formatted_metadata <- lapply(tipnames, Main)%>% 
  setNames(segnames)



Myfunc <- function(trees, metadatalist, x){
  
  trees[[x]][['tip.label']] <- metadatalist[[x]][['tipnames']]
  
  return(trees[[x]])
}

renamed_phylos <- lapply(segnames, Myfunc, trees = trees, metadatalist = formatted_metadata)

cleaned_phylos <- lapply(renamed_phylos, CleanPhylo)
mapply(write.tree, cleaned_phylos, paste0(segnames, '.tree'))

#test_main <- Main(tipnames$HA)


metadata <- read_csv('./data/metadata/h5_metadata.csv')
reassortants <- metadata %>% 
  select(c(isolate_id, contains('cluster'))) %>%
  rename_with(.fn = ~gsub('_', '.', .x))

# Filter out sequences with no date
CleanPhylo <- function(phylo){
  out <- drop.tip(phylo,
                  TipLabels(phylo)[grep('NA$', TipLabels(phylo))])
  
  return(out)
}

#rooted_phylos <- lapply(renamed_phylos, 
                       # CleanPhylo) %>%
 # mapply(rtt, .,
       #  tipdates,
       #  SIMPLIFY = F) %>%
 # mapply(function(x,y) ggtree(x) + geom_tiplab(is.na(y$cluster.genome)),
        # ., 
        # formatted_metadata, 
        # SIMPLIFY = F)
  
###################
CleanandRoot <- function(phylo, metadata){
  require(ape)
  require(TreeTools)
  require(adephylo)
  
  stopifnot(Ntip.phylo(phylo) == nrow(metadata))
  
  phylo_dropped <- CleanPhylo(phylo = phylo)
  
  r2t <- rtt(phylo_dropped,
             objective = 'rms', 
             tip.dates = metadata$collection.datedecimal[!is.na(metadata$collection.datedecimal)])
  
  rttdist <- adephylo::distRoot(r2t)
  
  df <-  cbind.data.frame('dist' = rttdist, 'tipnames' = names(rttdist)) %>% 
    left_join(metadata)
  
  return(df)
  
}


# Infer best-fitting root and get rtt distances (returns list of dataframes)
rtt_dist <- mapply(CleanandRoot,
                   phylo = renamed_phylos,
                   metadata = formatted_metadata,
                   SIMPLIFY = FALSE) %>%
  setNames(segnames) #%>%
  #bind_rows(.id = 'segment')

#pb2_tree_dropped <- drop.tip(pb2_tree, TipLabels(pb2_tree)[grep('NA$', TipLabels(pb2_tree))])
#r2t <- rtt(pb2_tree_dropped )
#rttdist <- adephylo::distRoot(r2t)

#test_df <-  cbind.data.frame('dist' = rttdist, 'tipnames' = names(rttdist)) %>% 
#left_join(df )


rtt_reg <- function(x){
  my_model <- lm(dist ~ collection.datedecimal, data = x)
  rate = coef(my_model)[2]
  x_intercept <- -coef(my_model)[1] / coef(my_model)[2] 
  
  out <- cbind.data.frame(rate, x_intercept) %>% as_tibble()
  return(out)
}


t <- rtt_dist %>% bind_rows(., .id = 'segment') %>%
  left_join(ratesandintercept)

# Calculate gradient (evolutionary rate) and intercept (tmrca)
ratesandintercept <- lapply(rtt_dist, rtt_reg) %>%
  bind_rows(.id = 'segment') %>%
  left_join(ratesandintercept)

# Plot
ggplot(t, aes(x = collection.datedecimal, y = dist)) +
  geom_point(aes(colour = !is.na(cluster.profile))) + 
  scale_y_continuous('Genetic Distance', 
                     expand = expansion(mult = c(0,0.05))) + 
  scale_x_continuous('Time') +
  scale_colour_brewer(
    'Reassortment',
    palette = 'Paired', 
    labels = c('FALSE' = 'Non Reassortment',
               'TRUE' = 'Reassortment'))+
  geom_point(aes(x = x_intercept, y = 0)) + 
  geom_smooth(method = 'lm', 
              fullrange = T,
              colour = 'black') +
  coord_cartesian(ylim = c(0, NA)) +
  geom_text(aes(x = x_intercept,
                y = Inf,
                label = paste('x intercept =', round(x_intercept, 2))),
            size = 3,
            vjust= +3, 
            hjust = -0.2, 
            data = ratesandintercept) +
  geom_text(aes(x = x_intercept,
                y = Inf,
                label = paste('rate =', formatC(rate, format = "e", digits = 2))), 
            size = 3,
            vjust= +5,
            hjust = -0.27,
            data = ratesandintercept) + 
  facet_wrap(.~segment, scales = 'free') +
  my_theme + 
  theme(legend.position = 'bottom')


###################
#t <- rooted_phylos %>%
 # lapply(., as_tibble) %>%
 # {lapply(seq_along(along.with = .), 
  #        function(i)  left_join( x = .[[i]], 
  #                                y = formatted_metadata[[i]], 
  #                                join_by(label == tipnames)))} %>% 
  #lapply(., as.treedata) %>%
 # lapply(., function(x) ggtree(x) + geom_tippoint(aes(colour = is.na(cluster.genome))) +
  #         scale_colour_brewer(
  #           'Reassortment',
  #           palette = 'Paired', 
  #           direction = -1,
 #           labels = c('TRUE' = 'Non Reassortment',
  #                      'FALSE' = 'Reassortment')) + theme(legend.position = 'none')) %>%
 # cowplot::plot_grid(plotlist = ., align = 'hv', labels = segnames, ncol = 3)


#t_mat <- rooted_phylos %>%
#  lapply(., as_tibble) %>%
#  {lapply(seq_along(along.with = .), 
#          function(i)  left_join( x = .[[i]], 
 #                                 y = formatted_metadata[[i]], 
 #                                 join_by(label == tipnames)))} %>% 
#  lapply(., as.treedata) %>%
#  lapply(., function(x) ggtree(x) + 
 #          theme(legend.position = 'none') %>%
#           gheatmap(., is.na(x@data$cluster.genome),
 #                   offset=8, 
#                    width=0.6, 
#                    colnames=FALSE) ) %>%
#  cowplot::plot_grid(plotlist = ., align = 'hv', labels = segnames, ncol = 3)
###################
# Remove tip names identified as problems







seqstoremove <- list('HA' = 'H5N8|EPI_ISL_4061715|Falconiformes|japan|10|2021-01-28',
                     'MP' = c('NA|EPI_ISL_13157449|Anseriformes|bangladesh|NA|2022-02-25',
                              'H9N2|EPI_ISL_6785027|Anseriformes|south_korea|NA|2020-10-07',
                              'H10N3|EPI_ISL_1096144|environment|china|NA|2019-06-18',
                              'H3N8|EPI_ISL_16166669|Anseriformes|china|NA|2019-10-08',
                              'H3N6|EPI_ISL_16212299|Anseriformes|mongolia|NA|2021-08-01',
                              'H4N6|EPI_ISL_6784656|Anseriformes|south_korea|NA|2020-10-13',
                              'H5N2|EPI_ISL_1760455|Anseriformes|china|286|2020-11-10',
                              'H5N8|EPI_ISL_7381026|Galliformes|china|1|2021-03',
                              'H5N6|EPI_ISL_6758007|Anseriformes|china|61|2021-04',
                              'H5N6|EPI_ISL_6772739|Galliformes|china|127|2021-08',
                              'H5N8|EPI_ISL_7380623|Anseriformes|china|1|2021-03',
                              'H5N6|EPI_ISL_6757631|Anseriformes|china|61|2021-04',
                              'H5N6|EPI_ISL_6772736|Galliformes|china|127|2021-08'),
                     'NA1' = NULL,
                     'NA2' = c('H10N2|EPI_ISL_18414927|Anseriformes|mongolia|NA|2022-05-01',
                               'H10N2|EPI_ISL_18414926|Anseriformes|mongolia|NA|2022-05-01',
                               'H3N2|EPI_ISL_337398|Anseriformes|russia|NA|2018-09-15',
                               'H4N2|EPI_ISL_503356|Anseriformes|mongolia|NA|2019-09-20'),
                     
                     'NA3' = 'NA|EPI_ISL_328982|Charadriiformes|netherlands|NA|2016-07-15',
                     
                     'NA6' = c('H5N6|EPI_ISL_6914131|Galliformes|china|61|2021-06-10',
                               'H5N6|EPI_ISL_17261964|environment|china|NA|2022-03-01',
                               'H5N6|EPI_ISL_17262043|environment|china|NA|2022-03-01',
                               'H5N6|EPI_ISL_17262174|environment|china|NA|2023-03-01' ),
                     'NA8' = 'H5N8|EPI_ISL_7380813|Anseriformes|china|1|2021-03',

                     'NP'= c(
                       'H5N1|EPI_ISL_18132572|Anseriformes|united_states|NA|2022-04-27',
                       'H5N1|EPI_ISL_17051418|Falconiformes|canada|NA|2022-10-28',
                       'H5N1|EPI_ISL_17051394|Pelecaniformes|canada|9|2022-08-12',
                       'H5N1|EPI_ISL_16023374|Carnivora|canada|9|2022-07-19',
                       'H5N1|EPI_ISL_18132829|Galliformes|united_states|NA|2022-04-19',
                       'H5N1|EPI_ISL_17051393|Pelecaniformes|canada|9|2022-08-12',
                       'H5N1|EPI_ISL_17964876|Charadriiformes|united_states|9|2023-02-14',
                       'H5N6|EPI_ISL_6772898|Anseriformes|china|209|2021-08',
                       'H5N1|EPI_ISL_17964946|Galliformes|united_states|9|2023-04-19',
                       'H5N1|EPI_ISL_17051423|Anseriformes|canada|9|2022-11-09',
                       'H5N1|EPI_ISL_17051343|Pelecaniformes|canada|NA|2022-06-02',
                       'H5N1|EPI_ISL_16023375|Carnivora|canada|9|2022-07-19',
                       'H5N1|EPI_ISL_16183504|Galliformes|canada|9|2022-05-16',
                       'H5N1|EPI_ISL_15647835|Anseriformes|south_korea|9|2022-10-19',
                       'H5N1|EPI_ISL_17051342|Pelecaniformes|canada|NA|2022-06-02',
                       'H5N1|EPI_ISL_17964894|Anseriformes|united_states|NA|2023-02-27',
                       'H5N1|EPI_ISL_17051473|Passeriformes|canada|NA|2022-12-30',
                       'H5N1|EPI_ISL_16271855|Anseriformes|united_states|9|2022-11-09',
                       'H5N1|EPI_ISL_17051341|Pelecaniformes|canada|9|2022-06-02',
                       'H5N1|EPI_ISL_18132418|Anseriformes|united_states|NA|2022-04-14',
                       'H5N1|EPI_ISL_17964944|Accipitriformes|united_states|9|2023-03-31',
                       'H5N1|EPI_ISL_17051458|Anseriformes|canada|NA|2022-12-12',
                       'H5N1|EPI_ISL_17690225|Anseriformes|united_states|9|2022-09-23',
                       'H5N1|EPI_ISL_17260681|Carnivora|united_states|9|2022-06-15',
                       'H5N1|EPI_ISL_18132826|Galliformes|united_states|NA|2022-04-19',
                       'H5N1|EPI_ISL_16023376|Carnivora|canada|9|2022-05-31',
                       'H5N1|EPI_ISL_17051490|Strigiformes|canada|NA|2023-01-16',
                       'H5N1|EPI_ISL_17051351|Anseriformes|canada|NA|2022-05-31',
                       'H5N1|EPI_ISL_17964858|Accipitriformes|united_states|9|2023-02-14',
                       'H5N8|EPI_ISL_7381183|Anseriformes|china|279|2021-04',
                       'H5N1|EPI_ISL_17424649|Anseriformes|united_states|9|2022-09-01',
                       'H5N1|EPI_ISL_17051485|Anseriformes|canada|NA|2023-01-16'),
                     'NS'= NULL,
                     'PA'= NULL,
                     'PB1' = NULL,
                     'PB2' = NULL) 


dropped_phylos <- mapply(drop.tip,
                             renamed_phylos,
                             seqstoremove,
                             SIMPLIFY = F)

dropped_metadata <- mapply(function(data,y) data %>% filter(!tipnames %in% y), 
                           formatted_metadata, 
                           seqstoremove, 
                           SIMPLIFY = F)
MyFunc5 <- function(data, n){
  out <- data %>%
    mutate(collection.subdiv1.code = gsub('^NA$', '', collection.subdiv1.code)) %>% 
    mutate(best_location_code = coalesce(collection.subdiv1.code, collection.country.code)) %>%
    group_by(best_location_code, host.family, cluster.genome, collection.tipdate) %>%
    slice_sample(n = n)
  
  return(out)
}


subset_metadata <- formatted_metadata %>% 
  lapply(., MyFunc5 , n = 1)


# Subset alignments
init_filenames <- list.files('./data/alignments/init_alignments', full.names = T,  pattern="\\.")
init_alignments <- lapply(init_filenames, read.dna, format = 'fasta', as.matrix = T)

# change alignment seqnames

subset_alignment <- init_alignments[[1]][rownames(init_alignments[[1]]) %in% subset_metadata[[1]]$tipnames,]

ReNameAlignment <- function(alignment, data){
  isolates <- str_extract(rownames(alignment), "([^|]*)\\|")%>%
    str_replace_all("\\|", "")
  
  new_seqnames <- sapply(isolates, function(x) data$tipnames[data$isolate.id %in% x]) %>% 
    as.vector() 

  rownames(alignment) <-  new_seqnames

  return(alignment)
}



SubsetAlignment <- function(alignment, data){
  out <- alignment[rownames(alignment) %in% data$tipnames,]
  return(out)
}


subsetted_alignments <- init_alignments %>%
  mapply(ReNameAlignment, 
         .,
         formatted_metadata, 
         SIMPLIFY = FALSE) %>%
  mapply(SubsetAlignment,
         .,
         subset_metadata,
         SIMPLIFY = FALSE) 
  
cleaned_filenames <- init_filenames %>%
  gsub('.trim.fasta', '.trim_subsampled.fasta',. ) %>%
  gsub('/init_alignments/', '/subsampled_alignments/',. )


mapply(write.dna, 
       subsetted_alignments , 
       cleaned_filenames,
       format = 'fasta')
  
