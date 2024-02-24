# update missing data from phylogenetic clusters

# read tree
treefiles <- list.files(path = './data/alignments/subsampled_alignments/2024Jan24/',
                        pattern = '.treefile',
                        recursive = FALSE,
                        include.dirs = FALSE, 
                        full.names = TRUE)

trees <- lapply(treefiles, read.newick) 


# Segment names
segnames <- str_split(treefiles,  '_') %>% 
  lapply(., tail, n = 1) %>% 
  unlist() %>%
  toupper() %>%
  gsub('NA', 'N', .)

names(trees) <- segnames

# Import existing metadata


####################################################################################################
# Extract metadata and format tipnames
phylotipnames <- lapply(trees, TipLabels) %>% 
  setNames(segnames) %>%
  lapply(., function(x)   gsub('\\.', '|', x))
# Only required because we are working backwards from tipnames


######## To be incorporated into main Subsample phylo (106) ########

ImputeCladeandCluster <- function(metadata, seqnames){
  out <- metadata %>%
    mutate(tipnames = ordered(tipnames, levels = seqnames)) %>%
    arrange(tipnames) %>%
    
    # Impute clade
    mutate(clade = case_when(is.na(clade) & na.locf0(clade, fromLast = TRUE) == na.locf0(clade, fromLast = FALSE) ~ na.locf0(clade, fromLast = TRUE), 
                             .default = clade)) %>%
    
    # Impute cluster (column is dependent on segment)
    rename(profile = cluster.profile) %>%
    rename(genome = cluster.genome) %>%
    pivot_longer(contains('cluster'), values_to = 'cluster.number', names_to = 'cluster.segment') %>%
    mutate(cluster.segment = gsub('.*\\.', '', cluster.segment)) %>%
    mutate(cluster.segment = case_when(grepl('^N[:0-9:]', segment) & cluster.segment == 'na' ~ tolower(segment),
                                       .default = cluster.segment)) %>% 
    filter(tolower(segment) == cluster.segment) %>%
    mutate(cluster.number = case_when(is.na(cluster.number) & na.locf0(cluster.number, fromLast = TRUE) == na.locf0(cluster.number, fromLast = FALSE) ~ na.locf0(cluster.number, fromLast = TRUE), 
                                      .default = cluster.number)) 
  
  return(out)
}

test <- mapply(ImputeCladeandCluster, reformatted, phylotipnames, SIMPLIFY = FALSE) %>%
  lapply(., function(x) x %>% mutate(tipnames = gsub( '\\|', '\\.', tipnames))) %>%
  lapply(., function(x) x %>% mutate(cluster.number = paste0('profile', str_pad(cluster.number, 3, pad = "0")))) %>%
  lapply(., function(x) x %>% select(c(tipnames, virus.subtype, host.class, host.order, collection.region.name, clade, tiplocation_lat, tiplocation_lon, cluster.number))) %>%
  lapply(., function(x) x %>% select(where(~n_distinct(., na.rm = TRUE) > 1)))



mapply(write_delim, 
       delim = '\t',
       quote= 'needed',
       test, 
       metadatafiles_tsv)


reformatted[[6]] %>%
  mutate(tipnames = ordered(tipnames, levels = phylotipnames[[6]])) %>%
  arrange(tipnames) %>%
  
  # Impute clade
  mutate(clade = case_when(is.na(clade) & na.locf0(clade, fromLast = TRUE) == na.locf0(clade, fromLast = FALSE) ~ na.locf0(clade, fromLast = TRUE), 
                           .default = clade)) %>%
  
  # Impute cluster (column is dependent on segment)
  rename(profile = cluster.profile) %>%
  rename(genome = cluster.genome) %>%
  pivot_longer(contains('cluster'), values_to = 'cluster.number', names_to = 'cluster.segment') %>%
  mutate(cluster.segment = gsub('.*\\.', '', cluster.segment)) %>%
  mutate(cluster.segment = case_when(grepl('^N[:0-9:]', segment) & cluster.segment == 'na' ~ tolower(segment),
                                     .default = cluster.segment)) %>% 
  filter(tolower(segment) == cluster.segment) %>%
  mutate(cluster.number = case_when(is.na(cluster.number) & na.locf0(cluster.number, fromLast = TRUE) == na.locf0(cluster.number, fromLast = FALSE) ~ na.locf0(cluster.number, fromLast = TRUE), 
                           .default = cluster.number)) 
