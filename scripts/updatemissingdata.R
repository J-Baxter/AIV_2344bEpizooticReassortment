# update missing data from phylogenetic clusters

# read tree
treefiles <- list.files(path = './iqtree',
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
  setNames(segnames) # Only required because we are working backwards from tipnames


test <- reformatted[[1]] 

test$tipnames <- ordered(test$tipnames,  levels =  tipnames[[1]])

test <- test %>%
  arrange(tipnames)


test <- reformatted[[1]] %>%
  mutate(tipnames = ordered(tipnames, levels = phylotipnames[[1]])) %>%
  arrange(tipnames) %>%
  
  # Impute clade
  mutate(clade = case_when(is.na(clade) & na.locf0(clade, fromLast = TRUE) == na.locf0(clade, fromLast = FALSE) ~ na.locf0(clade, fromLast = TRUE), 
                           .default = clade)) %>%
  
  # Impute cluster (column is dependent on segment)
  rename(profile = cluster.profile) %>%
  pivot_longer(contains('cluster'), values_to = 'cluster.number', names_to = 'cluster.segment') %>%
  mutate(cluster.segment = gsub('.*\\.', '', cluster.segment)) %>%
  filter(tolower(segment) == cluster.segment) %>%
  mutate(cluster.number = case_when(is.na(cluster.number) & na.locf0(cluster.number, fromLast = TRUE) == na.locf0(cluster.number, fromLast = FALSE) ~ na.locf0(cluster.number, fromLast = TRUE), 
                           .default = cluster.number)) 
