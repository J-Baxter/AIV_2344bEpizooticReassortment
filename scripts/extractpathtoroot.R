####################################################################################################
####################################################################################################
# Extract the time spent in each host state for AIV reassortants


########################################## DEPENDENCIES ############################################
library(treeio)
library(ape)
library(tidyverse)

# Function infers tip-to-root paths for each tip, and the 'time' (time-scaled branch lengths)
# spent in each reassortant state 
HostPersistence <- function(treedata){
  
  # Set dependencies
  require(ape)
  require(treeio)
  stopifnot('Input object must be a treedata object' = class(treedata) == 'treedata')
  
  # Format treedata object as tibble, format quantitative traits as 
  # numeric vectors
  tree_tibble <- as_tibble(treedata) %>%
    mutate(across(c(branch.length,
                    height,
                    height_median,
                    length,
                    length_median),
                  .fns = ~ as.numeric(.x))) #%>%
  #{if(!'cluster_profile' %in% colnames(.)) mutate(.,cluster_profile = str_extract(label, '\\d_\\d_\\d_\\d_\\d_\\d_\\d_\\d'))} 
  
  # Extract most recent tipdate and convert to decimal 
  most_recent_tipdate <- tree_tibble %>%
    filter(height == 0) %>%
    pull(label) %>%
    str_extract(.,  "(?<=\\|)\\d{4}(?![[:lower:]]).*$") %>%
    ymd() %>%
    decimal_date()
  
  # Extract the root ID
  root_id <- tree_tibble %>%
    slice_max(., height) %>%
    pull(node)
  #select(any_of(c('cluster_profile', 'node'))) %>%
  #rename(cluster_root = node)
  
  # Infer required groupings for pipeline - reassortant should be used where applicable
  if('cluster_profile' %in% colnames(tree_tibble)){
    groups <- list(
      c('tip', 
        'cluster_profile'),
      c('tip', 
        'cluster_profile', 
        'host_simplifiedhost',
        'host_change'),
      c('cluster_profile',
        'host_simplifiedhost')
    )
  }else{
    groups <- list(
      'tip',
      c('tip', 
        'host_simplifiedhost',
        'host_change'),
      'host_simplifiedhost'
    )
  }
  
  # 
  out <- tree_tibble %>%
    # start path with only tips (ie labelled nodes)
    filter(!is.na(label))  %>%
    rename(tip = node) %>%
    
    # include values of interst (node (tip) ID, height, host and reassortant)
    select(tip) %>%
    
    # For each node (tip), extract the node-path back to the root
    rowwise() %>%
    mutate(path = list(tibble(node = ape::nodepath(treedata@phylo, tip, root_id)))) %>%
    as_tibble() %>% # cancel out rowwise()
    unnest(path) %>%
    
    # For each tip-path, join the median height, branch length and host
    left_join(x = . ,
              y = tree_tibble %>%
                select(any_of(c('node',
                                'parent',
                                'height_median',
                                'branch.length',
                                'host_simplifiedhost',
                                'reassortant'))),
              by = join_by(node == node)) %>%
    
    # calculate decimal dates from node heights
    #mutate(date = most_recent_tipdate - height_median) %>%
    
    # for each tip, determine whether the host state changes at each node along path to root
    group_by(groups[[1]]) %>%
    mutate(host_change = cumsum(host_simplifiedhost != lag(host_simplifiedhost, 
                                                           def = first(host_simplifiedhost)))) %>%
    ungroup() %>%
    
    # sum branch lengths for contiguous periods in the same host
    # ie. we do not combine time spent in the same host, interdispersed by a different host class
    summarise(persistence_host = sum(branch.length, na.rm = TRUE), 
              .by = groups[[2]]) %>%
    
    
    # calculate total, maximum and median persistence times for each host class for each reassortant 
    summarise(persistence_host_median = median(persistence_host),
              persistence_host_max = max(persistence_host),
              persistence_host_sum = sum(persistence_host),
              .by = groups[[3]])
  
  
  return(out)
  
}



############################################# DATA ################################################
# Import MCC trees (reassortant)
reassortant_trees <- list.files('./2024Aug18/reassortant_subsampled_outputs/traits_mcc',
                                full.names = T) %>%
  lapply(., read.beast)

region_trees <- list.files('./2024Aug18/region_subsampled_outputs/traits_mcc',
                           full.names = T) %>%
  lapply(., read.beast)


#m <- region_trees[[1]]  %>%
  #as_tibble() %>%
  #filter(height == 0) %>%
  #pull(label) %>%
  #str_extract(.,  "(?<=\\|)\\d{4}(?![[:lower:]]).*$") %>%
  #ymd() %>%
  #decimal_date() %>% 
  #max(na.rm = T)


GetReassortantMRCAs <- function(treedata){
  # Set dependencies
  require(ape)
  require(treeio)
  stopifnot('Input object must be a treedata object' = class(treedata) == 'treedata')
  
  # Format treedata object as tibble, format quantitative traits as 
  # numeric vectors
  tree_tibble <- as_tibble(treedata) %>%
    mutate(across(c(branch.length,
                    height,
                    height_median,
                    length,
                    length_median),
                  .fns = ~ as.numeric(.x))) 
  
  # format phylo object
  tree <- treedata@phylo
  
  tree$edge.length <-ifelse(tree$edge.length < 0,
                               0.0001, 
                               tree$edge.length)
  
  # Extract the cluster profiles from tipnames
  cluster_profiles <- str_split_i(tree$tip.label, '\\|', 6) 
  
  # Maximum likelihood ancestral state reconstruction for cluster_profile
  anc_state <- ace(cluster_profiles,
                  tree,
                  model = 'ER', 
                  type = 'discrete')
  
  # extract likelihood for each state and each node
  lik_anc <- anc_state$lik.anc
  
  # infer most likely state at each node 
  anc_cluster <- apply(lik_anc,
                       1,
                       function(x) names(x)[which.max(x)]) %>%
    
    #format output
    enframe(., 
            name = 'node', 
            value = 'cluster_profile') %>%
    mutate(node = as.integer(node))
  
  
  test_reassortanttmcra <- tree_tibble %>%
    select(node, 
           contains('height'), 
           label,
           parent) %>%
    
    # Extract reassortant types and sort by node order
    mutate(cluster_profile = str_split_i(label, '\\|', 6)) %>%
    arrange(node) %>%
    
    # Update with ancestral states
    rows_update(.,
                anc_cluster) %>%
    
    # Group by cluster and identify eldest node
    group_by(cluster_profile) %>%
    
    # If only one example in group, then take height of parent node instead of tip date
    mutate(node = case_when(n() == 1 ~ parent, 
                            TRUE ~ node)) %>%
    
    mutate(across(starts_with('height'),
                  .fns = ~ case_when(n() == 1 ~ NA,
                                     TRUE~ .x))) %>%
    rows_patch(tree_tibble %>% 
                 select(node,
                        contains('height')), 
               unmatched = 'ignore') %>%
    
    # Sample tallest (ie oldest) node for each reassortant
    slice_max(height_median) %>%
    ungroup() %>%
    select(cluster_profile, node) 
  
  
  return(test_reassortanttmcra)
}

library(parallel)

ncpus = parallel::detectCores()-2
cl = makeCluster(ncpus)


reassortant_mrca <- mclapply(region_trees,  
                           function (x){return(tryCatch(GetReassortantMRCAs(x), 
                                                        error=function(e) NULL))})

reassortant_mrca_nodoms <- reassortant_mrca %>%
  lapply(., function(x)
    return(tryCatch(x %>% 
                      filter(!cluster_profile %in% c(
                        "3_2_3_1_3_2_1_2",
                        "2_1_1_1_1_1_1_1",
                        "1_1_1_1_1_1_1_1",
                        "2_1_2_1_1_1_1_1",
                        "1_1_2_1_1_1_1_1",
                        "1_6_2_1_1_1_1_1",
                        "2_6_1_1_6_1_1_1",
                        "2_1_6_1_1_4_1_1",
                        "7_1_5_2_1_3_1_2",
                        "4_3_1_1_2_1_1_3",
                        "5_1_1_1_2_1_1_3",
                        "5_4_9_1_2_1_1_1",
                        'NA')), 
                    error=function(e) NULL)))


stopCluster(cl)
############################################## RUN ################################################
persistence_dataframe <- lapply(reassortant_trees, HostPersistence) %>%
  
  # Set names for each tree extraction
  set_names(list.files('./2024Aug18/reassortant_subsampled_outputs/traits_mcc') %>% 
              gsub('_subsampled.*', '', .)) %>%
  
  # concatenate to a single dataframe and separate segment and reassortant
  bind_rows(., .id = 'tree') %>%
  separate_wider_delim(tree, '_', names = c('segment', 'reassortant'))


