####################################################################################################
####################################################################################################
# Extract the time spent in each host state for AIV reassortants


########################################## DEPENDENCIES ############################################
library(treeio)
library(ape)
library(tidyverse)
library(parallel)

# Function infers tip-to-root paths for each tip, and the 'time' (time-scaled branch lengths)
# spent in each reassortant state 
HostPersistence <- function(tbl_tree){
  
  # Set dependencies
  require(ape)
  require(treeio)
  require(tidyverse)
  #stopifnot('Input object must be a tbl_tree object' = class(treedata) == "tbl_tree")
  
  # Format treedata object as tibble, format quantitative traits as 
  # numeric vectors
  tbl_tree %<>%
    mutate(across(c(branch.length,
                    height,
                    height_median,
                    length,
                    length_median),
                  .fns = ~ as.numeric(.x)))
  
  #%>%
  #{if(!'cluster_profile' %in% colnames(.)) mutate(.,cluster_profile = str_extract(label, '\\d_\\d_\\d_\\d_\\d_\\d_\\d_\\d'))} 
  
  # Extract most recent tipdate and convert to decimal 
  #most_recent_tipdate <- tree_tibble %>%
    #filter(height == 0) %>%
    #pull(label) %>%
    #str_extract(.,  "(?<=\\|)\\d{4}(?![[:lower:]]).*$") %>%
    #ymd() %>%
    #decimal_date()
  
  
  out <- tbl_tree %>%
    mutate(across(c(branch.length,
                    height,
                    height_median,
                    length,
                    length_median),
                  .fns = ~ as.numeric(.x))) %>%
    # start path with only tips (ie labelled nodes)
    filter(!is.na(label))  %>%
    rename(tip = node) %>%
    mutate(cluster_profile = str_split_i(label, '\\|', 6)) %>%
    rename(tip = node) %>%
    
    # include values of interst (node (tip) ID, height, host and reassortant)
    select(any_of(c('tip',
                    'cluster_profile'))) %>%
    
    # For each node (tip), extract the node-path back to the root
    rowwise() %>%
    
    # check the following works
    mutate(path = list(bind_rows(ancestor(.data = tree_tibble, .node = tip), 
                                 itself(.data = tree_tibble, .node = tip)))) %>% #list(tibble(node = ape::nodepath(treedata@phylo, tip, root_id)))
    as_tibble() %>% # cancel out rowwise()
    filter(!is.na(cluster_profile)) %>%
    
    unnest(path) %>%
    
    select(any_of(c('tip', 
                    'cluster_profile',
                    'node',
                    'parent',
                    'height_median',
                    'branch.length',
                    'host_simplifiedhost',
                    'reassortant'))) %>%
    group_by(cluster_profile, 
             tip) %>%
    mutate(host_change = cumsum(host_simplifiedhost != lag(host_simplifiedhost, 
                                                           def = first(host_simplifiedhost)))) %>%
    ungroup() %>%
    
    # sum branch lengths for contiguous periods in the same host
    # ie. we do not combine time spent in the same host, interdispersed by a different host class
    summarise(persistence_host = sum(branch.length, na.rm = TRUE), 
              .by = c(cluster_profile,
                      tip, 
                      host_simplifiedhost,
                      host_change)) %>%
    
    
    # calculate total, maximum and median persistence times for each host class for each reassortant 
    summarise(persistence_host_median = median(persistence_host),
              persistence_host_max = max(persistence_host),
              persistence_host_sum = sum(persistence_host),
              .by = host_simplifiedhost)
  
  
  return(out)
  
}


# 
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


# From tidytree package code
itself <- function(.data, .node) {
  if (is.numeric(.node)) {
    i <- which(.data$node == .node)
  } else {
    i <- which(.data$label == .node)
  }
  
  ## .data[which(.data$node == .node | .data$label == .node), ]
  return(.data[i, ])
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


# Extract reassortant common ancestor from region trees
ncpus = detectCores()-2
cl = makeCluster(ncpus)

reassortant_mrca <- mclapply(region_trees,  
                           function (x){return(tryCatch(GetReassortantMRCAs(x), 
                                                        error=function(e) NULL))})

stopCluster(cl)

# exclude dominant (ie multi-continent) reassortants
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


# infer subtrees for each reassortant within each region

ExtractReassortantSubtrees <- function(reassortant_anc, 
                                       treedata){
  tbl_tree <- as_tibble(treedata)
  
  cluster_profile <- reassortant_anc %>%
    pull(cluster_profile)
  
  out <- reassortant_anc %>%
    rename(anc_node = node) %>%
    rowwise() %>%
    mutate(tree_tibble = list(offspring(tbl_tree, 
                                        anc_node, 
                                        self_include = T))) %>%
    unnest(tree_tibble) %>%
    select(-anc_node) %>%
    group_split(cluster_profile, .keep = FALSE) %>%
    lapply(., as.treedata) %>%
    lapply(., as_tibble) %>%
    setNames(cluster_profile)
  
  return(out)
}


n <- list.files('./2024Aug18/region_subsampled_outputs/traits_mcc',
                full.names = T) %>%
  str_split_i(., '\\/', 5) %>%
  gsub('_subsampled_traits_mcc.tree', '', .)

reassortant_subtrees <-  mapply(function(x,y){return(tryCatch(ExtractReassortantSubtrees(x,y), 
                                                            error=function(e) NULL))},
                                reassortant_mrca_nodoms,
                                region_trees,
                                SIMPLIFY = F) %>%
  set_names(n)
  


persistence_dataframe_2 <- lapply(reassortant_subtrees, HostPersistence) 

############################################## RUN ################################################
persistence_dataframe <- lapply(reassortant_trees, HostPersistence) %>%
  
  # Set names for each tree extraction
  set_names(list.files('./2024Aug18/reassortant_subsampled_outputs/traits_mcc') %>% 
              gsub('_subsampled.*', '', .)) %>%
  
  # concatenate to a single dataframe and separate segment and reassortant
  bind_rows(., .id = 'tree') %>%
  separate_wider_delim(tree, '_', names = c('segment', 'reassortant'))


