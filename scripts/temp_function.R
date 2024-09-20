AncestralStateReconstruction <- function(tree, trait.matrix, model = 'ER', type = 'discrete'){
  require(ape)
  
  # format phylo object
  tree$edge.length <-ifelse(tree$edge.length < 0,
                            0.0001, 
                            tree$edge.length)
  
  # Run discrete trait analysis
  anc <- ape::ace(trait.matrix,
                  tree,
                  model = model, 
                  type = type)
  
  # extract likelihood for each state and each node
  anc_lik <- anc$lik.anc
  
  # infer most likely state at each node and 
  anc_state <- apply(anc_lik,
                       1,
                       function(x) names(x)[which.max(x)]) %>%
    #format output
    enframe(., 
            name = 'node', 
            value = 'cluster_profile') %>%
    mutate(node = as.integer(node))
  
  return(anc_state)
}

HostPersistence2 <- function(treedata, region_tree = FALSE){
  
  ################################### dependencies ################################### 
  require(ape)
  require(tidyverse)
  require(treeio)
  require(tidytree)
  
  ################################### format input ################################### 
  
  tbl_tree <- as_tibble(treedata)
  tree <- treedata@phylo
  reassortants <- str_split_i(tree$tip.label, '\\|', 6) 
  
  tbl_tree %<>%
    # Format numeric vectors
    mutate(across(c(branch.length,
                    height,
                    height_median,
                    length,
                    length_median),
                  .fns = ~ as.numeric(.x)))
  
  
  ################################### DTA for region trees ################################### 
  if(region_tree){
    anc_state <- AncestralStateReconstruction(tree,
                                              reassortants)
  } 
  
  
  ################################### path to ancestor ################################### 
  tbl_pathtoancestor <- tbl_tree %>%
    select(node, 
           contains('height'), 
           label,
           parent) %>%
    
    # Extract reassortant types and sort by node order
    mutate(cluster_profile = str_split_i(label, '\\|', 6)) %>%
    rename(tip = node) %>%
    
    # include values of interst (node (tip) ID, height, host and reassortant)
    select(any_of(c('tip',
                    'cluster_profile'))) %>%
    #arrange(node) %>%
    
    # For each node (tip), extract the node-path back to the root
    rowwise() %>%
    mutate(path = list(bind_rows(ancestor(.data = tbl_tree, .node = tip), 
                                 itself(.data = tbl_tree, .node = tip)))) %>% 
    as_tibble() %>% # cancel out rowwise()
    
    filter(!is.na(cluster_profile)) %>%
    
    unnest(path) %>%
    {
      if (region_tree){
        ##################################
        # Update with ancestral states if region tree
        rows_update(.,
                    anc_state) %>%
          filter(cluster_profile != 'NA')
          ##################################
      } else {
        .
      }
    }
  
  ################################### summarise host persistence ################################### 
  
  # For each tip, include only nodes that match the given tip cluster
  out <- tbl_pathtoancestor %>%
    group_by(tip) %>%
    
    filter(cluster_profile %in% str_split_i(label, '\\|', 6)) %>% # need to check that tip == 
    
    # If only one node/tip in group, then add parent node
    group_modify(~ {
      if (nrow(.x) == 1){
        bind_rows(.x, parent(tbl_tree, .node = .x$node)) %>%
          fill(cluster_profile, .direction = 'downup')
      }else{
        .x
      }
    }) %>%
    ungroup() %>%
    
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
              .by = c(cluster_profile,
                      host_simplifiedhost))
  
  return(out)
  
  
}


reassort

HostPersistence2(reassortant_trees[[1]], region_tree = F)



##################################
# Maximum likelihood ancestral state reconstruction for cluster_profile

if{region_tree}{
  anc_state <- AncestralStateReconstruction()
  

  ##################################
} else{
  
}









test_reassortanttmcra <- tree_tibble %>%
  select(node, 
         contains('height'), 
         label,
         parent) %>%
  
  # Extract reassortant types and sort by node order
  mutate(cluster_profile = str_split_i(label, '\\|', 6)) %>%
  rename(tip = node) %>%
  
  # include values of interst (node (tip) ID, height, host and reassortant)
  select(any_of(c('tip',
                  'cluster_profile'))) %>%
  #arrange(node) %>%
  
  # For each node (tip), extract the node-path back to the root
  rowwise() %>%
  mutate(path = list(bind_rows(ancestor(.data = tree_tibble, .node = tip), 
                               itself(.data = tree_tibble, .node = tip)))) %>% #list(tibble(node = ape::nodepath(treedata@phylo, tip, root_id)))
  as_tibble() %>% # cancel out rowwise()

  filter(!is.na(cluster_profile)) %>%
  
  unnest(path) %>%
  {
    if (region_tree){
      ##################################
      # Update with ancestral states
      rows_update(.,
                  anc_cluster) %>%
        filter(cluster_profile != 'NA') %>%
        ##################################
    } else {
      .
    }
  } %>%
  

  
  
  
