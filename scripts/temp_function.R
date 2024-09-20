AncestralStateReconstruction <- function(tree, 
                                         trait.matrix,
                                         model = 'ER',
                                         type = 'discrete'){
  require(ape)
  
  # format phylo object
  tree$edge.length <-ifelse(tree$edge.length < 0,
                            0.0001, 
                            tree$edge.length)
  
  
  anc_state <- ape::ace(trait.matrix,
                        tree,
                        model = model, 
                        type = type)
  
  return(anc_state)
}




tree_tibble <- as_tibble(region_trees[[2]]) %>%
  mutate(across(c(branch.length,
                  height,
                  height_median,
                  length,
                  length_median),
                .fns = ~ as.numeric(.x))) 



##################################
# Maximum likelihood ancestral state reconstruction for cluster_profile

if{region_tree}{
  anc_state <- AncestralStateReconstruction()
  
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
  

  # For each tip, include only nodes that match the given tip cluster
  group_by(tip) %>%
  filter(cluster_profile %in% str_split_i(label, '\\|', 6)) %>% # need to check that tip == 
  
  # If only one node/tip in group, then add parent node
  group_modify(~ {
    if (nrow(.x) == 1){
      bind_rows(.x, parent(tree_tibble, .node = .x$node)) %>%
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
  
  
