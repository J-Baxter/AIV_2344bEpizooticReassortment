test <- HostPersistence2(reassortant_trees[[10]], region_tree = F)


test %>% pivot_wider(names_from = host_simplifiedhost, values_from = starts_with('persistence'))



tbl_tree <- as_tibble(reassortant_trees[[10]])
tree <- reassortant_trees[[10]]@phylo
reassortants <- str_split_i(tree$tip.label, '\\|', 6) 

tbl_tree %<>%
  # Format numeric vectors
  mutate(across(c(branch.length,
                  height,
                  height_median,
                  length,
                  length_median),
                .fns = ~ as.numeric(.x))) %>%
  mutate(branch.length = if_else(branch.length < 0, 0.0001, branch.length))


################################### DTA for region trees ################################### 
anc_state <- AncestralStateReconstruction(tree,
                                            reassortants)


################################### path to ancestor ################################### 
reassortant_mrca <- tbl_tree %>%
  
  # tip label reassortant profiles
  mutate(cluster_profile = str_extract(label, "(\\d+_)+[^.|]*") ) %>%
  
  # Get ancestral estimates for intermediate nodes
 # rows_update(anc_state) %>%
  
  # inner join to get parental identities
  left_join(select(., parent_profile = cluster_profile,
                   parent_height_median = height_median, node),
            by = c("parent" = "node")) %>%
  
  # filter 'change' nodes (ie, nodes that differ in value to their reassortant)
  filter(cluster_profile != parent_profile) %>%
  
  # In some cases (major profiles), reassortants are not monophyletic due to incomplete sampling.
  # we therefore assume the oldest included sequence represents the 'true' root.
  group_by(cluster_profile) %>%
  slice_max(height_median, n =1) %>%
  ungroup() %>%
  
  # select key cols
  select(ends_with('profile'),
         ends_with('height_median')) 



tbl_pathtoancestor <- tbl_tree %>%
  select(node, 
         contains('height'), 
         label,
         parent) %>%
  
  # Extract reassortant types and sort by node order
  mutate(cluster_profile = str_split_i(label, '\\|', 6)) %>%
  rename(tip = node) %>%
  
  # include values of interest (node (tip) ID, height, host and reassortant)
  select(any_of(c('tip',
                  'cluster_profile'))) %>%
  #arrange(node) %>%
  
  # For each node (tip), extract the node-path back to the root
  rowwise() %>%
  mutate(path = list(bind_rows(ancestor(.data = tbl_tree, .node = tip), 
                               itself(.data = tbl_tree, .node = tip)))) %>% 
  as_tibble() %>% # cancel out rowwise()
  
  drop_na(cluster_profile) %>%
  
  unnest(path) %>%
  
  select(c(tip, 
           label,
           node, 
           cluster_profile,
           height,
           height_0.95_HPD,
           branch.length,
           host_simplifiedhost,
           host_simplifiedhost.prob)) %>%

      left_join(., anc_state %>% rename(cluster_number = cluster_profile)) 
      ##################################
      # Update with ancestral states if region tree
      #  rows_update(.,
      #             anc_state) %>%
      #   filter(cluster_profile != 'NA')
      ##################################



test_2 <- tbl_pathtoancestor %>%
  
  # copy tip cluster across into column
  mutate(cluster_number = case_when(!is.na(label) ~ cluster_profile,
                                    .default = cluster_number)) %>%
  
  # If there are missing states, fill from top down
  group_by(tip) %>%
  fill(cluster_number, .direction = 'down') %>%
  
  # keep only rows that are the same cluster as the tip
  mutate(cluster_number = case_when(cluster_number == cluster_profile ~ cluster_number,
                                    .default = NA_character_)) %>%
  drop_na(cluster_number) %>%
  select(c(tip, 
           label,
           node, 
           cluster_number, 
           height,
           height_0.95_HPD,
           branch.length,
           host_simplifiedhost,
           host_simplifiedhost.prob)) %>%
  group_by(tip) %>%

  arrange(height, .by_group = TRUE) %>%
  mutate(host_change = cumsum(host_simplifiedhost != lag(host_simplifiedhost, 
                                                         def = first(host_simplifiedhost)))) %>%
  # Calculate total persistence for each pathway
 mutate(persistence_total = sum(branch.length, na.rm = TRUE)) %>%
  #mutate(total_persistence = max(height) - min(height)) %>%

  
  # Calculate host-stratified persistence for each pathway
  group_by(host_change, .add = TRUE) %>%
  
 # mutate(persistence_host = sum(branch.length, na.rm = TRUE)) %>%
  summarise(persistence_host = sum(branch.length, na.rm = TRUE),
    #persistence_host = max(height) - min(height),
            persistence_total = mean(persistence_total),
            cluster_number = unique(cluster_number),
            host_simplified_host = unique(host_simplifiedhost))


test_3 <- test_2 %>%
  ungroup() %>%
  
  # Calculate proportions within 'lineages'
  rowwise() %>%
  mutate(persistence_host_prop = persistence_host/persistence_total) %>%
  as_tibble() %>%
  
  group_by(cluster_number) %>%
  mutate(persistence_total = mean(persistence_total)) %>%
  
  # Calculate medians from proportions
  group_by(cluster_number, host_simplified_host) %>%
  summarise(across(starts_with('persistence'), .fns = ~ median(.x), .names = '{.col}_median')) %>%
  rename(cluster_profile = cluster_number)




check <- expand_grid(cluster_profile = unique(test_2$cluster_number),
            host_simplified_host = unique(test_2$host_simplified_host)) %>%
  drop_na() %>%
  mutate(persistence_host_prop_median = NaN,
         persistence_host_median = NaN) %>%
  rows_patch(test_3,
             by = c('cluster_profile', 'host_simplified_host') ,
             unmatched = 'ignore') %>%
  replace_na(list(persistence_host_prop_median = 0)) %>%
  
  # Adjust zeros
  mutate(persistence_host_prop_median = if_else(persistence_host_prop_median <= 0, 1*10**(-6), persistence_host_prop_median)) %>%
  group_by(cluster_profile) %>%
  
  # Scale so that total sums to 1
  mutate(persistence_host_prop_median_normalised =  persistence_host_prop_median/sum(persistence_host_prop_median)) %>%
  ungroup() %>%
  
  # Format variables
  filter(host_simplified_host %in% c('anseriformes-wild', 'galliformes-domestic', 'charadriiformes-wild')) %>%
  select(-persistence_host_prop_median) %>%
  pivot_wider(names_from = host_simplified_host, 
              values_from = persistence_host_prop_median_normalised,
              names_prefix = 'scaledpersistence_') 


test_3 <- test_2 %>%
  ungroup() %>%
  group_by(cluster_number, host_simplified_host) %>%
  summarise(persistence_host_median = median(persistence_host)) %>%
  mutate(persistence_host_median = if_else(persistence_host_median == 0, 1*10**(-6), persistence_host_median)) %>%
  mutate(persistence_host_median_normalised =  persistence_host_median/sum(persistence_host_median)) %>%
  ungroup() %>%
  filter(host_simplified_host %in% c('anseriformes-wild', 'galliformes-domestic', 'charadriiformes-wild')) 







