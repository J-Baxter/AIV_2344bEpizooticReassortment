test <- HostPersistence2(region_trees[[4]], region_tree = TRUE)


test %>% pivot_wider(names_from = host_simplifiedhost, values_from = starts_with('persistence'))



tbl_tree <- as_tibble(region_trees[[4]])
tree <- region_trees[[4]]@phylo
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
anc_state <- AncestralStateReconstruction(tree,
                                            reassortants)


################################### path to ancestor ################################### 
reassortant_mrca <- tbl_tree %>%
  
  # tip label reassortant profiles
  mutate(cluster_profile = str_extract(label, "(\\d+_)+[^.|]*") ) %>%
  
  # Get ancestral estimates for intermediate nodes
  rows_update(anc_state) %>%
  
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
  #group_modify(~{
  #if (nrow(.x) == 1){
  # bind_rows(.x, tidytree::parent(.data = tbl_tree, .node = .x$node)) %>%
  #    fill(cluster_number, .direction = 'downup')
  #  }else{
  #   .x
  # }
  # }) %>%
  arrange(height, .by_group = TRUE) %>%
  mutate(host_change = cumsum(host_simplifiedhost != lag(host_simplifiedhost, 
                                                         def = first(host_simplifiedhost)))) %>%
  # Calculate total persistence for each pathway
  mutate(total_persistence = sum(branch.length, na.rm = TRUE)) %>%
  
  # Calculate host-stratified persistence for each pathway
  group_by(host_change, .add = TRUE) %>%
 # mutate(persistence_host = sum(branch.length, na.rm = TRUE)) %>%
  summarise(persistence_host = sum(branch.length, na.rm = TRUE),
            total_persistence = mean(total_persistence),
            cluster_number = unique(cluster_number),
            host_simplified_host = unique(host_simplifiedhost))


test_2 %>%
  ungroup() %>%
  pivot_wider(names)
  group_by(cluster_number, host_simplifiedhost) %>%
  summarise(persistence_host_median = median(persistence_host),
            persistence_host_max = max(persistence_host, na.rm= TRUE),
            persistence_host_sum = sum(persistence_host, na.rm= TRUE)) %>%
  ungroup()



