####################################################################################################
####################################################################################################
# Extract the time spent in each host state for AIV reassortants


########################################## DEPENDENCIES ############################################
library(treeio)
library(ape)
library(tidyverse)
library(parallel)


# Helper function extracted from tidytree package 
itself <- function(.data, .node) {
  if (is.numeric(.node)) {
    i <- which(.data$node == .node)
  } else {
    i <- which(.data$label == .node)
  }
  
  ## .data[which(.data$node == .node | .data$label == .node), ]
  return(.data[i, ])
}


# Wrapper function for ancestral state reconstruction
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
    mutate(node = as.integer(node),
           cluster_profile = gsub('NA', NA_character_, cluster_profile))
  
  return(anc_state)
}


# Function to calculate and summarise host-persistence per reassortant
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
    {
      if (region_tree){
        left_join(., anc_state %>% rename(cluster_number = cluster_profile)) 
        ##################################
        # Update with ancestral states if region tree
      #  rows_update(.,
       #             anc_state) %>%
       #   filter(cluster_profile != 'NA')
        ##################################
        } else {
          .
        }
      }
  
  ################################### summarise host persistence ################################### 
  
  cols <- c('cluster_number' = reassortants[1])
  tbl_pathtoancestor %<>% 
    add_column(!!!cols[setdiff(names(cols), names(.))])
  
  # For each tip, include only nodes that match the given tip cluster
  out <- tbl_pathtoancestor %>%
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
    group_by(host_change, .add = TRUE) %>%
    #mutate(persistence_host = sum(branch.length, na.rm = TRUE)) %>%
    mutate(persistence_host = max(height) - min(height)) %>%
    
    ungroup() %>%
    
    group_by(cluster_number, host_simplifiedhost) %>%
    summarise(persistence_host_median = median(persistence_host)) %>%
    mutate(persistence_host_median_normalised =  persistence_host_median/sum(persistence_host_median)) %>%
    ungroup()
    
    
    
    
   # tbl_pathtoancestor %>%
   # group_by(tip) %>%
    
   # filter(cluster_profile %in% str_split_i(label, '\\|', 6)) %>% ########### need to check that tip == 
    
    # If only one node/tip in group, then add parent node
   # group_modify(~ {
    #  if (nrow(.x) == 1){
    #    bind_rows(.x, parent(tbl_tree, .node = .x$node)) %>%
     #     fill(cluster_profile, .direction = 'downup')
    #  }else{
    #    .x
     # }
    #}) %>%
    #ungroup() %>%
    
    #select(any_of(c('tip', 
    #                'cluster_profile',
    #                'node',
     #               'parent',
     #               'height_median',
     #               'branch.length',
     #               'host_simplifiedhost',
     #               'reassortant'))) %>%
    #group_by(cluster_profile, 
       #      tip) %>%
    
    #mutate(host_change = cumsum(host_simplifiedhost != lag(host_simplifiedhost, 
    #                                                       def = first(host_simplifiedhost)))) %>%
    #ungroup() %>%
    
    # sum branch lengths for contiguous periods in the same host
    # ie. we do not combine time spent in the same host, interdispersed by a different host class
   # summarise(persistence_host = sum(branch.length, na.rm = TRUE), 
     #         .by = c(cluster_profile,
       #               tip, 
      #                host_simplifiedhost,
        #              host_change)) %>%
    
    
    # calculate total, maximum and median persistence times for each host class for each reassortant 
   # summarise(persistence_host_median = median(persistence_host),
    #          persistence_host_max = max(persistence_host),
    #          persistence_host_sum = sum(persistence_host),
     #         .by = c(cluster_profile,
      #                host_simplifiedhost))
  
  return(out)
  
  
}


############################################# DATA ################################################
# Import MCC trees (reassortant)
reassortant_trees <- list.files('./2024Aug18/reassortant_subsampled_outputs/traits_mcc',
                                pattern = 'ha_',
                                full.names = T) %>%
  lapply(., read.beast)

region_trees <- list.files('./2024Aug18/region_subsampled_outputs/traits_mcc',
                           full.names = T,
                           pattern = 'ha_') %>%
  lapply(., read.beast)


############################################## RUN ################################################


tbl_moderateminorpersistence <- mclapply(region_trees, HostPersistence2, region_tree = TRUE, mc.cores = 12) %>%
  # Set names for each tree extraction
  set_names(list.files('./2024Aug18/region_subsampled_outputs/traits_mcc', pattern = 'ha_') %>% 
              gsub('_subsampled.*', '', .)) %>%
  
  # concatenate to a single dataframe and separate segment and reassortant
  bind_rows(., .id = 'tree') %>%
  separate_wider_delim(tree, '_', names = c('segment', 'region'))%>% 
  rename(cluster_profile = cluster_number) %>%
  
  # remove dominant cluster profiles
  filter(!cluster_profile %in% c("3_2_3_1_3_2_1_2",
                                "2_1_1_1_1_1_1_1",
                                "1_1_1_1_1_1_1_1",
                                "1_1_1_1_1_1_1_1A",
                                "2_1_2_1_1_1_1_1",
                                "1_1_2_1_1_1_1_1",
                                "1_6_2_1_1_1_1_1",
                                "1_1_4_1_4_1_1_4",
                                "2_6_1_1_6_1_1_1",
                                "2_1_6_1_1_4_1_1",
                                "7_1_5_2_1_3_1_2",
                                "4_3_1_1_2_1_1_3",
                                "5_1_1_1_2_1_1_3",
                                "5_4_9_1_2_1_1_1"))

test4 <- tbl_moderateminorpersistence %>%
  group_by(cluster_profile, host_simplifiedhost) %>%
  slice_min(persistence_host_median)

tbl_major_persistence <- mclapply(reassortant_trees, HostPersistence2, mc.cores = 12
                                      ) %>%
  # Set names for each tree extraction
  set_names(list.files('./2024Aug18/reassortant_subsampled_outputs/traits_mcc',
                       , pattern = 'ha_') %>% 
              gsub('_subsampled.*', '', .)) %>%
  
  # concatenate to a single dataframe and separate segment and reassortant
  bind_rows(., .id = 'tree') %>%
  separate_wider_delim(tree, '_', names = c('segment', 'reassortant'))


stopCluster(cl)


############################################## WRITE ################################################
tbl_combined <- expand_grid(cluster_profile = unique(meta$cluster_profile),
                            host_simplified_host = c("galliformes-domestic",
                                                     "anseriformes-wild",
                                                     "charadriiformes-wild")) %>%
  drop_na() %>%
  mutate(scaled_persistence = NaN) %>%
  rows_patch(tbl_moderateminorpersistence %>% select(cluster_profile, 
                                                      host_simplified_host = host_simplifiedhost,
                                                     scaled_persistence = persistence_host_median_normalised),
              by = c('cluster_profile', 'host_simplified_host'),
             unmatched = 'ignore') %>%
  rows_patch(tbl_major_persistence %>% select(cluster_profile = cluster_number, 
                                                     host_simplified_host = host_simplifiedhost,
                                                     scaled_persistence = persistence_host_median_normalised),
             by = c('cluster_profile', 'host_simplified_host'),
             unmatched = 'ignore') %>%
  replace_na(list(scaled_persistence = 0)) %>%
  pivot_wider(names_from = host_simplified_host,
              values_from = scaled_persistence)


  select(-region) %>%
  bind_rows(tbl_major_persistence %>% select(-reassortant)) %>%
  #mutate(cluster_profile = coalesce(cluster_profile, cluster_number)) %>%
  #dplyr::select(-cluster_number) %>%
  mutate(across(starts_with('persistence'), .fns = ~ ifelse(.x < 0, 0, .x)))

write_csv(tbl_combined, './2024Aug18/treedata_extractions/reassortant_stratifiedpersistence.csv')

############################################## END #################################################
####################################################################################################
####################################################################################################