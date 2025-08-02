################################################################################
## Script Name:       Evolutionary time proportions  
## Purpose:            Extract the time spent in each host state for AIV 
##                     reassortants
## Author:             James Baxter
## Date Created:       2025-08-02
################################################################################

############################### SYSTEM OPTIONS #################################
options(
  scipen = 6,     # Avoid scientific notation
  digits = 7      # Set precision for numerical display
)
memory.limit(30000000)

############################### DEPENDENCIES ###################################
# Load required libraries
library(tidyverse)
library(magrittr)
library(treeio)
library(ape)
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
  
  ################################### summarise persistence ################################### 
  
  cols <- c('cluster_number' = reassortants[1])
  tbl_pathtoancestor %<>% 
    add_column(!!!cols[setdiff(names(cols), names(.))])
  
  persistence <- tbl_pathtoancestor %>%
    
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
    group_by(cluster_number) %>%
    summarise(persistence = max(height) - min(height)) %>%
    rename(cluster_profile = cluster_number)
  
  
  ################################### summarise host proportions ################################### 
  # For each tip, include only nodes that match the given tip cluster
  out <- tbl_pathtoancestor %>%
    
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
    group_by(cluster_number) %>%
    distinct(node, .keep_all = T) %>%
    
    mutate(path_total = sum(branch.length, na.rm = T)) %>%
    group_by(cluster_number, host_simplifiedhost) %>%
    summarise(path_total_host = sum(branch.length, na.rm = T),
              path_total = mean(path_total))  %>%
    rename(cluster_profile = cluster_number) %>%
    mutate(path_prop_host = path_total_host/path_total) %>%
    ungroup() %>%
    
    select(cluster_profile,
           host_simplifiedhost,
           path_prop_host) %>%
    pivot_wider(names_from = host_simplifiedhost,
                values_from = path_prop_host,
                names_prefix = 'path_prop_') %>%
    mutate(across(starts_with('path'), .fns = ~replace_na(.x, 0))) %>%
    left_join(persistence)
  
  return(out)
  
  
}



################################### DATA #######################################
# Read and inspect data
# Import MCC trees (reassortant)
reassortant_trees <- list.files('./2024Aug18/reassortant_subsampled_outputs/traits_mcc',
                                pattern = 'ha_',
                                full.names = T) %>%
  lapply(., read.beast)

region_trees <- list.files('./2024Aug18/region_subsampled_outputs/traits_mcc',
                           full.names = T,
                           pattern = 'ha_') %>%
  lapply(., read.beast)


################################### MAIN #######################################
# Main analysis or transformation steps
tbl_moderateminorpersistence <- mclapply(region_trees, 
                                         HostPersistence2,
                                         region_tree = TRUE, 
                                         mc.cores = 12) %>%
  # Set names for each tree extraction
  set_names(list.files('./2024Aug18/region_subsampled_outputs/traits_mcc', 
                       pattern = 'ha_') %>% 
              gsub('_subsampled.*', '', .)) %>%
  
  # concatenate to a single dataframe and separate segment and reassortant
  bind_rows(., .id = 'tree') %>%
  separate_wider_delim(tree, '_', names = c('segment', 'region'))%>% 
  #rename(cluster_profile = cluster_number) %>%
  
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
                                 "5_4_9_1_2_1_1_1")) %>%
  
  # one sequence in south america (or fewer)
  filter(!(cluster_profile == '6_4_8_1_2_1_1_1' & region == 'southamerica')) 

#test4 <- tbl_moderateminorpersistence %>%
#group_by(cluster_profile, host_simplifiedhost) %>%
#slice_min(persistence_host_median)

tbl_major_persistence <- mclapply(reassortant_trees, HostPersistence2, mc.cores = 12
) %>%
  # Set names for each tree extraction
  set_names(list.files('./2024Aug18/reassortant_subsampled_outputs/traits_mcc',
                       , pattern = 'ha_') %>% 
              gsub('_subsampled.*', '', .)) %>%
  
  # concatenate to a single dataframe and separate segment and reassortant
  bind_rows(., .id = 'tree') %>%
  separate_wider_delim(tree, '_', names = c('segment', 'reassortant'))


################################### OUTPUT #####################################
# Save output files, plots, or results
tbl_combined <- bind_rows(tbl_major_persistence %>%
                            select(-reassortant),
                          tbl_moderateminorpersistence  %>% 
                            select(-region))  %>%
  mutate(across(starts_with('path'), .fns = ~replace_na(.x, 0))) %>%
  filter(cluster_profile != 'NA') %>%
  drop_na() 


tbl_interest <- tbl_combined %>%
  select(cluster_profile, 
         ends_with(c("_galliformes-domestic",
                     "_anseriformes-wild",
                     "_charadriiformes-wild")),
         persistence) 


write_csv(tbl_combined, './2025Jun10/reassortant_stratifiedpersistence.csv')

#################################### END #######################################
################################################################################