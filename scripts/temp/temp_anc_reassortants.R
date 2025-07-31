# Packages
library(tidyverse)
library(magrittr)
library(RColorBrewer)
library(ggtree)
library(ggtreeExtra)
library(scales)
library(cowplot)
library(treeio)
library(TreeTools)
library(ggnewscale)
library(rnaturalearth)
library(rnaturalearthdata)
library(ggstream)
# User functions



combined_data <- read_csv('./2024Aug18/treedata_extractions/2024-09-20_combined_data.csv')
summary_data <- read_csv('./2024Aug18/treedata_extractions/summary_reassortant_metadata_20240904.csv') %>%
  select(-c(cluster_label,
            clade)) 

meta <- read_csv('./2024-09-09_meta.csv') 

new_tree <- read.beast('./2025Feb26/globalsubsample/ha_global_SRD06_relaxLn_constant_mcc.tree')
final_clusters <- read_csv('./2025Jun10/2025Jun10_reasssortantclusters.csv')

# reassortant MRCA (will run for 10-20 mins)
traits <- new_tree@phylo$tip.label %>%
  str_extract(., "(\\d+_)+[^.|]*") 

new_tree@phylo$edge.length <- ifelse(new_tree@phylo$edge.length < 0, 0, new_tree@phylo$edge.length)

ace_test <- ace(traits,
                new_tree@phylo,
                type = 'discrete',
                method = 'ML',
                model = 'ER')

#ace_test$lik.anc


# extract likelihood for each state and each node
anc_lik <- ace_test$lik.anc

# infer most likely state at each node and 
anc_state <- apply(anc_lik,
                   1,
                   function(x) names(x)[which.max(x)]) %>%
#format output
  enframe(., 
          name = 'node', 
          value = 'cluster_profile') %>%
  mutate(node = as.integer(node))



test <- as_tibble(new_tree) %>%
  
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
  
  # Count the number of segments changed 
  rowwise() %>%
  mutate(segments_changed = sum(as.integer(str_split_fixed(cluster_profile, '_', 8)) != as.integer(str_split_fixed(parent_profile, '_', 8)), na.rm = TRUE)) %>%
  as_tibble() %>%
  
  # In some cases (major profiles), reassortants are not monophyletic due to incomplete sampling.
  # we therefore assume the oldest included sequence represents the 'true' root.
  group_by(cluster_profile) %>%
  slice_max(height_median, n =1) %>%
  ungroup() %>%
  
  # select key cols
  select(ends_with('profile'),
         ends_with('height_median'),
         segments_changed) %>%
  
  # map present and ancestral class
  left_join(final_clusters)  %>%
  rename(cluster_class = .cluster) %>%
  left_join(select(., 
                   parent_class = cluster_class,
                   parent_label = cluster_label,
                   cluster_profile), 
            by = c("parent_profile" = "cluster_profile")) %>%
  
  mutate(time_since_parent = parent_height_median-height_median)

write_csv(test,'./reassortant_ancestral_changes.csv')

