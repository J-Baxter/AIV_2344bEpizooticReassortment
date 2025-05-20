####################################################################################################
####################################################################################################
## Script name:
##
## Purpose of script:
##
## Date created: 2025-05-14
##
##
########################################## SYSTEM OPTIONS ##########################################
options(scipen = 6, digits = 7) 
memory.limit(30000000) 

  
########################################## DEPENDENCIES ############################################
# Packages
library(tidyverse)
library(magrittr)
library(tidygraph)
library(ggraph)

# User functions


############################################## DATA ################################################

reassortant_ancestral_changes <- read_csv('./reassortant_ancestral_changes.csv')


############################################## MAIN ################################################
# Plot network using tidygraph
my_graph <- reassortant_ancestral_changes %>% 
  select(ends_with('label'), 
         segments_changed, 
         time_since_parent) %>% 
  drop_na() %>%
  relocate(parent_label,
           cluster_label) %>%
  
  # Make graph
  as_tbl_graph(.) %>%
  
  # Add node data
  activate(nodes) %>%
  left_join(reassortant_ancestral_changes %>% 
              select(name = cluster_label, cluster_class)) %>% mutate(importance = centrality_degree())

as_tibble(my_graph) %>% 
  summarise(as_tibble_row(quantile(importance)), .by = cluster_class)

ggraph(my_graph %>%
         mutate(component = group_components()) %>%
         filter(component == which.max(as.numeric(table(component)))),  layout = "fr", weights = segments_changed) +
  geom_edge_link() +
  geom_node_point(aes(color = cluster_class, size=10 )) +
  geom_node_label(aes(label = ifelse(cluster_class == 'major', name, '')), repel = TRUE) +
  scale_colour_brewer(palette = 'Set1') +
  theme_void() + facet_wrap(~cluster_region)

degree(d= my_graph, mode = 'out', loops = F)
############################################## WRITE ###############################################




############################################## END #################################################
####################################################################################################
####################################################################################################