####################################################################################################
####################################################################################################
## Script name: Treedata (mostly dominant) plots?read.
##
## Purpose of script:
##
## Date created: 2025-01-22
##
##
########################################## SYSTEM OPTIONS ##########################################
options(scipen = 6, digits = 7) 
memory.limit(30000000) 

  
########################################## DEPENDENCIES ############################################
# Packages
library(tidyverse)
library(magrittr)
library(beastio)

# User functions
GetRootInfo <- function(treedata){
  
  root_to_tip_distances <- adephylo::distRoot(treedata@phylo) %>%
    as_tibble(rownames = 'label')
  
  root_height <- root_to_tip_distances %>%
    mutate(tipdate = str_extract(label, '\\d{4}-\\d{2}-\\d{2}') %>%
             ymd()) %>%
    slice_max(tipdate) %>%
    select(-label) %>%
    mutate(tmrca = date(tipdate - dyears(value)))
  
  root_state <- as_tibble(treedata) %>%
    filter(parent == node) %>%
    pull(host_simplifiedhost)
  
  out <- root_height %>%
    mutate(host_simplifiedhost = root_state) %>%
    select(-c(value, tipdate))
  
  return(out)
}


PlotTMRCA <- function(dataframe){
  plot <- ggplot(dataframe) + 
    geom_density(aes(x = tmrca, fill = host_simplifiedhost, colour = host_simplifiedhost), alpha= 0.7, position = stack) + 
    
    scale_x_date() + 
    scale_y_continuous()  +
    
    global_theme()
  
  return(plot)
}


PlotPhyloGeo <- function(tree_tbl, basemap){
  # Format Nodes
  nodes <- tree_tbl %>%
    dplyr::select(node,height, location1, location2, host_simplifiedhost) %>%
    mutate(year = decimal_date(ymd(most_recent_date)) - as.numeric(height)) %>%
    st_as_sf(coords = c( 'location2', 'location1'), 
             crs = 4326)
  
  # Format Arrows
  arrows <-  tree_tbl %>%
    dplyr::select(node, parent, host_simplifiedhost) %>%
    left_join(mcc_tree_tbl %>%
                dplyr::select(node, location1, location2)) %>%
    left_join(mcc_tree_tbl %>%
                dplyr::select(node, location1, location2),
              by = join_by(parent== node)) %>%
    mutate(across(starts_with('location'), .fns = ~ as.numeric(.x))) %>%
    mutate(location1.y = case_when(location1.y == location1.x ~ location1.y + 0.00000001, 
                                   .default = location1.y),
           location2.y = case_when(location2.y == location2.x ~ location2.y + 0.00000001,
                                   .default = location2.y))
  
  # Plot in GGplot
  plot <- ggplot(basemap) +
    geom_sf() +
    
    # Plot Branches
    geom_curve(data = arrows, aes(x = location2.x, y = location1.x,
                                  xend = location2.y, yend = location1.y,
                                  colou r= host_simplifiedhost),
               lwd = 0.3,
               curvature = 0.2) + 
    
    # Plot nodes
    geom_sf(data = nodes, shape = 21 , size = 2, aes(fill = host_simplifiedhost))+
    
    # Set graphical scales
    scale_fill_manual(values = host_colours) + 
    scale_colour_manual(values = host_colours) +
    
    #coord_sf(ylim = c(35,60), xlim = c(-8, 33), expand = FALSE) +
    theme_void(base_size = 18) + 
    global_theme() + 
    
    
    theme(legend.background = element_rect(colour="white", fill="white"),
          legend.position =c(.7,.1),
          legend.title = element_text( vjust = 1),
          legend.direction="horizontal",
          legend.key.width  = unit(3, "lines"),
          legend.key.height = unit(1, "lines")) 
  
  return(plot)
}
############################################## DATA ################################################

# read in cross-continental tree files
treefiles <- c(list.files('./2024Aug18/reassortant_subsampled_outputs/traits_1000trees',
                        pattern = 'ha',
                        full.names = TRUE)[-10] ,
               './2024Sept16/reassortant_subsampled_outputs/traits_1000trees/ha_43112113_subsampled_traits_1000.trees')

mcc_treefiles <- c(list.files('./2024Aug18/reassortant_subsampled_outputs/traits_mcc',
                              pattern = 'ha',
                              full.names = TRUE)[-10] ,
                   './2024Sept16/reassortant_subsampled_outputs/traits_mcc/ha_43112113_subsampled_traits_mcc.tree')


mcc_trees <- lapply(mcc_treefiles, read.beast)


# Note the core requirement - this will take a long time to run in series.
df_list <- mclapply(treefiles, function(x) x %>% 
           read.beast(.) %>% 
           lapply(., GetRootInfo) %>%
           bind_rows(),
           mc.cores = 13)

names(df_list) <- gsub('.*ha_|_subsampled.*', '' ,treefiles[[i]])

host_tmrca <- df_list %>%
  bind_rows(df_list, .id = 'cluster_profile')

write_csv(host_tmrca, './2025Jan06/clusterprofile_host_tmrca.csv')
host_tmrca <- read_csv('./2025Jan06/clusterprofile_host_tmrca.csv')
############################################## MAIN ################################################

# Filter the desired reassortants


# Plot the density of TMRCA, fill by origin  host


    

ggplot(df_list[[13]]) + 
  geom_density(aes(x = tmrca, fill = host_simplifiedhost, colour = host_simplifiedhost), alpha = 0.7, position = 'stack')

# Phylogeography by colour
map <- ne_countries(scale = "medium", returnclass = "sf")




############################################## WRITE ###############################################


# For each reassortant, we need to most recent tipdate


# Infer the root age as a date

#



############################################## END #################################################
####################################################################################################
####################################################################################################