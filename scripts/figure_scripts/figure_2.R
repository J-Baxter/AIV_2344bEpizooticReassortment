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
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(cowplot)

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
  plot <- ggplot(dataframe,
                 aes(x  = tmrca,
                     colour = host_simplifiedhost,
                     #y = host_simplifiedhost,
                     fill = host_simplifiedhost)) + 
    geom_density(alpha = 0.7, position = 'stack') + 
    scale_x_date('Date', 
                 limits = ymd(c('2017-01-01', '2023-01-01')),
                 date_breaks = "1 year",
                 date_labels = "%Y",
                 expand = c(0,0)) + 
    scale_y_continuous('Probability Density', expand = c(0,0))  +
    
    
    scale_fill_manual(values = host_colours) + 
    scale_colour_manual(values = host_colours) + 
    global_theme
  
  
  
  return(plot)
}


PlotPhyloGeo <- function(treedata){
  require(rnaturalearthdata)
  require(tidytree)
  require(sf)
  require(tidyverse)
  
  
  # load base map
  map <- ne_countries(scale = "medium", returnclass = "sf")
  
  # Create tree-tibble from 
  tree_tbl <- as_tibble(treedata)
  
  most_recent_date <- tree_tbl %>%
    dplyr::select(label, height) %>%
    slice_min(height) %>%
    pull(label) %>%
    str_extract(., '\\d{4}-\\d{2}-\\d{2}') %>%
    ymd()
  
  nodes <- tree_tbl %>%
    dplyr::select(node,height, location1, location2, host_simplifiedhost) %>%
    mutate(year = decimal_date(most_recent_date) - as.numeric(height)) %>%
    st_as_sf(coords = c( 'location2', 'location1'), 
             crs = 4326)
  
  # Format Arrows
  arrows <-  tree_tbl %>%
    dplyr::select(node, parent, host_simplifiedhost) %>%
    left_join(tree_tbl %>%
                dplyr::select(node, location1, location2)) %>%
    left_join(tree_tbl %>%
                dplyr::select(node, location1, location2),
              by = join_by(parent== node)) %>%
    mutate(across(starts_with('location'), .fns = ~ as.numeric(.x))) %>%
    mutate(location1.y = case_when(location1.y == location1.x ~ location1.y + 0.00000001, 
                                   .default = location1.y),
           location2.y = case_when(location2.y == location2.x ~ location2.y + 0.00000001,
                                   .default = location2.y))
  
  # Plot in GGplot
  plot <- ggplot(map) +
    geom_sf() +
    
    # Plot Branches
    geom_curve(data = arrows, aes(x = location2.x, y = location1.x,
                                  xend = location2.y, yend = location1.y,
                                  colour= host_simplifiedhost),
               lwd = 0.3,
               curvature = 0.2) + 
    
    # Plot nodes
    geom_sf(data = nodes, shape = 21 , size = 1, aes(fill = host_simplifiedhost, colour = host_simplifiedhost))+
    
    # Set graphical scales
    scale_fill_manual('Host',
                      values = host_colours,
                      labels =  c('anseriformes-domestic' = 'Domestic Anseriformes',
                                  'anseriformes-wild' = 'Wild Anseriformes',
                                  'charadriiformes-wild' = 'Wild Charadriiformes',
                                  'galliformes-domestic' = 'Domestic Galliformes',
                                  'galliformes-wild' = 'Wild Galliformes',
                                  'environment' = 'Environment',
                                  'mammal' = 'Mammal',
                                  'other-bird' = 'Other',
                                  'human' = 'Human')) + 
    
    scale_colour_manual('Host',
                        values = host_colours,
                        labels =  c('anseriformes-domestic' = 'Domestic Anseriformes',
                                    'anseriformes-wild' = 'Wild Anseriformes',
                                    'charadriiformes-wild' = 'Wild Charadriiformes',
                                    'galliformes-domestic' = 'Domestic Galliformes',
                                    'galliformes-wild' = 'Wild Galliformes',
                                    'environment' = 'Environment',
                                    'mammal' = 'Mammal',
                                    'other-bird' = 'Other',
                                    'human' = 'Human')) +
    
    coord_sf(ylim = c(-60, 75), xlim = c(-180, 180), expand = TRUE) +
    scale_x_continuous(expand = c(0,0)) + 
    scale_y_continuous(expand = c(0,0)) +
    theme_void() + 
    
    
    theme(legend.position = 'none' ) 
  
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
                              full.names = TRUE)[-11] ,
                   './2024Sept16/south_america/ha_43112113_subsampled_traits_mcc.tree')


mcc_trees <- lapply(mcc_treefiles, read.beast)
names(mcc_trees) <- gsub('.*ha_|_subsampled_traits_mcc.tree', '', mcc_treefiles)

# Note the core requirement - this will take a long time to run in series.
# Change to futures -> multisession/multicore if you ever wish to try this again....
df_list <- mclapply(treefiles, function(x) x %>% 
           read.beast(.) %>% 
           lapply(., GetRootInfo) %>%
           bind_rows(),
           mc.cores = 13)

names(df_list) <- gsub('.*ha_|_subsampled.*', '' ,treefiles)

host_tmrca <- df_list %>%
  bind_rows(., .id = 'cluster_profile')

write_csv(host_tmrca, './2025Jan06/clusterprofile_host_tmrca.csv')
host_tmrca <- read_csv('./2025Jan06/clusterprofile_host_tmrca.csv')



dominant_reassortants <- read_csv('./2024Aug18/treedata_extractions/summary_reassortant_metadata_20240904.csv') %>%
  select(c(cluster_label,
            cluster_profile))  %>%
  distinct() %>% 
  mutate(cluster_label = gsub('_.*', '', cluster_label)) %>%
  filter(cluster_label %in% c('H5N8/2019/R7', 'H5N1/2020/R1','H5N1/2021/R1',
                              'H5N1/2021/R3','H5N1/2022/R7',  'H5N1/2022/R12' )) %>%
  pull(cluster_profile) %>%
  gsub('_', '', .) %>%
  as.double()
############################################## MAIN ################################################

# Filter the desired reassortants


# Plot the density of TMRCA, fill by origin  host
host_tmrca_list <- host_tmrca %>%

  filter(cluster_profile %in% dominant_reassortants) %>%
  base::split(.,f = .$cluster_profile) 

tmrca_plots <- lapply(host_tmrca_list, PlotTMRCA)


# Phylogeography (colour by host)
phylogeo <- lapply(mcc_trees[names(mcc_trees) %in% as.character(dominant_reassortants)], PlotPhyloGeo )


phylogeo[[6]]


# align plots vertically (so that each row corresponds to a reassortant)

plot_legend <- get_plot_component(phylogeo[[4]]+theme(legend.position = 'bottom'), 'guide-box-bottom', return_all = TRUE)

main <- cowplot::plot_grid(
  ggplot() + theme_void(),ggplot() + theme_void(),
  phylogeo[[1]],tmrca_plots[[1]], 
  phylogeo[[2]], tmrca_plots[[2]],
  phylogeo[[3]], tmrca_plots[[3]],
  phylogeo[[4]], tmrca_plots[[4]],
  phylogeo[[5]], tmrca_plots[[5]],
  phylogeo[[6]], tmrca_plots[[6]],
          nrow = 7,
          ncol = 2,
  rel_heights = c(0.1,1,1,1,1,1,1),
          rel_widths = c(1.25,0.75),
          align = 'vh',
          axis = 'rbt',
          labels = c('', '', 
                     "11111111",'',
                     "11211111" , '',
                     "11414114" , '',
                     "21111111", '', 
                     "32313212",'',
                     "43112113",''),
          vjust = -0.2)

cowplot::plot_grid(main,
                   plot_legend, 
                   rel_heights = c(1,0.1),
                   ncol = 1)

ggsave('~/Downloads/figure2.jpeg',
       height = 30,
       width = 25,
       units = 'cm',
       dpi = 360
       )

############################################## WRITE ###############################################


# For each reassortant, we need to most recent tipdate


# Infer the root age as a date

#



############################################## END #################################################
####################################################################################################
####################################################################################################