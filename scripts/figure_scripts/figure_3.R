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
#library(beastio)
library(treeio)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(cowplot)

# User functions
FormatPhyloGeo <- function(mcc_file, posterior_file){
  require(rnaturalearthdata)
  require(tidytree)
  require(sf)
  require(tidyverse)
  require(seraphim)
  
  # Import tree data
  mcc_tree <- read.beast(mcc_file)
  tree_tbl <- as_tibble(mcc_tree)
  
  # Guess most_recent_date
  most_recent_date <- tree_tbl %>%
    slice_min(height_median, n = 1, with_ties = FALSE) %>%
    pull(label) %>%
    str_extract(., '\\d{4}-\\d{2}-\\d{2}') %>%
    ymd() %>%
    decimal_date() %>%
    unique()
  
  
  # Extract TMRCA
  start_date <- tree_tbl %>%
    slice_max(height_median, n = 1, with_ties = FALSE) %>% 
    pull(height_median) %>% 
    as.numeric() %>% 
    subtract(most_recent_date,.)
  print(start_date)

  # Scan posterior tree file
  allTrees <- scan(file = posterior_file,
                   what = '',
                   sep = '\n',
                   quiet = T)
  
  
  localTreesDirectory = "./2025Jun10/temp_tree_dir/"
  do.call(file.remove, list(list.files("./2025Jun10/temp_tree_dir/", full.names = TRUE)))
  
  burnIn <- 0
  randomSampling <- TRUE
  nberOfTreesToSample <- 100
  #mostRecentSamplingDatum <- most_recent_date
  coordinateAttributeName <- "location"
  treeExtractions(localTreesDirectory,
                  allTrees,
                  burnIn, 
                  randomSampling, 
                  nberOfTreesToSample, 
                  most_recent_date,
                  coordinateAttributeName,
                  nberOfCores = 8)
  
  # Step 4: Estimating the HPD region for each time slice ----
  
  # Format Polygons
  polygons <- suppressWarnings(spreadGraphic2(localTreesDirectory, 
                                              nberOfExtractionFiles = nberOfTreesToSample,
                                              prob = 0.95, 
                                              startDatum = start_date, 
                                              precision = 0.08))
  
  # Conver polgons to GEOMETRY, set coordinate system and add values
  polygons_sf <- lapply(polygons, st_as_sf)  %>%
    bind_rows() %>%
    pivot_longer(cols = starts_with('2'), names_to = 'year', values_to = 'value') %>%
    filter(value == 1) %>%
    mutate(year = as.numeric(year)) %>%
    dplyr::select(-value) %>%
    st_set_crs(4326) 
  
  
  # Format Nodes
  nodes_sf <- tree_tbl %>%
    dplyr::select(node,height_median, location1, location2, host_simplifiedhost, label) %>%
    mutate(height_median= as.numeric(height_median)) %>%
    replace_na(list(height_median = 0)) %>%
    mutate(year = most_recent_date - height_median) %>%
    
    # Convert to POINT & set coordinate system
    st_as_sf(coords = c( 'location2', 'location1'), 
             crs = 4326)
  
  # Format Arrows
  edges <-  tree_tbl %>%
    dplyr::select(node, parent) %>%
    left_join(tree_tbl %>%
                dplyr::select(node, location1, location2)) %>%
    left_join(tree_tbl %>%
                dplyr::select(node, location1, location2),
              by = join_by(parent== node)) %>%
    mutate(across(starts_with('location'), .fns = ~ as.numeric(.x))) %>%
    mutate(location1.y = case_when(location1.y == location1.x ~ location1.y + 0.00000001, 
                                   .default = location1.y),
           location2.y = case_when(location2.y == location2.x ~ location2.y + 0.00000001,
                                   .default = location2.y)) %>%
    rename(start_lat = location1.x,
           start_lon = location2.x,
           end_lat = location1.y,
           end_lon = location2.y)
  
  edges_sf <- edges %>%
    rowid_to_column(var = 'id') %>%
    pivot_longer(starts_with(c('start', 'end')),
                 names_to = c("type", ".value"),
                 names_sep = "_") %>%
    
    # Convert coordinate data to sf POINT
    st_as_sf( coords = c("lon", "lat"),
              crs = 4326) %>%
    
    # Convert POINT geometry to MULTIPOINT, then LINESTRING
    group_by(id) %>% 
    summarise(do_union = FALSE) %>% 
    st_cast("LINESTRING") %>% 
    
    # Convert rhumb lines to great circles
    st_segmentize(units::set_units(20, km)) %>%
    
    # Wrap dateline correctly
    st_wrap_dateline()
  
  
  out <- list(polygons = polygons_sf,
              edges = edges_sf,
              nodes = nodes_sf)
  
  
  return(out)
}


PlotPhyloGeo <- function(phylogeo_list){
  polygons_sf <- phylogeo_list[['polygons']]
  edges_sf <- phylogeo_list[['edges']]
  nodes_sf <- phylogeo_list[['nodes']]
  
  # load base map
  map <- ne_countries(scale = "medium", returnclass = "sf")
  
  # Plot in GGplot
  plot <- ggplot() +
    geom_sf(data = map) +
    
    # Plot HPD polygons
    geom_sf(data = polygons_sf, 
            aes(fill = year), 
            lwd = 0, 
            alpha = 0.04) + 
    
    # Plot Branches
    geom_sf(data = edges_sf,
            lwd = 0.2) + 
    
    # Plot nodes
    geom_sf(data = nodes_sf,
            size = 1.5, 
            aes(fill = year, colour = year, shape = is.na(label)))+
    
    scale_shape_manual(values = c(19,1),
    ) + 
    
    # Set graphical scales - must be fixed
    scale_fill_viridis_c('Year',
                         limits = c(2018, 2024.4),
                         #breaks=c(2019, 2020,2021,2022,2023,2024),
                         direction = -1,
                         option = 'C')+
    
    scale_colour_viridis_c('Year',
                           limits = c(2018, 2024.4),
                           #breaks=c(2019, 2020,2021,2022,2023,2024),
                           direction = -1,
                           option = 'C')+
    
    coord_sf(ylim = c(-60, 80),
             xlim = c(-185, 185),
             expand = TRUE) +
    
    scale_x_continuous(expand = c(0,0)) + 
    scale_y_continuous(expand = c(0,0)) +
    
    guides(colour = guide_colourbar(
      theme = theme(
        legend.key.height  = unit(0.75, "lines"),
        legend.key.width = unit(10, "lines")),
      #title.position = 'left',
      title.vjust = 1,
      position = 'bottom'), 
      
      fill = guide_colourbar(
        theme = theme(
          legend.key.height  = unit(0.75, "lines"),
          legend.key.width = unit(10, "lines")),
        title.vjust =1,
        #title.position = 'left',
        position = 'bottom'), 
      
      shape = 'none') +
    
    theme_void() + 
    theme(plot.margin=grid::unit(c(0,0,0,0), "mm"),
          legend.text = element_text(size = 8),
          legend.position = 'none', 
          panel.spacing = unit(2, "lines"), 
          strip.background = element_blank()) 
  
  return(plot)
}

############################################## DATA ################################################

# read in cross-continental tree files
treefiles <- list.files('./2024Aug18/reassortant_subsampled_outputs/traits_1000trees',
                        pattern = 'ha',
                        full.names = TRUE)

mcc_treefiles <-list.files('./2024Aug18/reassortant_subsampled_outputs/traits_mcc',
                              pattern = 'ha',
                              full.names = TRUE)


dominant_reassortants <- read_csv('./2024Aug18/treedata_extractions/summary_reassortant_metadata_20240904.csv') %>%
  dplyr::select(c(cluster_label,
            cluster_profile))  %>%
  distinct() %>% 
  mutate(cluster_label = gsub('_.*', '', cluster_label)) %>%
  filter(cluster_label %in% c('H5N8/2019/R7', 'H5N1/2020/R1','H5N1/2021/R1',
                             'H5N1/2022/R7',  'H5N1/2022/R12' )) %>%
  pull(cluster_profile) %>%
  gsub('_', '', .) %>%
  as.double()

mcc_treefiles <-  c(mcc_treefiles[sapply(dominant_reassortants, function(x) which(grepl(x, mcc_treefiles)))][-3],
                    #'~/Downloads/ha_43112113_subsampled_traits_mcc(1).tree',
                    '~/Downloads/32313212/ha_32313212_subsampled_traits_mcc.tree')

posterior_treefiles <- c( treefiles[sapply(dominant_reassortants, function(x) which(grepl(x, treefiles)))][-3],
                         # '~/Downloads/ha_43112113_subsampled_traits_1000.trees',
                          '~/Downloads/32313212/ha_32313212_subsampled_traits_1000.trees')


############################################## MAIN ################################################

# Phylogeography (colour by node date)
formatted_phylogeos <- mapply(FormatPhyloGeo, 
                              mcc_treefiles,
                              posterior_treefiles,
                              SIMPLIFY = FALSE)

#test <- FormatPhyloGeo(mcc_treefiles[3], posterior_treefiles[3])
phylogeo <- lapply(formatted_phylogeos, PlotPhyloGeo)
#phylogeo[[1]]
#north_america_ha <- read.beast('./2024Aug18/region_subsampled_outputs/traits_mcc/ha_northamerica_subsampled_traits_mcc.tree')
#mooflu_data <- as_tibble(north_america_ha) %>% 
 # select(label, height, host_simplifiedhost, location1, location2) %>%  mutate(across(c(location1, location2, height), .fns = ~as.numeric(.x)))

#mooflu <- tree_subset(north_america_ha, 622)  

#south_america <- read.beast( '~/Downloads/ha_43112113_subsampled_traits_mcc(1).tree')

#PlotPhyloGeo(south_america)  +  coord_sf(ylim = c(-60, 75), xlim = c(-150, -30), expand = TRUE)

# align plots vertically (so that each row corresponds to a reassortant)

plot_legend <- get_plot_component(phylogeo[[4]] + theme(legend.position = 'bottom') , 'guide-box-bottom', return_all = TRUE)

#main <- cowplot::plot_grid(
  #ncol = 2,
  #phylogeo[[3]],
  #phylogeo[[2]], 
  #phylogeo[[6]], 
  #phylogeo[[4]], 
  #phylogeo[[1]], 
  #phylogeo[[5]], 
  #labels = c("H5N8/2019/R7",
             #"H5N1/2020/R1", 
             #"H5N1/2021/R1" , 
             #"H5N1/2021/R3",
             #"H5N1/2022/R7",
            # "H5N1/2022/R12"  ),
  #label_size #= 8)

cowplot::plot_grid(
  ncol = 2,
  phylogeo[[5]],
  phylogeo[[2]], 
  phylogeo[[4]],
  phylogeo[[3]], 
  phylogeo[[1]],
  plot_legend,
  labels = c( "H5N8/2019/R7",
              "H5N1/2020/R1",
             "H5N1/2021/R1",
             "H5N1/2022/R12",
             "H5N1/2022/R7" ),
  label_size = 8)

# Chronological order
cowplot::plot_grid(
  ncol = 2,
  phylogeo[[5]],
  alt_plot,
  phylogeo[[4]], 
  phylogeo[[3]],
  phylogeo[[1]], 
  plot_legend,

  labels = c( "H5N8/2019/R7",
              "H5N1/2020/R1 (AIV07, C)",
              "H5N1/2021/R1 (AIV09, AB)",
             "H5N1/2022/R7 (B3.2)" ,
             "H5N1/2022/R12 (AIV48, BB)"
            ),
  label_x = 0, hjust = 0,
  label_size = 8) 
  




cowplot::plot_grid(main,
                   plot_legend, 
                   rel_heights = c(1,0.1),
                   ncol = 1)

ggsave('~/Downloads/flu_plots/phylogeo.jpeg',
       height = 15,
       width = 20,
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