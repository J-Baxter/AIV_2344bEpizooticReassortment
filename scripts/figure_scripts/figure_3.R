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


FormatPhyloGeo <- function(mcc_file, posterior_file){
  require(rnaturalearthdata)
  require(tidytree)
  require(sf)
  require(tidyverse)
  require(seraphim)
  
  
  mcc_tree <- read.beast(mcc_file)
  tree_tbl <- as_tibble(mcc_tree)
  
  # Guess most_recent_date
  most_recent_date <- tree_tbl %>%
    dplyr::select(label, height) %>%
    slice_min(height) %>%
    pull(label) %>%
    str_extract(., '\\d{4}-\\d{2}-\\d{2}') %>%
    ymd()
  
  # Scan posterior tree file
  allTrees <- scan(file = posterior_file,
                   what = '',
                   sep = '\n',
                   quiet = T)
  
  
  localTreesDirectory = "./temp_tree_dir/"
  do.call(file.remove, list(list.files("./temp_tree_dir/", full.names = TRUE)))
  
  burnIn <- 0
  randomSampling <- TRUE
  nberOfTreesToSample <- 100
  mostRecentSamplingDatum <- decimal_date(most_recent_date%>% ymd())
  coordinateAttributeName <- "location"
  treeExtractions(localTreesDirectory,
                  allTrees,
                  burnIn, 
                  randomSampling, 
                  nberOfTreesToSample, 
                  mostRecentSamplingDatum,
                  coordinateAttributeName,
                  nberOfCores = 6)
  
  mcc_tab <- mccExtractions(readAnnotatedNexus(TREEFILE), mostRecentSamplingDatum)
  
  # Step 4: Estimating the HPD region for each time slice ----
  nberOfExtractionFiles <- nberOfTreesToSample
  prob <- 0.95
  precision <- 0.08 # time interval that will be used to define the successive time slices
  startDatum <- min(mcc_tab[,"startYear"])
  
  
  # Format Polygons
  polygons <- suppressWarnings(spreadGraphic2(localTreesDirectory, 
                                              nberOfExtractionFiles,
                                              prob, 
                                              startDatum, 
                                              precision))
  
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
    dplyr::select(node,height, location1, location2, host_simplifiedhost, label) %>%
    mutate(year = decimal_date(most_recent_date) - as.numeric(height)) %>%
    
    # Convert to POINT & set coordinate system
    st_as_sf(coords = c( 'location2', 'location1'), 
             crs = 4326)
  
  # Format Arrows
  arrows <-  tree_tbl %>%
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
  
  arrows_sf <- arrows %>%
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
              arrows = arrows_sf,
              nodes = nodes_sf)
  
  
  return(out)
}



############################################## DATA ################################################

# read in cross-continental tree files
treefiles <- c(list.files('./2024Aug18/reassortant_subsampled_outputs/traits_1000trees',
                        pattern = 'ha',
                        full.names = TRUE))

mcc_treefiles <- c(list.files('./2024Aug18/reassortant_subsampled_outputs/traits_mcc',
                              pattern = 'ha',
                              full.names = TRUE))


mcc_trees <- lapply(mcc_treefiles, read.beast)
names(mcc_trees) <- gsub('.*ha_|_subsampled_traits_mcc.tree', '', mcc_treefiles)
#n <- names(mcc_trees) 

#names(mcc_trees) <- n
# Note the core requirement - this will take a long time to run in series.
# Change to futures -> multisession/multicore if you ever wish to try this again....
#df_list <- mclapply(treefiles, function(x) x %>% 
 #          read.beast(.) %>% 
#           lapply(., GetRootInfo) %>%
#           bind_rows(),
#           mc.cores = 13)
#
#names(df_list) <- gsub('.*ha_|_subsampled.*', '' ,treefiles)

#host_tmrca <- df_list %>%
 # bind_rows(., .id = 'cluster_profile')

#write_csv(host_tmrca, './2025Jan06/clusterprofile_host_tmrca.csv')
#host_tmrca <- read_csv('./2025Jan06/clusterprofile_host_tmrca.csv')



dominant_reassortants <- read_csv('./2024Aug18/treedata_extractions/summary_reassortant_metadata_20240904.csv') %>%
  dplyr::select(c(cluster_label,
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
#host_tmrca_list <- host_tmrca %>%

 # filter(cluster_profile %in% dominant_reassortants) %>%
  #base::split(.,f = .$cluster_profile) 

#tmrca_plots <- lapply(host_tmrca_list, PlotTMRCA)


# Phylogeography (colour by host)
phylogeo <- mapply(PlotPhyloGeo, 
                   mcc_treefiles[sapply(dominant_reassortants, function(x) which(grepl(x, mcc_treefiles)))][-1],
                   treefiles[sapply(dominant_reassortants, function(x) which(grepl(x, treefiles)))][-1],
                   )
  
  
  lapply(, PlotPhyloGeo )


phylogeo[[6]]


north_america_ha <- read.beast('./2024Aug18/region_subsampled_outputs/traits_mcc/ha_northamerica_subsampled_traits_mcc.tree')
mooflu_data <- as_tibble(north_america_ha) %>% 
  select(label, height, host_simplifiedhost, location1, location2) %>%  mutate(across(c(location1, location2, height), .fns = ~as.numeric(.x)))

mooflu <- tree_subset(north_america_ha, 622)  

 south_america <- read.beast( '~/Downloads/ha_43112113_subsampled_traits_mcc(1).tree')

PlotPhyloGeo(south_america)  +  coord_sf(ylim = c(-60, 75), xlim = c(-150, -30), expand = TRUE)
phylogeo[[10]] 
# align plots vertically (so that each row corresponds to a reassortant)

plot_legend <- get_plot_component(phylogeo[[4]]+theme(legend.position = 'bottom'), 'guide-box-bottom', return_all = TRUE)

main <- cowplot::plot_grid(
  ggplot() + theme_void(),ggplot() + theme_void(),
  phylogeo[[5]],tmrca_plots[[5]], 
  phylogeo[[4]], tmrca_plots[[4]],
  phylogeo[[2]], tmrca_plots[[2]],
  phylogeo[[1]], tmrca_plots[[1]],
  phylogeo[[6]], tmrca_plots[[6]],
  phylogeo[[3]], tmrca_plots[[3]],
          nrow = 7,
          ncol = 2,
  rel_heights = c(0.1,1,1,1,1,1,1),
          rel_widths = c(1.25,0.75),
          align = 'vh',
          axis = 'rbt',
          labels = c('', '', 
                     "H5N8/2019/R7",'',
                     "H5N1/2020/R1", '', 
                     "H5N1/2021/R1" , '',
                     "H5N1/2021/R3",'',
                     "H5N1/2022/R7",'',
                     "H5N1/2022/R12" , '' ),
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