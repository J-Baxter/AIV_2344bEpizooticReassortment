require(rnaturalearthdata)
require(tidytree)
require(sf)
require(tidyverse)
require(seraphim)

# Import tree data
mcc_tree <- read.beast(mcc_treefiles[2])
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
allTrees <- scan(file = posterior_treefiles[2],
                 what = '',
                 sep = '\n',
                 quiet = T)


localTreesDirectory = "./temp_tree_dir/"
do.call(file.remove, list(list.files("./temp_tree_dir/", full.names = TRUE)))

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
  #dplyr::select(node,height_median, location1, location2, host_simplifiedhost, label) %>%
  mutate(height_median= as.numeric(height_median)) %>%
  replace_na(list(height_median = 0)) %>%
  mutate(year = most_recent_date - height_median) %>%
  filter(! node %in% c(1164, 1198, 1165, 1169, 1166, 1163, 1162, 1121, 1123, 1124)) %>%
  #group_by(parent) %>% filter(all(location2 > 0)|all(location2<0)) %>% ungroup() %>%
  # Convert to POINT & set coordinate system
  st_as_sf(coords = c( 'location2', 'location1'), 
           crs = 4326)

# Format Arrows
edges <-  tree_tbl %>%
  dplyr::select(node, parent) %>%
  filter(node %in% nodes_sf$node & parent %in% nodes_sf$node) %>%
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


replacement <- tree_tbl %>% 
  filter(!node %in% edges$node) %>% 
  drop_na(label) %>%
  arrange(parent) %>%
  mutate(vis_pair = c(1,2,2,3,3,4,4,5,5)) %>%
  group_by(vis_pair) %>%
  arrange(desc(height),.by_group = TRUE) %>%
  dplyr::select(node, label, location1, location2, vis_pair) %>%
  mutate(position = case_when(row_number()==1 ~ 'start', 
                              .default = 'end')) %>%
  filter(n()>1) %>%
  pivot_wider(names_from = 'position', values_from = c(node, label, location1, location2)) %>%
  mutate(across(starts_with('location'), .fns = ~ as.numeric(.x))) %>%
  rename(start_lat = location1_start,
         start_lon = location2_start,
         end_lat =  location1_end,
         end_lon = location2_end,
         node = node_start,
         parent = node_end) %>%
  ungroup() %>%
  dplyr::select(-c(vis_pair, starts_with('label')))
  
edges %<>% bind_rows(replacement)
  
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





map <- ne_countries(scale = "medium", returnclass = "sf")

 alt_plot <- ggplot() +
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
   
   #geom_sf_text(data = nodes_sf,
               # size = 1.5, 
                #aes( label = node)) +
  
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