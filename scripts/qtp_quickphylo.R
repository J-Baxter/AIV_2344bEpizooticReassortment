library(tidyverse)
library(magrittr)
library(ggmap)
library(elevatr)
library(ape)
library(lubridate)
library(auk)
library(sf)
library(seraphim)
library(giscoR)
library(rnaturalearth)
library(rnaturalearthhires)
library(diagram)
library(treeio)
library(TreeTools)
library(ggtree)
library(ggtreeExtra)
library(ggnewscale)

# User functions
source('./scripts/FormatBirds.R')
birds <- read_csv('bird_taxonomy.csv')



GroupSequences <- function(aln, snp_threshold = 0){
  require(igraph)
  
  # Ensure alignment is correctly formatted
  if(class(aln) != 'PhyDat'){
    aln_formatted <- as.phyDat(aln) 
  }else{
    aln_formatted <- aln
  }
  
  # Calculate hamming distance
  hd_normalised <- dist.hamming(aln_formatted) %>%
    as.matrix()
  hd_raw <- hd_normalised * ncol(aln)
  
  # Obtain groups of sequences for which HD < SNP threshold
  
  if( any(hd_raw[lower.tri(hd_raw, diag = FALSE)] <= snp_threshold)){
    groups <- which(hd_raw <= snp_threshold, 
                    arr.ind = TRUE) %>%
      dplyr::as_tibble(rownames = 'tipnames') %>% 
      filter(row !=col) %>%
      dplyr::select(-tipnames) %>%
      
      # Infer network from HDs
      igraph::graph_from_data_frame(.,
                                    directed = F) %>%
      
      components() %>%
      getElement('membership') %>%
      stack() %>%
      as_tibble() %>%
      mutate(ind = as.numeric(as.character(ind))) %>%
      mutate(tipnames = map_chr(ind, ~ rownames(aln)[.x])) %>%
      dplyr::select(c(tipnames, values)) %>%
      dplyr::distinct() %>%
      dplyr::rename(sequence_group = values)
  }else{
    print('all sequences above threshold.')
    groups <- tibble(tipnames = rownames(aln)) %>%
      rowid_to_column(., var = 'sequence_group')
  }
  
  
  out <- tibble(tipnames = rownames(aln)) %>%
    left_join(groups) %>%
    mutate(sequence_group = 
             ifelse(is.na(sequence_group), 
                    max(sequence_group, na.rm = T) + row_number() + 1, 
                    sequence_group))
  
  
  return(out)
}


# Down sample formatted QTP data
qinghai_ha <- read.dna('./new_data/2024Dec23_CAG_h5nX_ha_formatted_aligned.fasta',
                       format = 'fasta',
                       as.matrix = T)


ha_groups <- GroupSequences(qinghai_ha, snp_threshold = 3) 

qtp_data <- read_csv('./new_data/cag_data_2025Jan06.csv') %>%
  left_join(ha_groups)


# subsample by location, month-year and sequence group
subsample <- qtp_data %>%
  filter(tipnames %in% rownames(qinghai_ha) ) %>%
  mutate(collection_yearmonth = format.Date(collection_date, '%Y-%m')) %>%
  group_by( collection_subdiv2, collection_yearmonth, sequence_group) %>%
  slice_sample(n = 1)

# subsample alignment
qinghai_ha_subsample <- qinghai_ha[rownames(qinghai_ha) %in% subsample$tipnames,]


write.FASTA(qinghai_ha_subsample, 
            './new_data/2025Feb12_qinghai_ha_subsample.fasta')



### post blast
# Locally run BLASTN, using all H5 2344b HPAI (note this works for HA, will have to update in future) 
# on GISAID between 2019-01-01 and 2024-12-31

# Import BLASTN CSV
blast_results <- read_csv('~/Downloads/results.csv') %>%
  mutate(Isolate_Id = gsub('EPI\\d+\\|', '' ,subject) %>%
           gsub('\\|.*', '', .)) 

# context alignments and metadata
gisaid_meta <- read_csv('~/Downloads/gisaid_epiflu_isolates.csv')
context_seqs <- read.dna('~/Downloads/gisaid_epiflu_sequence.fasta', format = 'fasta')

gisaid_meta %<>%
  select(Isolate_Id,
         Isolate_Name,
         Subtype,
         Genotype,
         Lineage,
         Clade,
         Location,
         Host,
         Collection_Date) %>%
  separate_wider_delim(Location, delim = ' / ', names = paste0('location', '_' ,1:5) , too_few = 'align_start')
# left join results with metadata (by subject)

context_sequences <- blast_results %>%
  filter(bitscore > 50) %>%
  left_join(gisaid_meta) %>%
  group_by(query) %>% # maybe add location as a grouping
  slice_sample(n = 20) %>%
  ungroup() %>%
  distinct(Isolate_Id)

context_subsample <- sapply(context_sequences$Isolate_Id, function(x) context_seqs[grepl(x, names(context_seqs))], simplify = T) %>%
  do.call(c, .)

context_meta <- gisaid_meta %>%
  filter(Isolate_Id %in% context_sequences$Isolate_Id) %>%
  
  
  # Coordinates
  unite(location_query, sep = ' ', location_2, location_3, location_4, remove = FALSE, na.rm = TRUE) %>%
  
  # geocode to obtain lat-lon coordinates
  mutate(location = geocode(location_query)) %>%  
  unnest(location) %>%
  st_as_sf(coords = c('lon', 'lat'), crs = 4326, remove = FALSE) %>%
  
  
  
  # obtain elevations for inferred coordinates
  rowwise() %>%
  get_elev_point(prj = 4326, src = "aws")

context_meta_2 <- context_meta %>%
  rowwise() %>%
  mutate(primary_com_name = FormatBirds(tolower(Host))) %>%
  left_join(birds)

# Make tipnames
context_meta_2 %<>%
  mutate(tipnames = paste(gsub('*. \\/ ', '', Subtype), 
                          gsub('\\.', '', Clade),
                          Isolate_Id,
                          tolower(order),
                          tolower(location_2),
                          NA,
                          Collection_Date,
                          sep = '|')) 
context_meta_3 <- context_meta_2 %>%
  mutate(tipnames = gsub(' |,|_', '_', tipnames))

old_names <- names(context_subsample)
new_seqnames <- c()
for (i in 1:length(old_names)){
  
  new_name <-  context_meta_3 %>% 
    filter(Isolate_Id %in% str_extract(old_names[[i]] , 'EPI_ISL_\\d+')) %>%
    pull(tipnames)
  
  new_seqnames[i] <- new_name

  
}

names(context_subsample) <- new_seqnames


# save outputs (unaligned sequences and metadata)
all_sequences <- c(as.list(qinghai_ha_subsample), context_subsample)

write.FASTA(all_sequences[!duplicated(names(all_sequences))], '~/Downloads/2025Feb12_qinghai_ha_subsample_withBLASTcontext.fasta')


#write continuous phylo traits file
aln <- all_sequences[!duplicated(names(all_sequences))]
traits <- context_meta_3 %>% filter(tipnames %in% names(aln)) %>%
  select(tipnames, lon, lat) %>%
  st_drop_geometry() %>%
  bind_rows(qtp_data %>% filter(tipnames %in% names(aln)) %>%
              select(tipnames, geometry) %>%
              mutate(geometry = gsub('c\\(|\\)', '', geometry)) %>%
              separate_wider_delim(geometry, delim = ', ', names = c('lon', 'lat')) %>%
              mutate(across(starts_with('l'), .fns = ~ as.numeric(.x))))


temp <- read.dna('~/Downloads/2025Feb12_qinghai_ha_subsample_withBLASTcontext_aligned_2.fasta', format = 'fasta')


write_delim(traits %>% distinct(tipnames, .keep_all = TRUE) %>% filter(tipnames %in% rownames(temp))
, '~/Downloads/2025Feb12_qinghai_ha_subsample_withBLASTcontext_aligned_2.txt',
            delim = '\t',
            quote = 'needed')


##### 
# Read MCC tree

library(treeio)
library(TreeTools)


beast_tree <- read.beast('~/Downloads/2025Feb12_qinghai_ha_subsample_withBLASTcontext_mcc.tree')

library(ggtree)

# MCC tree with hosts, countries and QTP sequences highlighted
mrsd <- '2024-12-13'

data_in_tree <-  context_meta_3 %>% filter(tipnames %in% names(aln)) %>%
  select(tipnames, location_2, order) %>%
  st_drop_geometry() %>%
  rename(host_order = order,
         collection_country = location_2) %>%
  bind_rows(qtp_data %>% filter(tipnames %in% names(aln))) %>%
  mutate(is_new = case_when(!is.na(collection_subdiv1_hanzi) ~ 'yes', .default = 'no')) %>%
  select(tipnames, collection_country, host_order, is_new, elevation, elevation_class) %>%
  rename(label = tipnames) %>%
  distinct(label, , .keep_all = TRUE)


temp_tree <- beast_tree %>%
  left_join(data_in_tree)  %>%
  
ggtree(mrsd = mrsd) + 
  
  
  # tip colour + shape = new sequences
  geom_tippoint(aes(colour = is_new,
                    shape = is_new),
                size = 3, 
                alpha = 0.9) +
  
  geom_tiplab(aes(colour = is_new),
              #align = TRUE, 
              size = 0) +
  
  scale_shape_manual(values = c("yes" = 18),
                     'QTP Sequences',
                     guide = 'none') +
  scale_colour_manual(values = c("yes" = 'blue'), 
                      'QTP Sequences',
                      guide = 'none') + 
  new_scale_fill()+
  geom_fruit(geom = geom_tile,
             mapping = aes(fill = collection_country),
             width = 0.5,
             colour = "white",
             pwidth = 1.2,
             offset = 0.1) +
  scale_fill_discrete('Country',
                      #alpha = 0.99, 
                      guide = guide_legend(keywidth = 1.5, keyheight = 1, ncol = 2, order = 1)) +
  #scale_fill_d3(name = 'Country', 
                #palette ='category20', 
                #alpha = 0.99, 
                #guide = guide_legend(keywidth = 1.5, keyheight = 1, ncol = 2, order = 1)) +
  
  new_scale_fill()+
  geom_fruit(geom = geom_tile,
             mapping = aes(fill = host_order),
             width = 0.5,
             colour = "white",
             pwidth = 1.2,
             offset = 00.1) +
  scale_fill_brewer('Host Order',
                    guide = guide_legend(keywidth = 1.5, keyheight = 1, ncol = 1, order = 3)) + 
  
  theme_tree2()


# phylogeography
# load base map
map <- ne_countries(continent = "asia", scale = "medium", returnclass = "sf")

# Create tree-tibble from 
tree_tbl <- as_tibble(beast_tree)

most_recent_date <- tree_tbl %>%
  dplyr::select(label, height) %>%
  slice_min(height) %>%
  pull(label) %>%
  str_extract(., '\\d{4}-\\d{2}-\\d{2}') %>%
  ymd()

nodes <- tree_tbl %>%
  dplyr::select(node,height, location1, location2) %>%
  mutate(year = decimal_date(most_recent_date) - as.numeric(height)) %>%
  st_as_sf(coords = c( 'location1', 'location2'), 
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
                                 .default = location2.y))

# Plot in GGplot
plot <- ggplot(map) +
  geom_sf() +
  
  # Plot Branches
  geom_curve(data = arrows, aes(y = location2.x, x = location1.x,
                                yend = location2.y, xend = location1.y),
             lwd = 0.3,
             curvature = 0.2) + 
  
  # Plot nodes
  geom_sf(data = nodes, shape = 21 , size = 1, aes(fill = year, colour = year))+
  
  
  coord_sf(ylim = c(10, 50), xlim = c(40, 160), expand = TRUE) +
  scale_x_continuous(expand = c(0,0)) + 
  scale_y_continuous(expand = c(0,0)) +
  theme_void() + 
  
  
  theme(legend.position = 'none' ) 
