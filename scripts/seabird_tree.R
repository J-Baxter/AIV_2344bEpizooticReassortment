# seabird tree

# base tree


# requirements: is_seabird (tip colour)
# coloumn = location
# column = reassortant class


new_tree <- read.beast('./2025Feb26/globalsubsample/ha_global_SRD06_relaxLn_constant_mcc.tree')
meta <- read_csv('./2024-09-09_meta.csv') 

seabirds <- read_csv('~/Downloads/seabird_H5Nx_meta.csv') %>% 
  mutate(is_seabird = 1) %>%
  select(tipnames, is_seabird) 

metadata_in_tree <- meta  %>%
  filter(isolate_id %in% str_extract(new_tree@phylo$tip.label, "EPI_ISL_(china_){0,1}\\d+[^.|]*")) %>%
  dplyr::select(tipnames,
                collection_regionname, 
                isolate_id,
                host_simplifiedhost,
                cluster_profile,
                host_class ) %>%
  left_join(summary_data %>% select(cluster_profile, group2)) %>%
  mutate(collection_regionname = case_when(grepl('europe', collection_regionname) ~ 'europe',
                                           grepl('africa', collection_regionname) ~ 'africa',
                                           grepl('asia', collection_regionname) ~ 'asia',
                                           grepl('(central|northern) america|caribbean', collection_regionname) ~ 'central & northern america',
                                           grepl('south america|southern ocean', collection_regionname) ~ 'south america',
                                           grepl('australia|melanesia', collection_regionname) ~ 'australasia',
                                           .default = collection_regionname)) %>%
  left_join(seabirds) %>%
  drop_na(collection_regionname)

drop <- new_tree@phylo$tip.label[str_extract(new_tree@phylo$tip.label, "EPI_ISL_(china_){0,1}\\d+[^.|]*") %in% (meta  %>%
                                   filter(isolate_id %in% str_extract(new_tree@phylo$tip.label, "EPI_ISL_(china_){0,1}\\d+[^.|]*")) %>%
                                   filter(is.na(collection_regionname)) %>%
                                   pull(isolate_id))]

new_names <- cbind.data.frame(old = new_tree@phylo$tip.label,
                              new = str_extract(new_tree@phylo$tip.label, "EPI_ISL_(china_){0,1}\\d+[^.|]*"))
 new_tree %>% 
  drop.tip(drop) %>%
  rename_taxa( new_names, old, new) %>%
  left_join(metadata_in_tree %>% 
              rename(label = isolate_id),
            by = 'label') %>%

  ggtree(mrsd =  "2024-03-18") + 
  theme_tree2(plot.margin = unit(c(1,1,1,1), units = "cm"),
              axis.text.x = element_text(size = 8),
              axis.title.x = element_text(size = 10),
              legend.title = element_text(size = 10),
              legend.text = element_text(size = 8)) +
  
   scale_x_continuous(
     #limits = c(2000, 2023),
     'Time',
     breaks = seq(2016, 2024, 1)) +
  
  # tip colour + shape = new sequences
  geom_tippoint(aes(colour = !is.na(is_seabird),
                    shape = !is.na(is_seabird)),
                size = 3, 
                alpha = 0.9) +
  
  geom_tiplab(aes(colour = !is.na(is_seabird)),
              #align = TRUE, 
              size = 0) +
  
  scale_shape_manual(values = c("TRUE" = 18),
                     'Seabird',
                     guide = 'none') +
  scale_colour_manual(values = c("TRUE" = 'blue'), 
                      'Seabird',
                      guide = 'none') + 
  
  
  # node colour to show pp support
 # new_scale_colour()+
  #geom_nodepoint(aes(colour = posterior), alpha = 0.7) +
  #scale_color_distiller(palette = 'YlOrRd', direction = 1, 'Posterior Support',
                        #guide = guide_colourbar(order = 4)) + 
  
  geom_fruit(geom = geom_tile,
             mapping = aes(fill = collection_regionname),
             #width = 4,
             #colour = "white",
             #pwidth = 1.2,
             offset = 0.03) + 
  scale_fill_manual('Continent', values = region_colours, labels = str_to_title,
                    guide = guide_legend(keywidth = 1.5, keyheight = 1, ncol = 1, order = 1)) + 
  
  new_scale_fill()+
  geom_fruit(geom = geom_tile,
             mapping = aes(fill = host_simplifiedhost),
             #width = 4,
             #colour = "white",
             #pwidth = 1.2,
             offset = 0.03) +
  scale_fill_manual('Host', values = host_colours, labels = str_to_title,
                    guide = guide_legend(keywidth = 1.5, keyheight = 1, ncol = 1, order = 2)) +
   
   new_scale_fill()+
   geom_fruit(geom = geom_tile,
              mapping = aes(fill = group2),
              #width = 4,
              #colour = "white",
              #pwidth = 1.2,
              offset = 0.03) +
   scale_fill_brewer('Cluster Class', palette = 'Set1',labels = str_to_title,
                     guide = guide_legend(keywidth = 1.5, keyheight = 1, ncol = 1, order = 2)) +
  theme(legend.position = c(0.2,0.6),
        # legend.position = "bottom",       # Place legends at the bottom
        legend.box = "vertical",
        legend.direction = 'vertical'
  )
  
  