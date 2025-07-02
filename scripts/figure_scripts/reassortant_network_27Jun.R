my_graph <- reassortant_ancestral_changes %>%
  dplyr::select(ends_with('label'), 
                segments_changed, 
                time_since_parent)%>% 
  drop_na() %>%
  relocate(parent_label,
           cluster_label) %>%
  
  # Make graph
  as_tbl_graph(.) %>%
  
  # Add node data
  activate(nodes) %>%
  left_join(reassortant_ancestral_changes %>% 
              dplyr::select(name = cluster_label, cluster_class)) %>% mutate(importance = centrality_degree()) %>%
  left_join(updated %>% dplyr::select(name = cluster_label, cluster_region, cluster_tmrca)) %>%
  
  mutate(is_key = if_else(importance > 5, name, NA_character_)) %>%
  mutate(is_key = gsub('_.*', '', is_key))


plt_a <- my_graph %>%
  mutate(component = group_components()) %>%
  filter(component == which.max(as.numeric(table(component)))) %>%
  activate(edges) %>%
  ggraph( layout = "dendrogram", length = segments_changed) +
  geom_edge_elbow(edge_width = 0.5) +
  geom_node_point(aes(size= importance, colour = is_key)) +
  scale_size('Node Degree', guide = "legend", range = c(1, 6)) + 
  scale_colour_discrete('Nodes with Degree > 5', guide = "legend")+ 
  guides(col = guide_legend(nrow = 2, theme = theme(legend.byrow = TRUE)),
         size = guide_legend(nrow = 1, theme = theme(legend.byrow = TRUE ))) +
  annotate("rect", 
           ymin = c(0,2,4,6,8,10,12,14,16,18), 
           ymax =  c(1,3,5,7,9,11,13,15,17,19), 
           xmin = -Inf, xmax = Inf, alpha = 0.2, fill = "grey") +
  #scale_colour_brewer(palette = 'Set1') + 
  coord_flip() +
  scale_x_reverse() +
  scale_y_reverse() + 
  theme_void() +
  theme(legend.position = 'inside',
        legend.position.inside = c(0.95, 0),
        legend.justification = c(1, 0),
        legend.title.position = "top",
        legend.text = element_text(size = 9),
        legend.title = element_text(size = 10, face = 'bold'))


#all_colors <- brewer.pal(9, "Greys")
#color_subset <- colorRampPalette(all_colors)(100)[0:70]  # middle 50%


plt_b <-updated %>%
  dplyr::select(time_since_parent, parent_class, parent_label, cluster_label, cluster_region) %>% 
  filter(parent_label %in% (my_graph %>% 
                              as_tibble() %>%
                              filter(importance > 5)  %>%
                              pull(name))) %>% 
  ggplot()  +
  
  stat_function(
    fun = dlnorm,
    geom = "polygon",
    args = list(meanlog = -1.6267804, sdlog = 0.9225096),
    fill = "grey",
    alpha = 0.5
  )  +

 # geom_segment(aes( yend = Inf, y = -Inf, x = x, colour = col, alpha = col), 
 #              inherit.aes = F, 
               #data = data.frame(x = seq(0,2, by = 0.001), col = dlnorm(seq(0,2, by = 0.001), meanlog = -1.6267804, sdlog = 0.9225096))) +
  #scale_colour_gradientn(colours = color_subset) + 
  new_scale_colour()+
  
  geom_point(aes(y = parent_label, 
                 x = time_since_parent,
                 colour = parent_label), 
             #shape = 1,
             alpha = 0.5,
             size = 1.5) +
  
  scale_x_continuous('Inter-reassortant Interval', limits = c(0, 2), expand= c(0,0)) + 
  scale_y_discrete(labels = function(x) str_to_upper(x) %>% gsub('_.*', '', .)) + 

  theme_void() + 
  theme(panel.grid.major.y = element_line(colour = 'grey'),
        axis.text.y = element_text(size = 9), 
        axis.title.x = element_text(size = 10), 
        axis.line.x = element_line(), 
        axis.text.x = element_text(size = 9), 
        axis.ticks.x = element_line(), 
        plot.margin = unit(c(5.5, 5.5, 5.5, 5.5), "points"),
        legend.position = 'none') 

cowplot::plot_grid(plt_a,plt_b, nrow = 2, rel_heights = c(0.6, 0.4), labels = 'AUTO', label_size = 9)



# Segment Proportion Plots
plt_c <- segments_change %>%
  dplyr::select(ends_with('switch')) %>%
  pivot_longer(everything(), values_to = 'switch', names_to = 'segment') %>%
  mutate(segment = gsub('_switch', '', segment),
         switch = as.numeric(switch)) %>%
  summarise(switch = sum(switch), .by = segment) %>%
  mutate(segment = reorder(segment, switch)) %>%
  ggplot() +
  scale_x_discrete('Segment' ,expand = c(0,0)) + 
  scale_y_continuous('Reassortants (n)', expand = c(0,0), limits = c(0,115)) + 
  scale_fill_distiller(palette = 'YlOrRd', na.value = NA, direction = 1) +
  geom_bar(aes(x = segment, y = switch, fill = switch), stat = 'identity', colour = 'black') + 
  geom_bracket(
    xmin = "NP", xmax = "PB2", y.position = 110,
    label = "vRNP Complex", tip.length = c(0.02, 0.02)
  ) + 
  theme_classic() + 
  theme(axis.text = element_text(size = 9),
        axis.title = element_text(size = 10),
        legend.text = element_text(size = 9),
        legend.title = element_text(size = 10),
        legend.position = 'none')


plt_d <- observed_jaccards %>%
  ggplot() + 
  geom_point(aes(x = seg1, y = seg2, fill = jaccard, size = jaccard), shape = 21) +
  scale_size( range = c(0.5, 10)) + 
  #geom_tile(aes(x = seg1, y = seg2, fill = jaccard)) +
  scale_fill_distiller(palette = 'YlOrRd', na.value = NA, direction = 1) +
  scale_x_discrete('Segment' ) + 
  scale_y_discrete('Segment' ) + 
  theme_classic()+ 
  theme(axis.text = element_text(size = 9),
        axis.title = element_text(size = 10),
        legend.text = element_text(size = 9),
        legend.title = element_text(size = 10),
        legend.position = 'none')



# Combine
ggdraw() +
  draw_plot(plt_a, x = 0, y = 0.5, width = 1, height = .5) +
  draw_plot(plt_b, x = 0, y = 0.25, width = 1, height = .25) +
  draw_plot(plt_c, x = 0, y = 0, width = 0.5, height = 0.25) +
  draw_plot(plt_d, x = 0.5, y = 0, width = 0.5, height = 0.25) +
  draw_plot_label(label = c("A", "B", "C", 'D'), size = 12,
                  x = c(0, 0, 0, 0.5), y = c(1, 0.5, 0.25, 0.25))
ggsave( '~/Downloads/flu_plots/figure_network_time.jpeg', height = 15, width = 20, units = 'cm')
