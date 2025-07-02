### D 
# Plot network using tidygraph
class_colours <- c('major'=        '#A92523'
                   ,'moderate'= '#FA9F42',
                   'minor' = '#2B4162')


graph_palette <- colorRampPalette(c('#2D00F7', '#6A00F4', '#8900F2', '#A100F2', '#B100E8', '#BC00DD', '#D100D1', '#DB00B6', '#E500A4','#F20089'))(100 )

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
  
  mutate(is_key = if_else(importance > 5, name, NA_character_))

my_graph %>%
  as_tibble() %>% filter(importance >0)%>%
  ggplot(aes(x = importance)) +
  geom_histogram(binwidth = 1, aes(y = after_stat(density))) + 
  geom_function(
    fun = dllogis,
    args = list(shape = 2.106323, scale = 1.587460),
    colour = "red"
  )  scale_x_continuous('Descendant Reassortants') 

out_deg <- my_graph %>% as_tibble() %>% filter(importance >0) %>%pull(importance) 
out_deg.ln <- fitdist(out_deg, "lnorm")
out_deg.exp <- fitdist(out_deg, "exp")

out_deg.ll <- fitdist(out_deg, "llogis", start = list(shape = 1, scale = 500))
out_deg_pareto <- fitdist(out_deg, "pareto", start = list(shape = 1, scale = 500))
cdfcomp(list(out_deg.ln, out_deg.ll, out_deg_pareto, out_deg.exp), xlogscale = TRUE, 
        ylogscale = TRUE, legendtext = c("lognormal", "loglogistic", "Pareto", 'Exp'))


plt_3a <- ggraph(my_graph %>%
                  mutate(component = group_components()) %>%
                  filter(component == which.max(as.numeric(table(component)))),  layout = 'kk') +
  geom_edge_link() +
  geom_node_point(aes(size= importance )) +
  #scale_colour_manual() + 
  coord_flip() +
  scale_x_reverse() +
  scale_y_reverse() +
  #geom_node_label(aes(label = ifelse(cluster_class == 'major', name, '')), repel = TRUE) +
  #scale_colour_brewer(palette = 'Set1') +
  #scale_colour_manual(values = class_colours, 'Cluster Class',  labels = str_to_title,) +
  #scale_colour_distiller(palette = 'Greens', direction = -1) + 
  #scale_colour_gradientn(colours = graph_palette )+
  scale_size_continuous('Offspring') + 
  theme_void() +

  theme(strip.placement  = 'inside',
        strip.text = element_text(face = 'bold', size = 10),
        strip.background = element_blank(),
        legend.title = element_text(size = 9),
        legend.text = element_text(size = 8),
        plot.margin = unit(c(0,0,1,0), "cm"),
        legend.position = 'bottom',
        legend.box = "vertical")


ggraph(my_graph %>%
         mutate(component = group_components()) %>%
         filter(component == which.max(as.numeric(table(component)))) %>% 
         activate(edges) %>%
         rename(length = segments_changed),  'unrooted', length = length) +
  geom_edge_link() +
  geom_node_point(aes(color = as.factor(cluster_region), size= importance )) +
  #geom_node_label(aes(label = ifelse(cluster_class == 'major', name, '')), repel = TRUE) +
  #scale_colour_brewer(palette = 'Set1') +
  #scale_colour_manual(values = class_colours, 'Cluster Class',  labels = str_to_title,) +
  scale_colour_manual(values = region_colours) + 
  scale_size_continuous('Offspring') + 
  theme_void() +
  
  theme(strip.placement  = 'inside',
        strip.text = element_text(face = 'bold', size = 10),
        strip.background = element_blank(),
        legend.title = element_text(size = 9),
        legend.text = element_text(size = 8),
        plot.margin = unit(c(0,0,1,0), "cm"),
        legend.position = 'bottom',
        legend.box = "vertical")

### E (Class probabilities)
plt_3b <-avg_predictions(ordinal_model, by = 'cluster_region') %>%
  get_draws(shape = "rvar") |>
  #filter(group == 'major') %>%
  ggplot(aes(x = cluster_region, ydist = rvar, colour = group)) + 
  stat_interval(aes(color_ramp = after_stat(level)), 
                position = position_dodgejust(width = 0.5, justification = 0),
                point_interval = "median_hdci")+
  
  scale_colour_manual(values = class_colours) + 
  
  scale_y_continuous('Pr(Class = k)' ,
                     breaks = seq(0,1,by=0.25),
                     labels = seq(0,1,by=0.25),
                     expand = c(0,0)) +
  
  scale_x_discrete('Continent', expand = c(0,0), labels =  function(x) str_to_title(x) ) + 
  global_theme + 
  theme(strip.placement  = 'inside',
        strip.text = element_text(face = 'bold', size = 10),
        strip.background = element_blank(),
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 9),
        legend.text = element_text(size = 8))



### F (Segment)
plt_3c <-avg_predictions(ordinal_model,  variables = list('segments_changed' = 1:7)) %>%
  as_tibble() %>%
  ggplot(aes(ymin = conf.low, ymax = conf.high, y = estimate, x = segments_changed, fill = group,
             colour = group)) + 
  
  geom_ribbon(aes(colour = NULL), alpha = 0.2) + 
  geom_line(linewidth = 1) + 
  
  scale_fill_manual(values = class_colours) + 
  scale_colour_manual(values = class_colours) + 
  
  scale_y_continuous('Pr(Class = k)' ,
                     breaks = seq(0,1,by=0.25),
                     labels = seq(0,1,by=0.25),
                     expand = c(0,0)) +
  
  scale_x_continuous('Segments changed from previous (N)', expand = c(0,0)) + 
  global_theme + 
  theme(strip.placement  = 'inside',
        strip.text = element_text(face = 'bold', size = 10),
        strip.background = element_blank(),
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 9),
        legend.text = element_text(size = 8))


### G (Time)



ggdraw() +
  draw_plot(plt_3a, x = 0, y = .5, width = 1, height = .5) +
  draw_plot(plt_3b, x = 0, y = 0, width = .5, height = .5) +
  draw_plot(plt_3c, x = 0.5, y = 0, width = 0.5, height = 0.5) +
  draw_plot_label(label = c("A", "B", "C"), size = 10,
                  x = c(0, 0, 0.5), y = c(1, 0.5, 0.5))



ggsave('~/Downloads/flu_plots/figure_clusters.jpeg', height = 18, width = 25, units = 'cm', dpi = 360)




plt_supa <- avg_predictions(ordinal_model, variables = list('time_since_last_major' = seq(0, 2, by  = 0.25))) %>%
  as_tibble() %>%
  ggplot(aes(ymin = conf.low, ymax = conf.high, y = estimate, x = time_since_last_major, fill = group,
             colour = group)) + 
  
  geom_ribbon(aes(colour = NULL), alpha = 0.2) + 
  geom_line(linewidth = 1) + 
  
  scale_fill_manual(values = class_colours) + 
  scale_colour_manual(values = class_colours) + 
  scale_y_continuous('Pr(Class = k)' ,
                     breaks = seq(0,1,by=0.25),
                     labels = seq(0,1,by=0.25),
                     expand = c(0,0)) +
  
  scale_x_continuous('Interval since last major reassortant', expand = c(0,0)) + 
  global_theme + 
  theme(strip.placement  = 'inside',
        strip.text = element_text(face = 'bold', size = 10),
        strip.background = element_blank(),
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 9),
        legend.text = element_text(size = 8))





# Generate a truncated version of the palette
all_colors <- brewer.pal(9, "Reds")
color_subset <- colorRampPalette(all_colors)(100)[0:70]  # middle 50%

plt_supb <- updated %>%
  dplyr::select(time_since_parent, parent_class, parent_label, cluster_label, cluster_region) %>% 
  filter(parent_label %in% (my_graph %>% 
                              as_tibble() %>%
                              filter(importance > 5)  %>%
                              pull(name))) %>% 
  
  ggplot()  +
  geom_segment(aes( yend = Inf, y = -Inf, x = x, colour = col, alpha = col), 
               inherit.aes = F, 
               data = data.frame(x = seq(0,2, by = 0.001), col = dlnorm(seq(0,2, by = 0.001), meanlog = -1.6267804, sdlog = 0.9225096))) +
  geom_point(aes(y = parent_label, 
                 x = time_since_parent), 
             shape = 21,
             size = 2.5) +
  
  scale_x_continuous('Inter-reassortant Interval', limits = c(0, 2), expand= c(0,0)) + 
  scale_y_discrete(labels = function(x) str_to_upper(x) %>% gsub('_.*', '', .)) + 
  scale_colour_gradientn(colours = color_subset) + 
  

  
  #facet_wrap(~parent_class,
             #nrow = 2,
            # scales = 'free_y', 
             #space = 'free_y',
             #switch = 'y',
             #labeller = as_labeller(str_to_title)
#  ) +
  theme_void() + 
  theme(panel.grid.major.y = element_line(colour = 'grey'),
        axis.text.y = element_text(size = 8), 
        axis.title.x = element_text(size = 9), 
        axis.line.x = element_line(), 
        axis.text.x = element_text(size = 8), 
        axis.ticks.x = element_line(), 
        strip.placement  = 'outside',
        plot.margin = unit(c(5.5, 5.5, 5.5, 5.5), "points"),
        strip.text = element_text(face = 'bold', size = 9, angle = 90),
        legend.position = 'none') 


# All change ~ Number of offspring
updated %>% 
  left_join(diffusion_data %>% dplyr::select(starts_with('path_prop'), cluster_profile, persistence)) %>% 

  
  
  mutate(across(starts_with('path_prop'), .fns = ~replace_na(.x, 0))) %>% 
  rowwise() %>%
  mutate(max_host = which.max(c_across(path_prop_anseriformes_wild:path_prop_galliformes_domestic))) %>% 
  ungroup() %>% left_join(dplyr::select(., parent_profile = cluster_profile, parent_host = max_host, parent_persist = persistence)) %>% 
  ggplot(aes(x = parent_persist, y = time_since_parent)) + geom_point()
  mutate(all_change = case_when(max_host != parent_host | region_changed_from_previous > 0~ 'change', .default = 'no_change')) %>%
  left_join(as_tibble(my_graph) %>% dplyr::select(cluster_label = name, importance)) %>% 
  #dplyr::select(all_change, importance) %>%
  ggplot(aes(x = as.factor(max_host != parent_host ), y = importance)) + geom_boxplot()

cowplot::plot_grid(plt_3a,plt_supb, nrow = 2)

ggsave('~/Downloads/flu_plots/figure_network_time.jpeg', height = 30, width = 25, units = 'cm', dpi = 360)
