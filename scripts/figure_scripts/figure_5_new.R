####################################################################################################
####################################################################################################
## Script name: 
##
## Purpose of script: 
##
## Date created: 2025-05-17
##
##
########################################## SYSTEM OPTIONS ##########################################
options(scipen = 6, digits = 7) 
memory.limit(30000000) 


########################################## DEPENDENCIES ############################################
# Packages
library(tidyverse)
library(magrittr)

############################################## DATA ################################################

reassortant_ancestral_changes <- read_csv('./reassortant_ancestral_changes.csv')


#e41a1c
#377eb8
#4daf4a
class_colours <- c('major' = '#e41a1c',
                   'moderate' = '#377eb8',
                     'minor' = '#4daf4a')
region_colours <- c('europe' = '#1b9e77',
                    'asia' ='#d95f02',
                    'africa' ='#7570b3',
                    'australasia' = '#e7298a',
                    'central & northern america' ='#66a61e',
                    'south america' ='#e6ab02')


############################################## MAIN ################################################
test_pred %>%
  group_by(.draw, collection_regionname) %>% 
  summarise(avg_epred = mean(.epred), .groups = "drop") %>% 
  group_by(collection_regionname, .epred) %>%
  median_hdci(avg_epred)

test_pred %>%
  ggplot() +
  geom_histogram(aes(x = .epred, fill = collection_regionname, colour = collection_regionname),
                 binwidth = 0.2,
                 alpha = 0.7) + 
  scale_colour_manual(values = region_colours)+
  scale_fill_manual(values = region_colours) + 
  
  facet_grid(
    cols = vars(collection_regionname),
    labeller =  labeller(collection_regionname=str_to_title),
    scales = 'free_y') +
  scale_y_continuous('Probability Density' ,
                     breaks = seq(0,1.5,by=0.5),
                     labels = seq(0,1.5,by=0.5),
                     expand = c(0,0)) +
  #geom_vline(aes(xintercept = emmean, colour = collection_regionname), data = averages, linetype = 'dashed') +
  #geom_text(aes(label =  paste0("E*'('*X*'|'*X*'>'*0*') = '*", label, "~km^2"), 
  #      colour = collection_regionname),
  #   parse = T,
  # x = 17.5, 
  #  y = 0.6,
  # size = 2.5,
  # data = averages) + 
  global_theme + 
  theme(strip.placement  = 'inside',
        strip.text = element_text(face = 'bold', size = 10),
        strip.background = element_blank(),
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 9),
        legend.text = element_text(size = 8))



####### Reassortant Class Plots #########
### D 
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
              select(name = cluster_label, cluster_class)) 

plt_5d <- ggraph(my_graph %>%
         mutate(component = group_components()) %>%
         filter(component == which.max(as.numeric(table(component)))),  layout = "fr", weights = segments_changed) +
  geom_edge_link() +
  geom_node_point(aes(color = cluster_class, size=10 )) +
  geom_node_label(aes(label = ifelse(cluster_class == 'major', gsub('_.*', '', name), '')), repel = TRUE) +
  scale_size(guide = NULL) + 
  scale_colour_brewer(palette = 'Set1', 'Cluster Class') + 
  theme_void() + 
  theme(legend.text = element_text(size = 8),
        legend.position = 'inside',
        legend.position.inside = c(0.1,0.9))

colour = .cluster
### E (Class probabilities)
plt_5e <-avg_predictions(ordinal_model, by = 'cluster_region') %>%
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
  
  scale_x_discrete('Continent', expand = c(0,0), labels = str_to_title) + 
  global_theme + 
  theme(strip.placement  = 'inside',
        strip.text = element_text(face = 'bold', size = 10),
        strip.background = element_blank(),
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 9),
        legend.text = element_text(size = 8))



### F (Segment)
plt_5f <-avg_predictions(ordinal_model, by = 'segments_changed', newdata = 'balanced') %>%
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
  
  scale_x_continuous('Segments Changed (N)', expand = c(0,0)) + 
  global_theme + 
  theme(strip.placement  = 'inside',
        strip.text = element_text(face = 'bold', size = 10),
        strip.background = element_blank(),
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 9),
        legend.text = element_text(size = 8))


### G (Time)
plt_5g <-avg_predictions(ordinal_model, by = 'time_since_last_major', newdata = 'balanced') %>%
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
  
  scale_x_continuous('Interval since last Major', expand = c(0,0)) + 
  global_theme + 
  theme(strip.placement  = 'inside',
        strip.text = element_text(face = 'bold', size = 10),
        strip.background = element_blank(),
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 9),
        legend.text = element_text(size = 8))


plot_grid(plt_5d, plt_5e, plt_5f, plt_5g, ncol = 2, align = 'hv')
############################################## WRITE ###############################################


############################################## END #################################################
####################################################################################################
####################################################################################################p