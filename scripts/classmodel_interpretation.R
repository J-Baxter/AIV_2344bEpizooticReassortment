####################################################################################################
####################################################################################################
## Script name:
##
## Purpose of script:
##
## Date created: 2025-03-28
##
##
########################################## SYSTEM OPTIONS ##########################################
options(scipen = 6, digits = 7) 
memory.limit(30000000) 

  
########################################## DEPENDENCIES ############################################
# Packages
library(tidyverse)
library(magrittr)
library(ggraph)
library(tidygraph)


# User functions

# User functions
scientific_10 <- function(x) {
  parse(text=gsub("e\\+*", " %*% 10^", scales::scientific_format()(x)))
  }

############################################## DATA ################################################
host_colours <- c(
  'anseriformes-domestic' = '#a6cee3',
  'anseriformes-wild' = '#1f78b4',
  'galliformes-domestic' = '#b2df8a',
  'galliformes-wild' = '#33a02c',
  'mammal' = '#fb9a99',
  'human' = '#e31a1c',
  'charadriiformes-wild' = '#fdbf6f',
  'other-bird' = '#ff7f00',
  'unknown' = '#cab2d6',
  'environment' = '#6a3d9a')


region_colours <- c('europe' = '#1b9e77',
                    'asia' ='#d95f02',
                    'africa' ='#7570b3',
                    'australasia' = '#e7298a',
                    'central & northern america' ='#66a61e',
                    'south america' ='#e6ab02')


global_theme <- theme_classic()+
  theme(
    #text = element_text(size=10),
    axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0), size = 10),
    axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0), size = 10),
    axis.text = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0), size = 8),
    legend.text = element_text(size = 8),
    legend.position = 'none', 
    panel.spacing = unit(2, "lines"), 
    strip.background = element_blank()
  )


############################################## MAIN ################################################

# Data to simulate from
# expectation of the posterior
class_model %>%
  epred_draws(newdata = new_data) %>%
  #filter(.epred >0.01) %>%
  ggplot(aes(x  = .epred, y = collection_regionname,
             slab_colour = collection_regionname,
             slab_fill = collection_regionname)) +
  stat_halfeye(p_limits = c(0.001, 0.999),
               point_interval = "median_hdci",
               # expand = T,
               #density = 'histogram',
               .width =  0.95,
               normalize = "xy"
  ) +
  facet_grid(cols = vars(.category)) 



###### Proportions (full posterior draw) ############
class_model %>%
  predicted_draws(newdata = new_data %>% filter(previous_class == 'minor')) %>%
  count(predicted_class = .prediction) %>%
  mutate(prop = n / sum(n)) %>%
  ggplot(aes(x = factor(predicted_class), y = n / sum(n))) +
  geom_col() + 
  facet_grid(cols = vars(collection_regionname))


avg_predictions(class_model, by = 'collection_regionname')
avg_predictions(class_model, by = 'collection_regionname', newdata = 'balanced')

avg_predictions(class_model,
                variables = list("time_since_previous" = 0.08),
                by = "time_since_previous",
                newdata = datagrid(time_since_previous = c(0.08, 0.25),
                                   grid_type = 'counterfactual'))

avg_predictions(class_model, by = 'collection_regionname', newdata = 'balanced')




class_model |> 
  avg_comparisons(variables = list("time_since_previous" = 0.08), by = "time_since_previous",
                  newdata = datagrid(time_since_previous = c(0.08, 0.25), grid_type = 'counterfactual'),  type = 'response')


avg_predictions(class_model, by = 'previous_class', newdata = datagrid(time_since_previous = 0.08, 
                                                                   previous_class=unique(class_data$previous_class),
                                                                   grid_type = 'counterfactual'))


plot_predictions(class_model, condition = c("time_since_previous", "collection_regionname")) + facet_wrap(.~group)
###### Expectation of the posterior (linear predictors) ############
# collection_regionname ~ previous class, at the median time between reassortants
#
new_data = class_data %>% 
  dplyr::select(collection_regionname, previous_class) %>% 
  filter(previous_class != 'none') %>%
  distinct() %>% 
  mutate(time_since_previous = mean(class_data$time_since_previous))

epreds <- class_model %>%
  epred_draws(newdata = new_data) %>%
  median_hdci() %>%
  mutate(across(c(.epred, .lower, .upper), ~ ifelse(.epred < 0.005, NA_real_, .x))) %>%
  mutate(label =  case_when(!is.na(.epred) ~ paste0(round(.epred, digits = 3), 
                                                    ' (', round(.lower, digits = 3),
                                                    '-', round(.upper, digits = 3), ')'))) %>%
  mutate(across(c(.category, previous_class), .fns = ~ str_to_title(.x)))



graph <- epreds %>%
  rename(from = previous_class, 
         to = .category) %>% 
  drop_na(.epred) %>%
  relocate(from, to) %>%
  mutate(across(c(from, to), .fns = ~paste0(.x, "_", collection_regionname))) %>%  # Create unique node IDs for each region
  as_tbl_graph(directed = TRUE) %>% 
  activate(nodes) %>%
  
  mutate(collection_regionname = gsub('.*_', '', name),
         name = gsub('_.*', '', name))   %>%
  
  mutate(lab = case_when(name == 'Major'~ 'Mj',
                         name == 'Moderate' ~ "Md",
                         name == 'Minor' ~ 'Mn'))


layout <- create_layout(graph, layout = "kk") %>%
  mutate(x = case_when(name == 'Minor' ~ 1,
                       name == 'Moderate' ~ 2,
                       name == 'Major'  ~ 3),
         y = case_when(name == 'Minor' ~ 1.5,
                       name == 'Moderate' ~ 2, 
                       name == 'Major'  ~ 1))

nodes <- layout %>% 
  dplyr::select(x,y,name) %>%
  pivot_longer(cols = c(x,y), values_to = 'coord', names_to = 'axis') %>% 
  pull(name)

self_moderate_span <- ifelse(nodes == "Minor", 90, -90)
self_moderate_direction <- ifelse(nodes == "Moderate", 45, 270)


layout %>%
  # activate(nodes) %>%
  # left_join(layout, by = "name") %>%
  ggraph() + 
  geom_edge_bend(arrow = arrow(length = unit(4, 'mm')), 
                 strength = 0.3, 
                 start_cap = circle(12, 'mm'),
                 end_cap = circle(12, 'mm'),
                 angle_calc = 'along',
                 label_dodge	= unit(3, 'mm'),
                 aes(label = label,
                     edge_width = .epred ),
                 label_size = 3) +
  geom_edge_loop(arrow = arrow(length = unit(4, 'mm')), 
                 end_cap = circle(12, 'mm'),
                 start_cap =  circle(12, 'mm'),
                 vjust	= 2,
                 aes(span = self_moderate_span,
                     direction = self_moderate_direction,
                     strength = 1,
                     label = label,
                     edge_width = .epred),
                 label_size = 3) +
  geom_node_circle(aes(x0 = x, y0 = y, r = 0.16, fill = collection_regionname, colour = collection_regionname, alpha = name)) +
  
  #geom_node_point(size = 15, aes(colour = collection_regionname, alpha = name)) +
  geom_node_text(aes(label = lab),
                 colour = 'white',
                 size = 6,
                 fontface = 'bold') +
  facet_grid(cols = vars(collection_regionname),
             labeller = labeller(collection_regionname= str_to_title)) + 
  scale_colour_manual(values = region_colours) + 
  scale_fill_manual(values = region_colours) + 
  
  scale_alpha_manual( values = c('Major' = 1, 'Moderate' = 0.7, 'Minor' = 0.5),
                      labels = str_to_title) + 
  scale_edge_width_continuous(range = c(0.1, 1.5)) + 
  theme_void() + 
  theme(legend.position = 'none', 
        strip.background = element_blank(),
        strip.text = element_text(size = 10, face = 'bold')) +
  scale_x_continuous(expand = c(0,0), limits = c(0,3.5)) + 
  scale_y_continuous(expand = c(0,0), limits = c(0.5,2.5))



###### Conditional effect of Continent on the probability of a 'major' reassortant' ############
group_effects <- class_model %>%
  avg_comparisons(newdata =  ungroup(new_data),
                  variables = list('collection_regionname' = 'pairwise') , type = 'response') %>%
  get_draws() %>%
  as_tibble()


group_effects_hdci <- group_effects %>%
  group_by(term, group, contrast) %>%
  median_hdci(draw, .exclude = c('drawid')) %>%
  rename(median = draw) %>%
  as_tibble()



group_effects %>% 
  filter(contrast == 'europe - central & northern america') %>%
  mutate(group = fct_relevel(group, c('Minor', 'Moderate', 'Major'))) %>%
  left_join(group_effects_hdci) %>%
  ggplot(aes(x = draw,
             y = str_to_title(group),
             slab_colour = median,
             slab_fill = median)) +
  stat_halfeye(slab_alpha = 0.7,
               p_limits = c(0.05, 0.95),
               point_interval = "median_hdci",
               linewidth = 1.5,
               .width =  0.95) +
  geom_vline(xintercept = 0, linetype = 'dashed') + 
  scale_x_continuous(expand = c(0,0),
                     limits = c(-0.6, 0.6),
                     'Percentage Point Change') + 
  scale_y_discrete('Reassortant Class') + 
  scale_fill_distiller(palette = 'RdYlBu', aesthetics = 'slab_fill') +
  scale_colour_distiller(palette = 'RdYlBu', aesthetics = 'slab_colour') +
  theme_classic()




 


############################################## WRITE ###############################################




############################################## END #################################################
####################################################################################################
####################################################################################################