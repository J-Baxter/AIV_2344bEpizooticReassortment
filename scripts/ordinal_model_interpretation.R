####################################################################################################
####################################################################################################
## Script name: Ordinal Model Interpretation
##
## Purpose of script:
##
## Date created: 2025-05-13
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
library(brms)
library(broom)
library(broom.mixed)
library(tidybayes) # missing
library(bayesplot)
library(emmeans) # missing
library(marginaleffects) # missing
library(magrittr)
library(ggmcmc) # missing
library(bayestestR)
library(modelbased)


# User functions


############################################## DATA ################################################
# average posterior predictions 
avg_predictions(ordinal_model)
avg_predictions(ordinal_model, by = 'cluster_region') # empirical distribution (ie mean of predictions)


#percentage point change of reassortant class ~ previous class (check this actuall is)
avg_comparisons(ordinal_model, variables = list("parent_class" = 'revpairwise'), newdata = 'balanced')


# average prediction of time since last major (1, 2, 3yr interals)
avg_predictions(ordinal_model, variables = list('time_since_last_major' = c(0.5,1,3,5)))

# slops of time since last major (1, 2, 3yr interals)
ordinal_model %>%
  avg_comparisons(variables = list('time_since_last_major' = 1), 
                  by = "time_since_last_major",
                  newdata = datagrid(time_since_last_major = c(0.5, 1, 3), grid_type = 'counterfactual'), 
                  type = 'response')



# average prediction for 1,3 segmetns changing
avg_predictions(ordinal_model, variables = list('segments_changed' = c(1,2,4)))


# effect of eac hadditional segment
ordinal_model %>%
avg_comparisons(variables = list('segments_changed' = 1), 
                by = "segments_changed",
                type = 'response')


# contrasts for continent, all else equal
avg_comparisons(ordinal_model, variables = list("cluster_region" = 'pairwise'), newdata = 'balanced')

############################################## MAIN ################################################

plt_d <- avg_predictions(ordinal_model, by = 'cluster_region') %>%
  get_draws(shape = "rvar") |>
  ggplot(aes(y = cluster_region, xdist = rvar, slab_colour = group,
             slab_fill = group)) + 
  stat_slabinterval(position = 'dodgejust',
                    point_interval = "median_hdci",
                    slab_alpha = 0.7 ,
                    p_limits = c(0.01, 0.99),
                    normalize = 'xy',
                    .width = 0.95) + 
  scale_fill_brewer(palette = 'Set1', aesthetics = 'slab_fill') + 
  scale_colour_brewer(palette = 'Set1', aesthetics = 'slab_colour') + 
  global_theme



new_data = class_data %>% 
  dplyr::select(collection_regionname, previous_class) %>% 
  filter(previous_class != 'none') %>%
  distinct() %>% 
  mutate(time_since_previous = mean(class_data$time_since_previous))

epreds <- ordinal_model %>%
  epred_draws(newdata = new_data) %>%
  median_hdci() %>%
  mutate(across(c(.epred, .lower, .upper), ~ ifelse(.epred < 0.005, NA_real_, .x))) %>%
  mutate(label =  case_when(!is.na(.epred) ~ paste0(round(.epred, digits = 3), 
                                                    ' (', round(.lower, digits = 3),
                                                    '-', round(.upper, digits = 3), ')'))) %>%
  mutate(across(c(.category, previous_class), .fns = ~ str_to_title(.x)))

graph <- avg_predictions(ordinal_model, by = 'parent_class') %>%
  as_tibble() %>%
  mutate(across(c(estimate, conf.low, conf.high), ~ ifelse(estimate < 0.005, NA_real_, .x))) %>%
  mutate(label =  case_when(!is.na(estimate) ~ paste0(round(estimate, digits = 3), 
                                                      ' (', round(conf.low, digits = 3),
                                                      '-', round(conf.high, digits = 3), ')'))) %>%
  mutate(across(c(group, parent_class), .fns = ~ str_to_title(.x))) %>%
  
  rename(from = parent_class, 
         to = group) %>% 
  drop_na(estimate) %>%
  relocate(from, to) %>%
 # mutate(across(c(from, to), .fns = ~paste0(.x, "_", collection_regionname))) %>%  # Create unique node IDs for each region
  as_tbl_graph(directed = TRUE) %>% 
  activate(nodes) %>%
  
  mutate(name = gsub('_.*', '', name))   %>%
  
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


plt_e <- layout %>%
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
                     edge_width = estimate ),
                 label_size = 3) +
  geom_edge_loop(arrow = arrow(length = unit(4, 'mm')), 
                 end_cap = circle(12, 'mm'),
                 start_cap =  circle(12, 'mm'),
                 vjust	= 2,
                 aes(#span = self_moderate_span,
                     #direction = self_moderate_direction,
                     strength = 1,
                     label = label,
                     edge_width = estimate),
                 label_size = 3) +
  geom_node_circle(aes(x0 = x, y0 = y, r = 0.16, fill = name, colour = name)) +
  
  #geom_node_point(size = 15, aes(colour = collection_regionname, alpha = name)) +
  geom_node_text(aes(label = lab),
                 colour = 'white',
                 size = 6,
                 fontface = 'bold') +
  scale_fill_brewer(palette = 'Set1') + 
  scale_edge_width_continuous(range = c(0.1, 1.5)) + 
  theme_void() + 
  theme(legend.position = 'none', 
        strip.background = element_blank(),
        strip.text = element_text(size = 10, face = 'bold')) +
  scale_x_continuous(expand = c(0,0), limits = c(0,3.5)) + 
  scale_y_continuous(expand = c(0,0), limits = c(0.5,2.5))





plt_f <- avg_predictions(ordinal_model, by = 'segments_changed') %>%
  as_tibble() %>%
  ggplot(aes(ymin = conf.low, ymax = conf.high, y = estimate, x = segments_changed, fill = group,
             colour = group)) + 

  geom_ribbon(aes(colour = NULL), alpha = 0.2) + 
  geom_line() + 
  
  scale_fill_brewer(palette = 'Set1') + 
  scale_colour_brewer(palette = 'Set1') + 
  global_theme


plt_g <- conditional_effects(ordinal_model, "time_since_last_major", categorical = TRUE) %>%plot()


cowplot::plot_grid(plt_d, plt_e, plt_f,plt_g, labels = 'AUTO', ncol = 2, align = 'hv')
############################################## WRITE ###############################################




############################################## END #################################################
####################################################################################################
####################################################################################################