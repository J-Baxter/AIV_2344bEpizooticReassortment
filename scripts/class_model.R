####################################################################################################
####################################################################################################
## Script name: Class Model
##
## Purpose of script: to model the probability that any given reassortant belongs to class X,
## stratified by continent
##
## Date created: 2025-03-24
##
##
########################################## SYSTEM OPTIONS ##########################################
options(scipen = 6, digits = 7) 
memory.limit(30000000) 

  
########################################## DEPENDENCIES ############################################
# Packages
library(tidyverse)
library(magrittr)


# User functions

# variables required
# - class of reassortant /
# - immediately preceeding reassortant /
# - time since last ressortant /
# - continent /
# (- diversity)
# - proxy incidence
# - number of genomes


############################################## DATA ################################################
combined_data <- read_csv('./2024Aug18/treedata_extractions/2024-09-20_combined_data.csv')
summary_data <- read_csv('./2024Aug18/treedata_extractions/summary_reassortant_metadata_20240904.csv') %>%
  dplyr::select(-c(cluster_label,
                   clade)) 

final_clusters <- read_csv('./final_clusters_2025Mar21.csv')


############################################## MAIN ################################################
class_data <- combined_data %>% 
  dplyr::select(-c(group2, col2)) %>%
  left_join(final_clusters) %>%
  
  dplyr::select(cluster_profile, 
                segment,
                TMRCA,
                collection_regionname,
                class) %>%
  
  mutate(collection_regionname = case_when(grepl('europe', collection_regionname) ~ 'europe',
                                           grepl('africa', collection_regionname) ~ 'africa',
                                           grepl('asia', collection_regionname) ~ 'asia',
                                           grepl('(central|northern) america', collection_regionname) ~ 'central & northern america',
                                           .default = collection_regionname
  )) %>%
  
  mutate(TMRCA = date_decimal(TMRCA) %>%
           round_date(unit = 'day')%>%
           as_date(.)) %>%
  
  filter(segment== 'ha') %>%
  arrange(TMRCA) %>%
  group_by(collection_regionname) %>%
  
  # class of most recent prior reassortant
  mutate(previous_class = lag(class),
         previous_class = replace_na(previous_class, 'none')) %>%
  
  # time since last reassortant
  # calculate as interval, then convert to decimal of years
  mutate(time_since_previous = interval(lag(TMRCA), TMRCA) %>%
           as.numeric('years')) %>%
  
  # revert formatting of dates for model
  # TMRCA to decimal, then subtract 2016 for scale
  mutate(tmrca_decimal = decimal_date(TMRCA),
         tmrca_decimal_adjusted = tmrca_decimal - 2016) %>%
  
  # select vars of interest
  dplyr::select(collection_regionname,
                class, 
                previous_class,
                tmrca_decimal_adjusted,
                time_since_previous) %>%
  
  # exclude the first reassortant from each continent
  drop_na(time_since_previous) %>%
  drop_na(class) %>%
  
  filter(collection_regionname != 'africa')


# Plot Data
ggplot(class_data) + geom_bar(aes(x = class, fill = previous_class), position = 'stack') + facet_wrap( ~ collection_regionname)
ggplot(class_data, aes(x = time_since_previous)) + geom_histogram() + 
  #geom_smooth(method = "lm", se = FALSE) +
  facet_wrap( class ~ collection_regionname)


# Formula
class_formula_priors <- get_prior(class ~  previous_class:collection_regionname + ,
                                  data = class_data,
                                  family = categorical(link ='logit')) 



# Priors
class_formula_priors$prior[c(1:13)] <- "normal(0,5)"
class_formula_priors$prior[c(16:27)] <- "normal(0,5)"


# Model

basic_model <- brm(
  class ~ previous_class:collection_regionname,
  prior = class_formula_priors,
  data = class_data,
  family = categorical(link ='logit'),
  chains = 2, iter = 10000, cores = 2
)


# Posterior Predictive
pp_check(basic_model, ndraws = 100, type = 'bars') 


# Data to simulate from
# expectation of the posterior
basic_model %>%
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
basic_model %>%
  predicted_draws(newdata = new_data %>% filter(previous_class == 'minor')) %>%
  count(predicted_class = .prediction) %>%
  mutate(prop = n / sum(n)) %>%
  ggplot(aes(x = factor(predicted_class), y = n / sum(n))) +
  geom_col() + 
  facet_grid(cols = vars(collection_regionname))



###### Expectation of the posterior (linear predictors) ############
# collection_regionname ~ previous class, at the median time between reassortants
#
new_data = class_data %>% 
  dplyr::select(collection_regionname, previous_class) %>% 
  filter(previous_class != 'none') %>%
  distinct() %>% 
  mutate(time_since_previous = mean(class_data$time_since_previous))

epreds <- basic_model %>%
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
group_effects <- basic_model %>%
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
