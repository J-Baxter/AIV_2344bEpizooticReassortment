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
library(ggdist)
library(tidygraph)
library(ggraph)
library(tidybayes)
library(cowplot)

############################################## DATA ################################################

reassortant_ancestral_changes <- read_csv('./reassortant_ancestral_changes.csv')

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

# A - Posterior distribution of latent counts 
test_pred <- as_draws_df(numbers_model_2$draws('N_rep') ) %>%
  pivot_longer(cols = starts_with("N_rep["), names_to = "row", values_to = ".epred") %>%
  mutate(row = as.integer(str_extract(row, "\\d+"))) %>%
  left_join(data_processed_2 %>% rowid_to_column('row'),
            by = 'row')

#test_pred %>%
 # group_by(.draw, collection_regionname) %>% 
  #summarise(avg_epred = mean(.epred), .groups = "drop") %>% 
  #group_by(collection_regionname, .epred) %>%
  #median_hdci(avg_epred)

plt_5a <- test_pred %>%
  ggplot() +
  geom_histogram(aes(x = .epred, fill = collection_regionname, colour = collection_regionname, y = after_stat(density)),
                 binwidth = 1, center = 0.5,
                 alpha = 0.7) + 
  scale_colour_manual(values = region_colours)+
  scale_fill_manual(values = region_colours) + 
  
  facet_grid(
    cols = vars(collection_regionname),
    labeller =  labeller(collection_regionname=str_to_title),
    scales = 'free_y') +
  scale_y_continuous('Probability Density' ,
                     breaks = seq(0,1,by=0.25),
                     labels = seq(0,1,by=0.25),
                     limits = c(0,1),
                     expand = c(0,0)) +
  scale_x_continuous('Reassortants (N)', expand= c(0,0) ,limits = c(0, 10),   breaks = seq(0,10,by=2),
                     labels = seq(0,10,by=2)) + 
  global_theme + 
  theme(strip.placement  = 'inside',
        strip.text = element_text(face = 'bold', size = 10),
        strip.background = element_blank(),
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 9),
        legend.text = element_text(size = 8))


# B - Posterior distribution latent counts
thetas <- c('continent_specific_theta[1]',
             'continent_specific_theta[2]',
             'continent_specific_theta[3]',
             'continent_specific_theta[4]')

plt_5b <- numbers_model_2 %>%
  gather_draws(., !!!syms(thetas)) %>%
  mutate(.variable = case_when(.variable == 'continent_specific_theta[1]' ~ 'africa',
                               .variable == 'continent_specific_theta[2]' ~ 'asia',
                               .variable == 'continent_specific_theta[3]' ~ 'central & northern america',
                               .variable == 'continent_specific_theta[4]' ~ 'europe')) %>%
  ggplot(aes(x= .value, 
             y = .variable,
             slab_colour = .variable, 
             slab_fill = .variable)) +
  stat_slabinterval(point_interval = "median_hdci",
                    slab_alpha = 0.7 ,
                    p_limits = c(0.025, 0.975),
                   # normalize = 'xy',
                    scale = 0.85,
                    .width = 0.95) + 
  
  scale_colour_manual(values = region_colours, aesthetics = 'slab_colour') + 
  scale_fill_manual(values = region_colours, aesthetics = 'slab_fill') + 
  
  scale_x_continuous('Zero-Inflation Probability',
                     breaks = seq(0,1,by=0.25),
                     labels = seq(0,1,by=0.25),
                     expand = c(0,0)) +
  
  scale_y_discrete('Continent',  labels = ~str_to_title(str_wrap(.x, width = 15))) + 
  global_theme + 
  theme(strip.placement  = 'inside',
        strip.text = element_text(face = 'bold', size = 10),
        strip.background = element_blank(),
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 9),
        legend.text = element_text(size = 8))



# C - observed ~ number of sequences
# Extract relevant parameters
inv_logit <- function(x){
  return(exp(x)/(1+exp(x)))
}

continent_specific_detection_samples <- numbers_model_2 %>%
  gather_draws(., continent_specific_detection[i]) 

beta_sequences_samples <-  numbers_model_2 %>%
  gather_draws(., beta_sequences) %>%
  ungroup() %>% 
  select(.draw,
         .chain, 
         .iteration, 
         beta_sequences = .value)

year_detection_samples <- numbers_model_2 %>%
  gather_draws(., year_detection[i]) %>% 
  select(.draw,
         .chain,
         .iteration,
         year_detection = .value)


detection_probabilities <- continent_specific_detection_samples %>%
  rename(continent_specific_detection = .value) %>%
  inner_join(beta_sequences_samples ,
             by = c(".draw", ".chain", ".iteration")) %>%
  inner_join(year_detection_samples, 
             by = c(".draw", ".chain", ".iteration"), 
             relationship = "many-to-many") %>%
  mutate(sequences = list(log1p(1:75))) %>%
  unnest(sequences) %>%
  mutate(
    probability = inv_logit(continent_specific_detection + beta_sequences * sequences + year_detection)
  )


#summary_results <- detection_probabilities %>%
  #group_by(sequences) %>% 
  #median_hdci(probability, .width = 0.95)

plt_5c <- detection_probabilities %>%
  ggplot(aes(x = expm1(sequences), y = probability)) +
  stat_lineribbon(point_interval = "median_hdci", 
                  alpha = 0.7) + 
  
  scale_y_continuous('P(Detection)',
                breaks = seq(0.2, 1, by = 0.2),
                labels = seq(0.2, 1, by = 0.2),
                expand = c(0,0))+ 
  scale_x_continuous('Sequences (n)',
                     expand = c(0,0),
                     breaks = seq(0,75,by=25)) +
  scale_fill_brewer() + 
  global_theme+ 
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
  geom_node_point(aes(color = cluster_class, size=10, fill = cluster_class ), shape = 21, alpha = 0.8) +
  geom_node_label(aes(label = ifelse(cluster_class == 'major', gsub('_.*', '', name), '')), repel = TRUE) +
  scale_size(guide = NULL) + 
  scale_colour_manual(values = class_colours, 'Cluster Class') + 
  scale_fill_manual(values = class_colours, 'Cluster Class') + 
  
  theme_void() + 
  theme(legend.text = element_text(size = 8),
        legend.position = 'inside',
        legend.position.inside = c(0.1,0.9))



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
plt_5f <-avg_predictions(ordinal_model,  variables = list('segments_changed' = 1:7)) %>%
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
  
  scale_x_continuous('Segments changed from previous(N)', expand = c(0,0)) + 
  global_theme + 
  theme(strip.placement  = 'inside',
        strip.text = element_text(face = 'bold', size = 10),
        strip.background = element_blank(),
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 9),
        legend.text = element_text(size = 8))


### G (Time)
plt_5g <-avg_predictions(ordinal_model, variables = list('time_since_last_major' = seq(0, 5, by  = 0.5))) %>%
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

plt_5_lh <- align_plots(plt_5a, plt_5b, plt_5d, plt_5f, align = 'v', axis = 'l')
plt_5_bttm <- plot_grid( plt_5_lh[[2]], plt_5c, plt_5_lh[[3]], plt_5e, plt_5_lh[[4]], plt_5g, ncol = 2, align = 'vh', axis = 'tbrl', labels = c('B', 'C', 'D', 'E', 'F', 'G'), label_size = 9)
plot_grid( plt_5_lh[[1]], plt_5_bttm,  nrow = 2, align = 'v', axis = 'r', labels = c('A', ''), label_size = 9, rel_heights = c(0.25,0.75))
############################################## WRITE ###############################################


############################################## END #################################################
####################################################################################################
####################################################################################################p