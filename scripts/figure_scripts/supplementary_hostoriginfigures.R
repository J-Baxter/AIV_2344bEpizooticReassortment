####################################################################################################
####################################################################################################
## Script name:
##
## Purpose of script:
##
## Date created: 2025-01-15
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


############################################## DATA ################################################



############################################## MAIN ################################################
plt_s1a <- combined_data %>% 
  select(collection_regionname, segment) %>%
  mutate(segment = factor(segment, levels = c('pb2', 'pb1', 'pa', 'ha', 'np', 'nx', 'mp', 'ns'))) %>%
  mutate(collection_regionname = case_when(grepl('europe', collection_regionname) ~ 'europe',
                                           grepl('africa', collection_regionname) ~ 'africa',
                                           grepl('asia', collection_regionname) ~ 'asia',
                                           grepl('(central|northern) america', collection_regionname) ~ 'central & northern america',
                                           grepl('south america|southern ocean', collection_regionname) ~ 'south america',
                                           grepl('australia', collection_regionname) ~ 'australasia',
                                           .default = collection_regionname)) %>%
  drop_na(collection_regionname) %>%
  summarise(n = n(), .by = c(collection_regionname, segment)) %>%
  ggplot() + 
  geom_bar(aes(x = collection_regionname, y = n, fill = collection_regionname), stat = 'identity') +
  scale_fill_manual(values = region_colours) + 
  scale_x_discrete('Origin Region', expand = c(0.15,0.1),labels = c('europe' = 'EUR',
                                                                      'africa' = 'AFR',
                                                                      'asia' = 'ASIA',
                                                                      'central & northern america' = 'AMR')) + 
  scale_y_continuous('Frequency', expand = c(0,0)) + 
  facet_wrap(.~segment,
             labeller = labeller(segment = toupper),
             nrow = 2) +
  global_theme+ 
  theme(legend.position = 'none')




test_sankey <- combined_data %>% 
  select(collection_regionname, segment, cluster_profile) %>%
  mutate(segment = factor(segment, levels = c('pb2', 'pb1', 'pa', 'ha', 'np', 'nx', 'mp', 'ns'))) %>%
  mutate(collection_regionname = case_when(grepl('europe', collection_regionname) ~ 'europe',
                                           grepl('africa', collection_regionname) ~ 'africa',
                                           grepl('asia', collection_regionname) ~ 'asia',
                                           grepl('(central|northern) america', collection_regionname) ~ 'central & northern america',
                                           grepl('south america|southern ocean', collection_regionname) ~ 'south america',
                                           grepl('australia', collection_regionname) ~ 'australasia',
                                           .default = collection_regionname)) %>%
  drop_na(collection_regionname) %>%

  pivot_wider(names_from = segment, values_from = collection_regionname, values_fn = ~.x[1]) %>%
  drop_na() %>%
  make_long(pb2, pb1, pa, ha, np, nx, mp, ns)


plt_s1b <-ggplot(test_sankey, aes(x = x, 
               next_x = next_x, 
               node = node, 
               next_node = next_node,
               fill = factor(node))) +
  geom_sankey(node.color = 1, flow.alpha = 0.5) +
  scale_x_discrete('Segment', labels = toupper) + 
  scale_fill_manual(values = region_colours, 'Origin Region') + 
  theme_sankey() + 
  theme(legend.position = 'bottom')

plt_s1 <- plot_grid(plt_s1a, plt_s1b, labels = 'AUTO', align = 'hv', axis = 'lftb', nrow = 2)
plt_s1

ggsave('~/Downloads/supfig_originbyregion.jpeg',
       height = 25,
       width = 20,
       units = 'cm',
       dpi = 360)




plt_s2a <- combined_data %>% 
  mutate(host_simplifiedhost = gsub('\\+.*', '', host_simplifiedhost)) %>%
  select(host_simplifiedhost, segment) %>%
  mutate(segment = factor(segment, levels = c('pb2', 'pb1', 'pa', 'ha', 'np', 'nx', 'mp', 'ns'))) %>%
 
  drop_na(host_simplifiedhost) %>%
  summarise(n = n(), .by = c(host_simplifiedhost, segment)) %>%
  ggplot() + 
  geom_bar(aes(x = host_simplifiedhost, y = n, fill = host_simplifiedhost), stat = 'identity') +
  scale_fill_manual(values = host_colours) + 
  scale_x_discrete('Origin Host', expand = c(0.1,0.1) ,labels = c('anseriformes-domestic' = 'ANS-dom',
                                                                    'anseriformes-wild' = 'ANS-wild',
                                                                    'charadriiformes-wild' = 'CHAR-wild',
                                                                    'galliformes-domestic' = 'GAL-dom',
                                                                    'galliformes-wild' = 'GAL-wild',
                                                                    'environment' = 'ENV',
                                                                    'mammal' = 'MAM',
                                                                    'other-bird' = 'OTHER')) + 
  scale_y_continuous('Frequency', expand = c(0,0)) + 
  facet_wrap(.~segment,
             labeller = labeller(segment = toupper),
             nrow = 2) +
  global_theme+ 
  theme(legend.position = 'none',
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

temp <- combined_data %>% 
  select(host_simplifiedhost, segment, cluster_profile) %>%
  mutate(segment = factor(segment, levels = c('pb2', 'pb1', 'pa', 'ha', 'np', 'nx', 'mp', 'ns'))) %>%
  group_by(cluster_profile) %>%
  mutate(check = n() == 8) %>%
  filter(check) %>%
  ungroup() %>%
  select(-check) %>%
  pivot_wider(names_from = segment, values_from = host_simplifiedhost, values_fn = ~.x[1])   %>%
  drop_na() %>%
  rowwise() %>%
  mutate(line = case_when(c_across(2:9) %>% n_distinct(na.rm = TRUE) == 1 ~ NA, .default = '1')) %>%
  select(cluster_profile, line)



test_sankey_2 <- combined_data %>% 
  select(host_simplifiedhost, segment, cluster_profile) %>%
  mutate(segment = factor(segment, levels = c('pb2', 'pb1', 'pa', 'ha', 'np', 'nx', 'mp', 'ns'))) %>%
  drop_na(host_simplifiedhost) %>%
  
  pivot_wider(names_from = segment, values_from = host_simplifiedhost, values_fn = ~.x[1]) %>%
  drop_na() %>%
  make_long(pb2, pb1, pa, ha, np, nx, mp, ns)


plt_s2b <-ggplot(test_sankey_2, aes(x = x, 
                                  next_x = next_x, 
                                  node = node, 
                                  next_node = next_node,
                                  fill = factor(node))) +
  geom_sankey(node.color = 1, flow.alpha = 0.5) +
  scale_x_discrete('Segment', labels = toupper) + 
  scale_fill_manual(values = host_colours, 'Origin Host') + 
  theme_sankey() + 
  theme(legend.position = 'bottom')

plt_s2<- plot_grid(plt_s2a, plt_s2b, labels = 'AUTO', align = 'hv', axis = 'lftb', nrow = 2)
plt_s2

ggsave('~/Downloads/supfig_originbyhost.jpeg',
       height = 25,
       width = 20,
       units = 'cm',
       dpi = 360)



combined_data %>% 
  select(prob_host_simplifiedhost,host_simplifiedhost, cluster_profile, segment) %>%
  mutate(segment = factor(segment, levels = c('pb2', 'pb1', 'pa', 'ha', 'np', 'nx', 'mp', 'ns'))) %>%
  inner_join(temp) %>%
  ggplot(aes(x = segment, y = prob_host_simplifiedhost, colour = host_simplifiedhost)) + 
  geom_point() + 
  geom_line(aes(alpha = line, group = cluster_profile), colour = 'black') + 
  #scale_colour_manual(values = host_colours)+
  scale_alpha_manual(values = 0.1, na.value = 0) +
  global_theme
############################################## WRITE ###############################################




############################################## END #################################################
####################################################################################################
####################################################################################################
