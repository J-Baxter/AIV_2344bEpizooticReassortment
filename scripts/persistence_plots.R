# evolutionary rate x persistence time

library(tidyverse)

data <- read_csv('asia_pb2_ori_mcc_Rcode_stats_joint.csv') %>%
  mutate(plot_HostType = case_when(grepl(',', Range_of_HostType) ~ 'mix',
                                   .default = Range_of_HostType))

persitence <- ggplot(data %>% filter(plot_HostType %in% c('mix', 'domestic-ave', 'wild-ave')) ,aes(x = persist.time, y = evoRate, colour = plot_HostType))+
  geom_point() + my_theme + 
  #draws a 95% confidence level for a multivariate t-distribution
  stat_ellipse(aes(color = plot_HostType), type = "t")+
  scale_y_continuous('Evolutionary Rate',
                     expand = c(0.001,0),
                     limits = c(0, 0.015)
  )+
  scale_x_continuous('Persistence Time', expand = c(0,0), limits = c(0,2.5)) +
  scale_colour_brewer(palette = 'OrRd', 'Host Type', labels = c('wild-ave' = 'Wild bird',
                                                                'domestic-ave' = 'Domestic bird', 
                                                                'mix' = 'Mixed')) +
  theme(legend.position = 'right') + 
  


ggsave('persitence_domwildsplitpb2.svg', device= 'svg',  height = 110, width = 140, units = 'mm')
Sys.sleep(0.5)
data %>%
  separate_wider_delim(.,Host_of_Earliest_Sample_Date,
                       delim = '-',
                       names = c('wild', 'ancestral_host')) %>% 
  filter(wild!='NA') %>%
  mutate(fade = case_when(persist.time <2 & evoRate < 0.01 ~ '1',
                          .default = '0')) %>%
  mutate(ancestral_host = case_when(ancestral_host %in% c('Passeriformes', 'Accipitriformes', 'Charadriiformes', 'unknown avian') ~ 'Other', .default = ancestral_host)) %>%
  ggplot(. ,aes(x = persist.time, y = evoRate, colour = ancestral_host, alpha = fade)) +
  geom_point() + my_theme + 
  #draws a 95% confidence level for a multivariate t-distribution
  stat_ellipse(aes(color = ancestral_host), type = "t")+
  scale_alpha_discrete(range = c(0.35, 0.9)) +
  scale_y_continuous('Evolutionary Rate',
                     expand = c(0.001,0),
                     limits = c(0, 0.015)
  )+
  scale_x_continuous('Persistence Time', expand = c(0,0), limits = c(0,2.5)) +
  scale_colour_discrete('Ancestral Host Type') +
 # scale_colour_brewer('Ancestral Host Type') +
  theme(legend.position = 'right')+
  facet_wrap(.~wild)  + 
  geom_smooth(method='lm', aes(group = wild,x = persist.time, y = evoRate), colour = 'black', data = data %>%
                separate_wider_delim(.,Host_of_Earliest_Sample_Date,
                                     delim = '-',
                                     names = c('wild', 'ancestral_host')) %>% 
                filter(wild!='NA') %>%
                filter(persist.time <2 & evoRate < 0.01), inherit.aes = FALSE )
dev.off()

ggsave('persitence.svg', device= 'svg',  height = 110, width = 140, units = 'mm')
Sys.sleep(0.5)
persitence
dev.off()


summary_tbl <- read_csv('summary_reassortant_metadata.csv') %>%
  rename(reassortant_class = group) %>%
  filter(reassortant_class %in% c('dominant', 'major', 'minor')) %>%
  select(reassortant_class, cluster_profile, colour) %>%
  mutate(reassortant_class = case_when(cluster_profile == '1_1_1_1_1_1_1_1' ~ 'major', .default = reassortant_class)) %>%
  distinct()
dispersal_col <- as.character(china_reassortants$colour)
names(dispersal_col) <- as.character(china_reassortants$cluster.profile)

china_reassortants %>% 
  #left_join(summary_tbl, by = join_by(cluster.profile == cluster_profile)) %>%
  ggplot() +
  geom_count(aes(x = cluster.profile, y = virus.subtype, colour = cluster.profile)) +
  scale_x_discrete( "Reassortant") + 
  scale_colour_manual(values = dispersal_col) +
  facet_grid(cols = vars(reassortant_class), scales = 'free_x', space='free', labeller = as_labeller(c('dominant' = 'Dominant', 'major' = 'Major', 'minor' = 'Minor'))) +
  my_theme +
  theme(axis.text.x = element_text(angle = 45, hjust=1))
  

