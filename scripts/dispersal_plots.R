#diffusion coefficient
library(tidyverse)

reassortant_classification <- summary_tbl %>% 
  filter(reassortant_class %in% c('dominant', 'major', 'minor')) %>%
  select(reassortant_class, cluster_profile, colour) %>%
  mutate(reassortant_class = case_when(cluster_profile == '1_1_1_1_1_1_1_1' ~ 'major', .default = reassortant_class)) %>%
  distinct()

dispersal_data <- read_csv('dispersalcoefficients.csv') %>%
  left_join(reassortant_classification) %>% 
  filter(reassortant_class %in% c('dominant', 'major', 'minor'))

dispersal_col <- as.character(dispersal_data$colour)
names(dispersal_col) <- as.character(dispersal_data$cluster_profile)


dispersal_summary = dispersal_data %>% 
  summarise(., median  = median(velocity_coefficent), .by = reassortant_class) 

dispersal <- ggplot(dispersal_data)+
  geom_bar(aes(x = cluster_profile, y = velocity_coefficent, fill = cluster_profile), stat = 'identity') +
  geom_errorbar(aes(x = cluster_profile, ymin = velocity_coefficent_95lower, ymax = velocity_coefficent_95upper)) +
  geom_hline(aes(yintercept = median), data = dispersal_summary, linetype = 'dotted') +
  scale_x_discrete(expand = c(0.15, 0), "Reassortant")+
  scale_y_continuous(expand = c(0,0),
                     #limits = c(0, 12500000),
                     'Dispersal Coefficient')+
  scale_fill_manual(values = dispersal_col) + 
  coord_cartesian(ylim=c(0,10000000))+
  
  facet_grid(cols = vars(reassortant_class), scales = 'free_x', space='free', labeller = as_labeller(c('dominant' = 'Dominant', 'major' = 'Major', 'minor' = 'Minor'))) +
  my_theme +
  theme(axis.text.x = element_text(angle = 45, hjust=1))

ggsave('dispersal.svg', device= 'svg',  height = 110, width = 170, units = 'mm')
Sys.sleep(0.5)
dispersal
dev.off()
