# Base plots
base_plot <- count_data %>%
  drop_na(collection_regionname) %>%
  as_tibble() %>%  
  filter(collection_datemonth > as_date('2018-06-01'))  %>%
  mutate(across(where(is.numeric), .fns = ~replace_na(.x, 0))) %>%
  pivot_longer(starts_with('woah'), names_to = 'woah_metric', values_to = 'woah_values') %>%
  
  ggplot() + 
  
  scale_fill_manual(values = region_colours) +
  scale_colour_manual(values = region_colours) +
  theme_classic() + 
  scale_x_date(limits = as_date(c('2019-01-01', '2024-05-01')),
               breaks = '2 year', 
               date_labels = "%Y", 'Date',
               expand = c(0,0)) + 
  
  facet_wrap(~collection_regionname,  
             scales = 'free_y',
             ncol = 2) + 
  
  geom_text(data = count_data %>% select(collection_regionname) %>% drop_na() %>% distinct(), 
            aes(label = str_wrap(str_to_title(collection_regionname), width = 20)),
            y = Inf, 
            x = as_date('2019-01-01'), 
            size = 3,
            fontface = 'bold',
            vjust = 'top',
            hjust = 'left') +
  theme(
    strip.background = element_blank(),
    strip.text = element_blank(),
    axis.text = element_text(size = 8),
    axis.title = element_text(size = 10))


base_plot + 
  geom_bar(data = sequences_month, 
           aes(x = collection_datemonth, 
               y = n_sequences*3000, 
               colour = NULL,
               fill = collection_regionname), 
           alpha = 0.3, 
           stat = 'identity',
           position = 'stack') +
  
  geom_line(aes(x = collection_datemonth, 
                y = woah_values/6, 
                colour= collection_regionname,
                linetype = woah_metric)) +
 
  scale_y_continuous(expand = c(0.01, 0),
                     'H5 HPAIV',
                     sec.axis = sec_axis(transform = ~ ./3000, 
                                         'GISAID Whole Genomes (n)')) +

  
  guides(fill = 'none',
         colour = 'none',
         linetype=guide_legend(keywidth = 3, keyheight = 1, title = NULL)) + 
  theme(legend.position = c(0.8, 0.1))


ggsave('~/Downloads/flu_plots/figureS1.jpeg', height = 17, width = 17, units = 'cm', dpi = 360)

  