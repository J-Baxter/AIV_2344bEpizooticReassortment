sequences_host_month <- meta %>%  
  drop_na(cluster_profile) %>%
  dplyr::select(starts_with('collection_date'),
         host_simplifiedhost,
         collection_regionname) %>%
  mutate(collection_regionname = case_when(grepl('europe', collection_regionname) ~ 'europe',
                                           grepl('africa', collection_regionname) ~ 'africa',
                                           grepl('asia', collection_regionname) ~ 'asia',
                                           grepl('(central|northern) america|caribbean', collection_regionname) ~ 'central & northern america',
                                           grepl('south america|southern ocean', collection_regionname) ~ 'south america',
                                           grepl('australia|melanesia', collection_regionname) ~ 'australasia',
                                           .default = collection_regionname)) %>%
  group_by(collection_datemonth, host_simplifiedhost) %>%
  summarise(n_sequences = n()) %>%
  ungroup() %>%
  mutate(collection_datemonth = ymd(paste0(collection_datemonth, '-01'))) %>%
  drop_na(host_simplifiedhost,collection_datemonth)




sequences_subtype_month <- meta %>%  
  drop_na(cluster_profile) %>%
  dplyr::select(starts_with('collection_date'),
         virus_subtype,
         collection_regionname) %>%
  mutate(collection_regionname = case_when(grepl('europe', collection_regionname) ~ 'europe',
                                           grepl('africa', collection_regionname) ~ 'africa',
                                           grepl('asia', collection_regionname) ~ 'asia',
                                           grepl('(central|northern) america|caribbean', collection_regionname) ~ 'central & northern america',
                                           grepl('south america|southern ocean', collection_regionname) ~ 'south america',
                                           grepl('australia|melanesia', collection_regionname) ~ 'australasia',
                                           .default = collection_regionname)) %>%
  group_by(collection_datemonth, virus_subtype) %>%
  summarise(n_sequences = n()) %>%
  ungroup() %>%
  mutate(collection_datemonth = ymd(paste0(collection_datemonth, '-01'))) %>%
  drop_na(virus_subtype,collection_datemonth)


plt_a <- ggplot(sequences_host_month%>% filter(collection_datemonth > as.Date('2019-12-01'))) + 
  geom_stream(aes(x = collection_datemonth, 
                  y = n_sequences, 
                  colour = host_simplifiedhost,
                  fill = host_simplifiedhost),
              alpha = 0.7) + 
  scale_x_date(limits = as_date(c('2019-01-01', '2024-05-01')), breaks = '1 year', date_labels = "%Y", 'Date') + 
  scale_fill_manual('Host Order', values = host_colours) +
  scale_colour_manual('Host Order', values = host_colours) +
  scale_y_continuous('GISAID Whole Genomes (n)' ,
                     labels = abs,
                     breaks = seq(-225, 225, by = 75),
                     limits = c(-230, 230)) + 
  theme_classic() + 
  theme(axis.text = element_text(size = 8),
        axis.title = element_text(size = 9),
        legend.position = 'inside',
        legend.justification.inside = c(0, 1),
        legend.key.size=  unit(0.4, 'cm'),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 9),
        legend.position.inside = c(0.01,1.05),
        legend.background = element_blank())

plt_b <- ggplot(sequences_subtype_month %>% filter(collection_datemonth > as.Date('2019-12-01'))) + 
  geom_stream(aes(x = collection_datemonth, 
                  y = n_sequences, 
                  colour = virus_subtype,
                  fill = virus_subtype),
              alpha = 0.7) + 
  scale_x_date(limits = as_date(c('2019-01-01', '2024-05-01')), breaks = '1 year', date_labels = "%Y", 'Date') + 
  scale_fill_brewer('Subtype', palette = 'OrRd', direction = -1) +
  scale_colour_brewer('Subtype', palette = 'OrRd', direction = -1) +
  scale_y_continuous('GISAID Whole Genomes (n)' ,
                     labels = abs,
                     breaks = seq(-225, 225, by = 75),
                     limits = c(-230, 230)) + 
  theme_classic() + 
  theme(
    axis.text = element_text(size = 8),
    axis.title = element_text(size = 9),
    legend.position = 'inside',
    legend.key.size=  unit(0.4, 'cm'),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 9),
    legend.justification.inside = c(0, 1),
    legend.position.inside = c(0.01,1.05),
    legend.background = element_blank())



plot_grid(plt_a, plt_b, labels = 'AUTO', label_size = 9, align = 'v', axis = 'lr', nrow = 2)


ggsave('~/Downloads/flu_plots/fig_s2.pdf', height = 130, width = 206, units = 'mm', dpi = 360)
