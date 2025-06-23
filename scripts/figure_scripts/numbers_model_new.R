plt_2a <- test_pred %>%
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

plt_2b <- expand_grid(tmrca_date = seq(ymd('2019-01-01'),ymd('2024-05-01'),by='day'),
                      collection_regionname = c('europe', 'africa', 'asia', 'central & northern america'),
                      cs = NA_real_) %>%
  rows_patch(t, by = c('collection_regionname', 'tmrca_date')) %>%
  group_by(collection_regionname) %>%
  fill(cs) %>%
  ungroup() %>%
  mutate(cs = replace_na(cs, 0)) %>%
  
  ggplot() + 
  geom_step(aes(x=tmrca_date, y=cs, color=collection_regionname), lwd = 0.8) +
  scale_colour_manual('Continent', values = region_colours, labels = str_to_title) +
  scale_y_continuous('Reassortans (Cumulative)', expand = c(0.01, 0)) + 
  scale_x_date(limits = as_date(c('2019-01-01', '2024-02-01')),
               breaks = '1 year', 
               date_labels = "%Y", 'Date (Years)',
               expand = c(0,0)) + 
  theme_classic() + 
  theme(legend.position = 'inside',
        legend.position.inside = c(0,1),
        legend.justification.inside = c(0,1),
        legend.background = element_blank(),
        legend.title = element_text(size = 9),
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 9),
        legend.text = element_text(size = 8))

plt_2c <- h5_diversity_sliding_window %>%
  filter(! continent %in%  c('South America', 'Antarctica')) %>%
  mutate(continent = str_to_lower(continent) %>%
           case_when(grepl('north america', .) ~ 'central & northern america',
                     .default = .)) %>%
  mutate(midpoint = (start_point + end_point)/2) %>%
  mutate(midpoint_date = date_decimal(midpoint)) %>%
  ggplot(aes(y = diversity, x = midpoint_date, colour = continent, fill = continent,label = str_to_title(continent))) +
  #geom_point() + 
  geom_smooth(method = 'gam' , formula = y ~ s(x, bs = 'cs'), se = TRUE, alpha = 0.1)  +
  #geom_labelsmooth(text_smoothing = 30, 
  # fill = "white",
  #formula = y ~ s(x, bs = 'cs'), 
  #method = "gam",
  #size = 4, 
  #linewidth = 0.1,
  #boxlinewidth = 0.3) +
  scale_fill_manual('Continent', values = region_colours, labels = str_to_title) + 
  scale_colour_manual('Continent', values = region_colours, labels = str_to_title) +
  scale_y_continuous('Numbers-Equivalent Shanon Entropy', expand = c(0.01, 0)) + 
  scale_x_datetime(limits = as_datetime(c('2020-06-01', '2024-02-01')),
                   breaks = '1 year', 
                   date_labels = "%Y", 'Date (Years)',
                   expand = c(0,0)) + 
  
  theme_classic()+ 
  theme(legend.position = 'none',
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 9),
        legend.text = element_text(size = 8))



plt_2d <- h5_diversity_similarity_sliding_window %>%
  filter(! continent %in%  c('South America', 'Antarctica')) %>%
  mutate(continent = str_to_lower(continent) %>%
           case_when(grepl('north america', .) ~ 'central & northern america',
                     .default = .)) %>%
  mutate(midpoint = (start_point + end_point)/2) %>%
  mutate(midpoint_date = date_decimal(midpoint)) %>%
  ggplot(aes(y = diversity, x = midpoint_date, colour = continent, fill = continent,label = str_to_title(continent))) +
  #geom_point() + 
  geom_smooth(method = 'gam' , formula = y ~ s(x, bs = 'cs'), se = TRUE, alpha = 0.1)  +
  #geom_labelsmooth(text_smoothing = 30, 
  # fill = "white",
  #formula = y ~ s(x, bs = 'cs'), 
  #method = "gam",
  #size = 4, 
  #linewidth = 0.1,
  #boxlinewidth = 0.3) +
  scale_fill_manual('Continent', values = region_colours, labels = str_to_title) + 
  scale_colour_manual('Continent', values = region_colours, labels = str_to_title) +
  scale_y_continuous('Diversity-Similarity', expand = c(0.01, 0)) + 
  scale_x_datetime(limits = as_datetime(c('2020-0-01', '2024-02-01')),
                   breaks = '1 year', 
                   date_labels = "%Y", 'Date (Years)',
                   expand = c(0,0)) + 
  
  theme_classic()+ 
  theme(legend.position = 'none',
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 9),
        legend.text = element_text(size = 8))

# Genetic distance between reassortants

# Combine
ggdraw() +
  draw_plot(plt_2a, x = 0, y = .5, width = 1, height = .5) +
  draw_plot(plt_2c, x = 0, y = 0, width = .5, height = .5) +
  draw_plot(plt_2d, x = 0.5, y = 0, width = 0.5, height = 0.5) +
  draw_plot_label(label = c("A", "B", "C"), size = 10,
                  x = c(0, 0, 0.5), y = c(1, 0.5, 0.5))


ggsave('~/Downloads/flu_plots/numbers_plot.jpeg', height = 15, width = 25, units = 'cm', dpi = 360)