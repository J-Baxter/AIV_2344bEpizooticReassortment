global_theme <- theme_classic()+
  theme(
    #text = element_text(size=10),
    #axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0), size = 10),
    #axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0), size = 10),
    #axis.text = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0), size = 8),
    legend.text = element_text(size = 8),
    legend.position = 'none', 
    panel.spacing = unit(2, "lines"), 
    strip.background = element_blank()
  )


plt_5a  <- diffusionmodel1_fit_gamma_19 %>%
  #  back-transformed linear predictive draws, equivalent to add_linpred (and in this case equviv to epred)
  # this uses the empirical data distribution by default
  avg_predictions(by = 'collection_regionname',  type = 'response') %>% 
  get_draws() %>%
  as_tibble() %>%
  ggplot() + 
  geom_histogram(aes(x = log1p(draw), y = after_stat(density), colour = collection_regionname, fill = collection_regionname), 
                 binwidth = 0.2,
                 alpha = 0.7) + 
  scale_colour_manual(values = region_colours)+
  scale_fill_manual(values = region_colours) + 
  scale_x_continuous(expression(paste('Predicted Weighted Diffusion Coefficient (',Km**2~year**-1, ')' )),
                     breaks =log(10^seq(5, 8, b = 1)),
                     labels = expression(1%*%10^5, 1%*%10^6,  1%*%10^7, 1%*%10^8),
                     limits = c(11.5, 19),
                     expand = c(0,0)
  )+
 
  facet_grid(
    cols = vars(collection_regionname),
    labeller =  labeller(collection_regionname=str_to_title),
    scales = 'free_y') +
  scale_y_continuous('Probability Density' ,
                     breaks = seq(0,1.5,by=0.5),
                     labels = seq(0,1.5,by=0.5),
                     expand = c(0,0)) +
  #geom_vline(aes(xintercept = emmean, colour = collection_regionname), data = averages, linetype = 'dashed') +
  #geom_text(aes(label =  paste0("E*'('*X*'|'*X*'>'*0*') = '*", label, "~km^2"), 
  #      colour = collection_regionname),
  #   parse = T,
  # x = 17.5, 
  #  y = 0.6,
  # size = 2.5,
  # data = averages) + 
  global_theme + 
  theme(strip.placement  = 'inside',
        strip.text = element_text(face = 'bold', size = 10),
        strip.background = element_blank(),
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 9),
        legend.text = element_text(size = 8))


# count_cross_species_log1p (single plot)
plt_5b <- avg_predictions(diffusionmodel1_fit_gamma_19, 
                variables = list('count_cross_species_log1p' = log1p(0:10)))%>%
  get_draws(shape = 'rvar') %>%
  ggplot(aes(x = expm1(count_cross_species_log1p), ydist = rvar)) +
  stat_lineribbon(point_interval = "median_hdci",
                  alpha = 0.7 ,
                  #p_limits = c(0.025, 0.975),
                  #normalize = 'xy'
                  ) +
  scale_y_log10(expression(paste('Predicted Diffusion Coefficient (',Km**2~year**-1, ')' )),
                breaks = 10^seq(5.5, 6.5, by =0.5),
                labels = expression(5%*%10^5, 1%*%10^6, 1%*%10^6.5),
                expand = c(0,0))+ 
  scale_x_continuous('Host Jumps',
                     expand = c(0,0),
                     breaks = seq(0,10,by=2)) +
  scale_fill_brewer() + 
  global_theme+ 
  theme(strip.placement  = 'inside',
        strip.text = element_text(face = 'bold', size = 10),
        strip.background = element_blank(),
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 9),
        legend.text = element_text(size = 8))



# collection_year (single plot)
plt_5c <- avg_predictions(diffusionmodel1_fit_gamma_19, by = 'collection_year') %>%
get_draws(shape = 'rvar') %>%
ggplot(aes(x = collection_year, ydist = rvar, slab_colour = NULL,
           slab_fill = estimate)) +
  stat_slabinterval(point_interval = "median_hdci",
                    slab_alpha = 0.7 ,
                    p_limits = c(0.025, 0.975),
                    normalize = 'xy',
                    scale = 0.85,
                    side = 'both',
                    .width = 0.95) +
  scale_y_log10(expression(paste('Predicted Diffusion Coefficient (',Km**2~year**-1, ')' )),
                breaks = 10^seq(4, 8, b = 1),
                labels = expression(1%*%10^4, 1%*%10^5, 1%*%10^6,  1%*%10^7, 1%*%10^8),
                expand = c(0,0))+ 
  scale_x_discrete('Year of Reassortant MRCA') +
  scale_fill_distiller(palette = 'Greens', aesthetics = 'slab_fill', direction = 1) + 
  scale_colour_distiller(palette = 'Greens', aesthetics = 'slab_colour', direction = 1) + 
  global_theme+ 
  theme(strip.placement  = 'inside',
        strip.text = element_text(face = 'bold', size = 10),
        strip.background = element_blank(),
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 9),
        legend.text = element_text(size = 8))


# median_anseriformes_wild_prop (strat by region)
plt_5d <- avg_predictions(diffusionmodel1_fit_gamma_19, 
                variables = list('path_prop_anseriformes_wild' = seq(0,1,by = 0.1)),
                by = c('collection_regionname', 'path_prop_anseriformes_wild')) %>%
  get_draws(shape = 'rvar') %>%
  ggplot(aes(x = path_prop_anseriformes_wild, ydist = rvar, fill = collection_regionname)) +
  stat_lineribbon(point_interval = "median_hdci", 
                  alpha = 0.5) +
 # geom_point(data = diffusion_data, aes(x = path_prop_anseriformes_wild, y = weighted_diff_coeff), inherit.aes = F) +
  
  facet_wrap(~collection_regionname, ncol = 4,
             labeller =  labeller(collection_regionname=str_to_title)) +
  scale_y_log10(expression(paste('Predicted Diffusion Coefficient (',Km**2~year**-1, ')' )),
                breaks = 10^seq(5, 8, by =1),
                labels = expression(1%*%10^5, 1%*%10^6, 1%*%10^7, 1%*%10^8),
                expand = c(0,0))+ 
  scale_x_continuous('Proportion of Circulation in Anseriformes',
                     expand = c(0,0),
                     breaks = seq(0,1,by=0.25)) +
  scale_fill_manual(values = region_colours)+
  global_theme+ 
  theme(strip.placement  = 'inside',
        strip.text = element_text(face = 'bold', size = 10),
        strip.background = element_blank(),
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 9),
        legend.text = element_text(size = 8))
  

# median_charadriiformes_wild_prop (strat by region)
plt_5e <- avg_predictions(diffusionmodel1_fit_gamma_19, 
                variables = list('path_prop_charadriiformes_wild' = seq(0,1,by = 0.1)),
                by = c('collection_regionname','path_prop_charadriiformes_wild')) %>%
  get_draws(shape = 'rvar') %>%
  ggplot(aes(x = path_prop_charadriiformes_wild, ydist = rvar, fill = collection_regionname)) +
  stat_lineribbon(point_interval = "median_hdci",
                  alpha = 0.5) +
 # geom_point(data = diffusion_data, aes(x = path_prop_charadriiformes_wild, y = weighted_diff_coeff), inherit.aes = F) +
  facet_wrap(~collection_regionname, ncol = 4,
             labeller =  labeller(collection_regionname=str_to_title)) +
  scale_y_log10(expression(paste('Predicted Diffusion Coefficient (',Km**2~year**-1, ')' )),
                breaks = 10^seq(5, 8, by =1),
                labels = expression(1%*%10^5, 1%*%10^6, 1%*%10^7, 1%*%10^8),
                expand = c(0,0))+ 
  scale_x_continuous('Proportion of Circulation in Charadriiformes',
                     expand = c(0,0),
                     breaks = seq(0,1,by=0.25)) +
  scale_fill_manual(values = region_colours)+
  coord_cartesian(ylim =c(10^4.5, 10^8.5)) +
  global_theme+ 
  theme(strip.placement  = 'inside',
        strip.text = element_text(face = 'bold', size = 10),
        strip.background = element_blank(),
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 9),
        legend.text = element_text(size = 8))

plt_5f <-avg_predictions(diffusionmodel1_fit_gamma_19, 
                variables = list('path_prop_galliformes_domestic' = seq(0,1,by = 0.1)),
                by = c('collection_regionname', 'path_prop_galliformes_domestic')) %>%
  get_draws(shape = 'rvar') %>%
  ggplot(aes(x = path_prop_galliformes_domestic, ydist = rvar, fill = collection_regionname)) +
  stat_lineribbon(point_interval = "median_hdci",
                  alpha = 0.5) +
  
 # geom_point(data = diffusion_data, aes(x = path_prop_galliformes_domestic, y = weighted_diff_coeff), inherit.aes = F) +
  facet_wrap(~collection_regionname, ncol = 4,
             labeller =  labeller(collection_regionname=str_to_title)) +
  scale_y_log10(expression(paste('Predicted Diffusion Coefficient (',Km**2~year**-1, ')' )),
                breaks = 10^seq(5, 8, by =1),
                labels = expression(1%*%10^5, 1%*%10^6, 1%*%10^7, 1%*%10^8),
                expand = c(0,0))+ 
  scale_x_continuous('Proportion of Circulation in Galliformes',
                     expand = c(0,0),
                     breaks = seq(0,1,by=0.25)) +
  scale_fill_manual(values = region_colours)+
  coord_cartesian(ylim =c(10^4.5, 10^8.5)) +
  global_theme+ 
  theme(strip.placement  = 'inside',
        strip.text = element_text(face = 'bold', size = 10),
        strip.background = element_blank(),
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 9),
        legend.text = element_text(size = 8))

plt5_lh <- align_plots(plt_5a, plt_5b, plt_5d, plt_5e, align = 'v', axis = 'l')
plt_5mid <- plot_grid(plt5_lh[[2]], plt_5c, ncol = 2, align='h', labels = c('B', 'C'), label_size = 9)


plot_grid(plt_5d, plt_5e, plt_5f, nrow = 3, labels = 'AUTO' ,label_size = 10)



ggdraw() +
  draw_plot(plt_5a, x = 0, y = .5, width = 1, height = .5) +
  draw_plot(plt_5b, x = 0, y = 0, width = .5, height = .5) +
  draw_plot(plt_5c, x = 0.5, y = 0, width = 0.5, height = 0.5) +
  draw_plot_label(label = c("A", "B", "C"), size = 10,
                  x = c(0, 0, 0.5), y = c(1, 0.5, 0.5))



ggsave('~/Downloads/flu_plots/figure6.jpeg', height = 20, width = 25, units = 'cm', dpi = 360)



plot_grid(plt_5d, plt_5e, plt_5f, nrow = 3, labels = 'AUTO' ,label_size = 10)
ggsave('~/Downloads/flu_plots/diffusion_hosteffects.jpeg', height = 25, width = 25, units = 'cm', dpi = 360)

