scientific_10 <- function(x) {   parse(text=gsub("e\\+*", " %*% 10^", scales::scientific_format()(x)))}


# Plot Diffusion Coefficient Distribution
plt_5a <- diffusion_data %>%
  ggplot() +
  geom_histogram(aes(x = log1p(weighted_diff_coeff), 
                     fill = weighted_diff_coeff>0,
                     colour = weighted_diff_coeff>0,
                     y = after_stat(density)), 
                 binwidth = 0.5,
                 alpha = 0.7) +
  scale_fill_brewer(palette = 'Accent', 'Is Zero') +
  scale_colour_brewer(palette = 'Accent', 'Is Zero') +
  #scale_x_continuous('Log1p Weighted Diffusion Coefficient' , expand = c(0.01,0.01)) +
  scale_x_continuous(expression(paste('Weighted Diffusion Coefficient (',Km**2~year**-1, ')' )),
                     breaks = log1p(c(0, 10^(seq(from = 1, to = 7, by = 2)))),
                     labels = expression(0,  1%*%10^1,  1%*%10^3, 1%*%10^5,  1%*%10^7),
                     expand = c(0.02,0.02))+

  scale_y_continuous('Probability Density',expand = c(0,0)) + 
  global_theme + 
  theme(legend.position = 'none')



# Plot Persistence distribuions
plt_5b <- diffusion_data %>%
  ggplot() +
  geom_histogram(aes(x = median_anseriformes_wild,
                     y = after_stat(density)), 
                 binwidth = 0.2,
                 fill = host_colours['anseriformes-wild'],
                 colour = host_colours['anseriformes-wild'],
                 alpha = 0.7) +
  scale_x_continuous('Persistence in wild Anseriformes', 
                     expand = c(0.02,0.02), 
                     breaks = seq(from = 0, to = 10, by = 2)) +
  scale_y_continuous('Probability Density',expand = c(0,0)) + 
  coord_cartesian(xlim = c(0, 10))  +
  global_theme + 
  theme(legend.position = 'none')

plt_5c  <- diffusion_data %>%
  ggplot() +
  geom_histogram(aes(x = median_charadriiformes_wild,
                     y = after_stat(density)),
                 binwidth = 0.2,
                 fill = host_colours['charadriiformes-wild'],
                 colour = host_colours['charadriiformes-wild'],
                 alpha = 0.7) +
  scale_x_continuous('Persistence in wild Charadriiformes', 
                     expand = c(0.02,0.02), 
                     breaks = seq(from = 0, to = 6, by = 2)) +
  scale_y_continuous('Probability Density',expand = c(0,0)) + 
  coord_cartesian(xlim = c(0, 6))  +
  global_theme + 
  theme(legend.position = 'none')




####################################################################################################
# Conditional prediction draws  + Conditional marginal means stratified by continent
averages <-  diffusionmodel1_fit %>%
  emmeans(~ collection_regionname,
          var = "weighted_diff_coeff",
          at = list(collection_regionname = unique(diffusion_data$collection_regionname)),
          # epred = TRUE, 
          dpar = "mu",
          re_formula = NA , 
          regrid = "response",
          #tran = "log", 
          type = "response",
          allow_new_levels = TRUE) %>%
  as_tibble() %>%
  mutate(label = round(expm1(emmean)))

plt_5d  <- diffusionmodel1_fit %>%
  predicted_draws(newdata = expand_grid(collection_regionname = unique(as.character(diffusion_data$collection_regionname)),
                                   season = unique(diffusion_data$season)) %>%
                    drop_na() %>%
      
                    mutate(median_anseriformes_wild = median(diffusion_data$median_anseriformes_wild),
                           median_charadriiformes_wild = median(diffusion_data$median_charadriiformes_wild)),
                  re_formula = NA) %>%
  drop_na(collection_regionname) %>%
  ggplot() + 
  geom_histogram(aes(x = log1p(.prediction), y = after_stat(density), colour = collection_regionname, fill = collection_regionname), 
                 alpha = 0.7) + 
  scale_colour_manual(values = region_colours)+
  scale_fill_manual(values = region_colours) + 
  scale_x_continuous(expression(paste('Predicted Weighted Diffusion Coefficient (',Km**2~year**-1, ')' )),
                     breaks = log1p(c(0, 10^(seq(from = 2, to = 10, by = 4)))),
                     labels = expression(0,  1%*%10^2,  1%*%10^6, 1%*%10^10),
                     #limits = c(-0.01, log1p(10^9)),
                     expand = c(0.02,0.02))+
  scale_y_continuous('Probability Density' ,
                     expand = c(0,0))+
  facet_grid(
    cols = vars(collection_regionname),
    labeller =  labeller(collection_regionname=str_to_title)) +
  geom_vline(aes(xintercept = emmean, colour = collection_regionname), data = averages, linetype = 'dashed') +
  geom_text(aes(label =  paste0("E*'('*X*'|'*X*'>'*0*') = '*", label, "~km^2"), 
                colour = collection_regionname),
            parse = T,
            x = 17.5, 
            y = 0.6,
            size = 2.5,
            data = averages) + 
  global_theme + 
  theme(strip.placement  = 'inside')

####################################################################################################



# Posterior Prediction HU ~ Season
plt_5e <- predict_hu_season %>%
  gather_emmeans_draws() %>%
  ggplot(aes(x  = 1-.value,
             y = season,
             slab_colour = season,
             slab_fill = season)) +
  stat_halfeye(slab_alpha = 0.7,
               p_limits = c(0.001, 0.999),
               point_interval = "median_hdi",
               linewidth = 1.5,
               .width =  0.95) +
  scale_colour_manual(values = c('#E8E1E9FF', '#C0A5AAFF', '#4D3944FF', '#7083A4FF'), aesthetics = 'slab_colour') +
  scale_fill_manual(values = c('#E8E1E9FF', '#C0A5AAFF', '#4D3944FF', '#7083A4FF'), aesthetics = 'slab_fill') +
  scale_x_continuous('P(Weighted Diffusion Coefficient > 0)') + 
  scale_y_discrete('Season',
                   labels = c('overwintering' = 'Overwintering',
                              'migrating_spring' = 'Spring Migration',
                              'migrating_autumn' = 'Autumn Migration',
                              'breeding' = 'Breeding')) + 
  global_theme 


# Posterior Prediction/Average HU ~ Region
plt_5f <-regional_average %>%
  gather_emmeans_draws() %>%
  ggplot(aes(x  = 1-.value,
             y = collection_regionname, 
             slab_colour = collection_regionname,
             slab_fill = collection_regionname)) +
  stat_halfeye(slab_alpha = 0.7,
               p_limits = c(0.001, 0.999),
               point_interval = "median_hdi",
               linewidth = 1.5,
               .width =  0.95) +
  scale_fill_manual(values = region_colours, aesthetics = 'slab_fill') +
  scale_colour_manual(values = region_colours, aesthetics = 'slab_colour') + 
  scale_x_continuous('P(Weighted Diffusion Coefficient > 0)') + 
  scale_y_discrete('Continent', 
                   labels = function(x) str_to_title(x) %>% str_wrap(., width = 10)) + 
  global_theme + 
  theme(legend.position = 'none')


# average marginal effect of anseriformes/ charadriiformes
plt_5g <- bind_rows(diffusionmodel1_fit %>%
  emtrends(~ median_anseriformes_wild,
           var = "median_anseriformes_wild",
           at = list(median_anseriformes_wild = c(0.25 , 0.5 , 1)),
            epred = TRUE, 
           dpar = "mu") %>% 
  gather_emmeans_draws() %>%
  mutate(var = 'median_anseriformes_wild') %>%
  rename(var_change = median_anseriformes_wild),
  
  diffusionmodel1_fit %>%
    emtrends(~ median_charadriiformes_wild,
             var = "median_charadriiformes_wild",
             at = list(median_charadriiformes_wild = c(0.25 , 0.5 , 1)),
             epred = TRUE, 
             dpar = "mu") %>% 
    gather_emmeans_draws() %>%
    mutate(var = 'median_charadriiformes_wild') %>%
    rename(var_change = median_charadriiformes_wild)) %>%
  
  mutate(var = gsub('median_', '', var) %>%
           gsub('_', '-', .)) %>%

  ggplot(aes(x  = .value,
             y = as.factor(var_change),
             slab_colour = var,
             slab_fill = var)) +
  stat_halfeye(slab_alpha = 0.7 ,
               p_limits = c(0.001, 0.999),
               point_interval = "median_hdi",
               linewidth = 1.5,
               .width =  0.95)  +
  geom_vline(aes(xintercept = 0), linetype = 'dashed')  +
  scale_x_continuous(expression(paste('Change in Weighted Diffusion Coefficient (',Km**2~year**-1, ')' )),
                     breaks = seq(from = -1*10**6, to = 5*10**6, by = 2*10**6),
                     labels = scientific_10,
                     limits = c(-1*10**6, 5*10**6),
                     expand = c(0.02,0.02))+
  scale_fill_manual(values = host_colours, aesthetics = 'slab_fill') +
  scale_colour_manual(values = host_colours, aesthetics = 'slab_colour') + 
  scale_y_discrete('Persistence in Host (Years)', 
                   labels = function(x) str_to_title(x) %>% str_wrap(., width = 10)) +
  global_theme + 
  theme(legend.position = 'none')


  
  
# Marginal effect of continent on MU
# Results are averaged over the levels of: season


plt_5h <-diffusionmodel1_fit %>%
  emmeans(~ collection_regionname,
          dpar = "mu",
          epred = TRUE) %>% 
  contrast(method = "revpairwise") %>% 
  gather_emmeans_draws() %>%
  filter(grepl('asia', contrast)) %>%
  mutate(contrast = gsub(' - .*', '' , contrast)) %>%
  ggplot(aes(x  = .value,
             y = contrast,
             slab_colour = as.factor(contrast),
             slab_fill = as.factor(as.factor(contrast)))) +
  geom_vline(aes(xintercept = 0), linetype = 'dashed')  +
  stat_halfeye(slab_alpha = 0.7,
               p_limits = c(0.001, 0.999),
               point_interval = "median_hdi",
               linewidth = 1.5,
               .width =  0.95)+
  scale_fill_manual(values = region_colours, aesthetics = 'slab_fill') +
  scale_colour_manual(values = region_colours, aesthetics = 'slab_colour') + 
  scale_x_continuous(expression(paste('Change in Weighted Diffusion Coefficient (',Km**2~year**-1, ')' )),
                     breaks = seq(from = -2*10**6, to = 5*10**6, by = 1*10**6),
                     labels = scientific_10,
                     expand = c(0.02,0.02))+
  scale_y_discrete('Continent', 
                   labels = function(x) str_to_title(x) %>% str_wrap(., width = 10)) + 
  global_theme + 
  theme(legend.position = 'none')
             
  
  
  
# Combine plots together
plt5_lh <- align_plots(plt_5a, plt_5d, plt_5e, plt_5g,align = 'v', axis = 'l')


plt5_top <- plot_grid(plt5_lh[[1]], plt_5b, plt_5c, align = 'h', nrow = 1, axis = 'tb', labels = 'AUTO')
plt5_bottom <- plot_grid(plt5_lh[[3]], plt_5f, plt5_lh[[4]], plt_5h, align = 'hv', 
                         ncol = 2, nrow = 2, axis = 'tblr', labels = c('E', 'F', 'G', 'H'))


plt5 <- plot_grid(plt5_top, plt5_lh[[2]], plt5_bottom, labels = c('', 'D', ''), nrow = 3, rel_heights = c(1,1,2))
plt5

ggsave('~/Downloads/figure6.jpeg', height = 40, width = 35, units = 'cm', dpi = 360)
