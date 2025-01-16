# Plot Diffusion Coefficient Distribution
plt_diffusiondata <- diffusion_data %>%
  ggplot() +
  geom_histogram(aes(x = log1p(weighted_diff_coeff), 
                     fill = weighted_diff_coeff>0,
                     colour = weighted_diff_coeff>0,
                     y = after_stat(density)), 
                 binwidth = 0.5,
                 alpha = 0.7) +
  scale_fill_brewer(palette = 'Set1', 'Is Zero') +
  scale_colour_brewer(palette = 'Set1', 'Is Zero') +
  scale_x_continuous('Log1p Weighted Diffusion Coefficient' , expand = c(0.01,0.01)) +
  scale_y_continuous('Probability Density',expand = c(0,0)) + 
  global_theme + 
  theme(legend.position = 'none')

# Plot Persistence distribuions
plt_anseriformespersistence <- diffusion_data %>%
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

plt_charadriiformespersistence <- diffusion_data %>%
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

plt_diffusionmodel_d <- diffusionmodel1_fit %>%
  predicted_draws(newdata = expand_grid(collection_regionname = unique(diffusion_data$collection_regionname),
                                   season = unique(diffusion_data$season)) %>%
                    drop_na() %>%
                    mutate(median_anseriformes_wild = median(diffusion_data$median_anseriformes_wild),
                           median_charadriiformes_wild = median(diffusion_data$median_charadriiformes_wild)),
                  re_formula = NA) %>%
  ggplot() + 
  geom_histogram(aes(x = log1p(.prediction), y = after_stat(density), colour = collection_regionname, fill = collection_regionname), 
                 alpha = 0.7) + 
  scale_colour_manual(values = region_colours)+
  scale_fill_manual(values = region_colours) + 
  scale_x_continuous(expression(paste('Predicted Branch-Weighted Diffusion Coefficient (',Km**2~year**-1, ')' )),
                     breaks = log1p(c(0, 10^(seq(from = 1, to = 10)))),
                     labels = expression(0, 10^1, 10^2, 10^3, 10^4, 10^5,  10^6, 10^7, 10^8, 10^9, 10^10),
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
plt_hu_season <- predict_hu_season %>%
  gather_emmeans_draws() %>%
  ggplot(aes(x  = 1-.value,
             y = season,
             slab_colour = season,
             slab_fill = season)) +
  stat_halfeye(slab_alpha = 0.7) +
  scale_colour_brewer(palette = 'GnBu', aesthetics = 'slab_colour') +
  scale_fill_brewer(palette = 'GnBu', aesthetics = 'slab_fill') +
  scale_x_continuous('P(Weighted Diffusion Coefficient > 0)') + 
  scale_y_discrete('Season',
                   labels = c('overwintering' = 'Overwintering',
                              'migrating_spring' = 'Spring Migration',
                              'migrating_autumn' = 'Autumn Migration',
                              'breeding' = 'Breeding')) + 
  global_theme 


# Posterior Prediction/Average HU ~ Region
regional_average %>%
  gather_emmeans_draws() %>%
  ggplot(aes(x  = 1-.value,
             y = collection_regionname, 
             slab_colour = collection_regionname,
             slab_fill = collection_regionname)) +
  stat_halfeye(slab_alpha = 0.7) +
  scale_fill_manual(values = region_colours, aesthetics = 'slab_fill') +
  scale_colour_manual(values = region_colours, aesthetics = 'slab_colour') + 
  scale_x_continuous('P(Weighted Diffusion Coefficient > 0)') + 
  scale_y_discrete('Continent', 
                   labels = function(x) str_to_title(x) %>% str_wrap(., width = 10)) + 
  global_theme + 
  theme(legend.position = 'none')


# average marginal effect of anseriformes/ charadriiformes
diffusionmodel1_fit %>%
  emtrends(~ median_anseriformes_wild,
           var = "median_anseriformes_wild",
           at = list(median_anseriformes_wild = c(0.25 , 0.5 , 1)),
            epred = TRUE, 
           dpar = "mu") %>% 
  gather_emmeans_draws() %>%
  ggplot(aes(x  = .value,
             slab_colour = as.factor(median_anseriformes_wild),
             slab_fill = as.factor(as.factor(median_anseriformes_wild)))) +
  stat_halfeye(slab_alpha = 0.7)  
  

diffusionmodel1_fit %>%
  emtrends(~ median_charadriiformes_wild,
           var = "median_charadriiformes_wild",
           at = list(median_charadriiformes_wild = c(0.25 , 0.5 , 1)),
           epred = TRUE, 
           dpar = "mu") %>% 
  gather_emmeans_draws() %>%
  ggplot(aes(x  = .value,
             slab_colour = as.factor(median_charadriiformes_wild),
             slab_fill = as.factor(as.factor(median_charadriiformes_wild)))) +
  stat_halfeye(slab_alpha = 0.7)  
  
  
  
# Results are averaged over the levels of: season
diffusionmodel1_fit %>%
  emmeans(~ collection_regionname,
          epred = TRUE) %>% 
  contrast(method = "revpairwise") %>% 
  gather_emmeans_draws()


  emtrends(~ collection_regionname,
           var = "collection_regionname",
           at = list(collection_regionname = unique(diffusion_data$collection_regionname)),
           epred = TRUE, 
           dpar = "mu") %>% 
  gather_emmeans_draws() %>%
  ggplot(aes(x  = .value,
             slab_colour = as.factor(median_anseriformes_wild),
             slab_fill = as.factor(as.factor(median_anseriformes_wild)))) +
  stat_halfeye(slab_alpha = 0.7)          
             
  
  
  
# Combine plots together
align_plots(align = 'v', axis = 'l')

