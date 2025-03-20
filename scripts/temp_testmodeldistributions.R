diffusionmodel1_fit %>%
  emmeans(~ collection_regionname,
          # epred = TRUE, 
          #dpar = "mu",
          re_formula = NA , 
          regrid = "response",
          tran  = "log", 
          type = "response",
          allow_new_levels = TRUE)#


diffusionmodel1_fit %>%
  emmeans(~ collection_regionname,
          re_formula = NA , 
          epred = TRUE,
          at = list(median_anseriformes_wild = 0.6224637,
                    median_charadriiformes_wild = 0.08553386))




diffusionmodel_simple3_fit <- brm(
  bf(weighted_diff_coeff ~ 0 + collection_regionname ) ,
  data = diffusion_data,
  family = Gamma(link = "log"),
  chains = CHAINS,
  cores = CORES, 
  threads = 2, 
  backend = "cmdstanr",
  iter = ITER,
  warmup = BURNIN,
  seed = SEED,
  control = list(adapt_delta = 0.95))



diffusionmodel_simple3_fit %>%
  linpred_draws(tibble(collection_regionname = unique(diffusion_data$collection_regionname))) %>%
  ggplot(aes(x  = .linpred,
             y = collection_regionname,
             slab_colour = collection_regionname,
             slab_fill = collection_regionname)) +
  stat_halfeye(slab_alpha = 0.7,
               p_limits = c(0.01, 0.99),
               point_interval = "median_hdi",
               linewidth = 1.5,
               .width =  0.95)+
  scale_x_continuous(expression(paste('Weighted Diffusion Coefficient (',Km**2~year**-1, ')' )),
                     breaks = log1p(c(0, 10^(seq(from = 1, to = 10, by = 1)))),
                     labels = expression(0, 1%*%10^1, 1%*%10^2, 1%*%10^3, 1%*%10^4,1%*%10^5,1%*%10^6,1%*%10^7,1%*%10^8, 1%*%10^9,1%*%10^10),
                     limits = c(-0.01, log1p(10^10.5)),
                     expand = c(0.02,0.02))

diffusionmodel_simple3_fit %>%
  epred_draws(tibble(collection_regionname = unique(diffusion_data$collection_regionname))) %>%
  ggplot(aes(x  = log(.epred),
             y = collection_regionname,
             slab_colour = collection_regionname,
             slab_fill = collection_regionname)) +
  stat_halfeye(slab_alpha = 0.7,
               p_limits = c(0.01, 0.99),
               point_interval = "median_hdi",
               linewidth = 1.5,
               .width =  0.95)+
  scale_x_continuous(expression(paste('Weighted Diffusion Coefficient (',Km**2~year**-1, ')' )),
                     breaks = log1p(c(0, 10^(seq(from = 1, to = 10, by = 1)))),
                     labels = expression(0, 1%*%10^1, 1%*%10^2, 1%*%10^3, 1%*%10^4,1%*%10^5,1%*%10^6,1%*%10^7,1%*%10^8, 1%*%10^9,1%*%10^10),
                     limits = c(-0.01, log1p(10^10.5)),
                     expand = c(0.02,0.02))

diffusionmodel_simple3_fit %>%
  predicted_draws(tibble(collection_regionname = unique(diffusion_data$collection_regionname))) %>%
  ggplot(aes(x  = log(.prediction),
             y = collection_regionname,
             slab_colour = collection_regionname,
             slab_fill = collection_regionname)) +
  stat_halfeye(slab_alpha = 0.7,
               p_limits = c(0.01, 0.99),
               point_interval = "median_hdi",
               linewidth = 1.5,
               .width =  0.95)+
  scale_x_continuous(expression(paste('Weighted Diffusion Coefficient (',Km**2~year**-1, ')' )),
                     breaks = log1p(c(0, 10^(seq(from = 1, to = 10, by = 1)))),
                     labels = expression(0, 1%*%10^1, 1%*%10^2, 1%*%10^3, 1%*%10^4,1%*%10^5,1%*%10^6,1%*%10^7,1%*%10^8, 1%*%10^9,1%*%10^10),
                     limits = c(-0.01, log1p(10^10.5)),
                     expand = c(0.02,0.02))


diffusionmodel1_fit_gamma_5 %>%
  avg_predictions( by = 'collection_regionname')

diffusionmodel1_fit_gamma_5 %>%
  emmeans(~ collection_regionname,
          epred = T)

diffusionmodel1_fit%>%
  emmeans(~ collection_regionname,
          # epred = TRUE, 
          #dpar = "mu",
          re_formula = NA , 
          regrid = "response",
          tran  = "log", 
          type = "response",
          allow_new_levels = TRUE) %>%
  gather_emmeans_draws() %>%
  ggplot(aes(x  = .value,
             y = collection_regionname,
             slab_colour = collection_regionname,
             slab_fill = collection_regionname)) +
  stat_halfeye(slab_alpha = 0.7,
               p_limits = c(0.01, 0.99),
               point_interval = "median_hdi",
               linewidth = 1.5,
               .width =  0.95)+
  scale_x_continuous(expression(paste('Weighted Diffusion Coefficient (',Km**2~year**-1, ')' )),
                     breaks = log1p(c(0, 10^(seq(from = 1, to = 10, by = 1)))),
                     labels = expression(0, 1%*%10^1, 1%*%10^2, 1%*%10^3, 1%*%10^4,1%*%10^5,1%*%10^6,1%*%10^7,1%*%10^8, 1%*%10^9,1%*%10^10),
                     limits = c(-0.01, log1p(10^10.5)),
                     expand = c(0.02,0.02))





diffusionmodel1_fit_gamma_5%>%
  emmeans(~ collection_regionname,
          epred = TRUE, 
          #dpar = "mu",
          re_formula = NA ) %>%
  gather_emmeans_draws() %>%
  ggplot(aes(x  = log(.value),
             y = collection_regionname,
             slab_colour = collection_regionname,
             slab_fill = collection_regionname)) +
  stat_halfeye(slab_alpha = 0.7,
               p_limits = c(0.01, 0.99),
               point_interval = "median_hdi",
               linewidth = 1.5,
               .width =  0.95)+
  scale_x_continuous(expression(paste('Weighted Diffusion Coefficient (',Km**2~year**-1, ')' )),
                     breaks = log1p(c(0, 10^(seq(from = 1, to = 10, by = 1)))),
                     labels = expression(0, 1%*%10^1, 1%*%10^2, 1%*%10^3, 1%*%10^4,1%*%10^5,1%*%10^6,1%*%10^7,1%*%10^8, 1%*%10^9,1%*%10^10),
                     limits = c(-0.01, log1p(10^10.5)),
                     expand = c(0.02,0.02))
  
diffusionmodel1_fit_gamma_5%>%
  emmeans(~ collection_regionname,
          epred = TRUE, 
          #dpar = "shape",
          re_formula = NA ) %>%
  pairs()


diffusion_data %>%
  add_residual_draws(diffusionmodel1_fit_gamma) %>%
  ggplot(aes(x = .row, y = .residual)) +
  stat_pointinterval() + 
  facet_grid(
    rows = vars(collection_regionname),
    labeller =  labeller(collection_regionname=str_to_title)) 


diffusionmodel1_fit_gamma_5 %>%
  predicted_draws(newdata = expand_grid(collection_regionname = unique(as.character(diffusion_data$collection_regionname)),
                                        collection_season = unique(diffusion_data$collection_season)) %>%
                    drop_na() %>%
                    
                    mutate(median_anseriformes_wild = median(diffusion_data$median_anseriformes_wild),
                           median_charadriiformes_wild = median(diffusion_data$median_charadriiformes_wild)),
                  re_formula = NA) %>%
  drop_na(collection_regionname) %>%
  ggplot() + 
  geom_histogram(aes(x = log1p(.prediction), y = after_stat(density), colour = collection_regionname, fill = collection_regionname), 
                 alpha = 0.7, binwidth = 0.5) + 
  scale_colour_manual(values = region_colours)+
  scale_fill_manual(values = region_colours) +
  scale_x_continuous(expression(paste('Predicted Weighted Diffusion Coefficient (',Km**2~year**-1, ')' )),
                     breaks = log1p(c(0, 10^(seq(from = 1, to = 10, by = 1)))),
                     labels = expression(0, 1%*%10^1, 1%*%10^2, 1%*%10^3, 1%*%10^4,1%*%10^5,1%*%10^6,1%*%10^7,1%*%10^8, 1%*%10^9,1%*%10^10),
                     limits = c(-0.01, log1p(10^10.5)),
                     expand = c(0.02,0.02)) + 
  scale_y_continuous('Probability Density' ,
                     expand = c(0,0))+
  facet_grid(
    rows = vars(collection_regionname),
    labeller =  labeller(collection_regionname=str_to_title)) +
  theme(legend.position = 'none')
  
  geom_vline(aes(xintercept = emmean, colour = collection_regionname), data = averages, linetype = 'dashed') +
  geom_text(aes(label =  paste0("E*'('*X*'|'*X*'>'*0*') = '*", label, "~km^2"), 
                colour = collection_regionname),
            parse = T,
            x = 17.5, 
            y = 0.1,
            size = 2.5,
            data = averages) + 
  global_theme + 
  theme(strip.placement  = 'inside')


posteriorpredictive <-pp_check(diffusionmodel1_fit_gamma, ndraws = 500,type = 'stat_grouped', group = 'collection_regionname' )
posteriorpredictive$data  %>%
  mutate(value = log1p(value)) %>%
  ggplot() + 
  geom_density(aes(x = value, 
                   group= rep_id, 
                   alpha = is_y,  
                   linewidth = is_y, 
                   colour = is_y), 
               key_glyph = draw_key_path) + 
  scale_alpha_manual(values = c('FALSE' = 0.00001, 
                                'TRUE' = 1)) + 
  scale_linewidth_manual(values = c('FALSE' = 0.1, 
                                    'TRUE' = 1),
                         labels  = c(expression(italic('y')['rep']),
                                     expression(italic('y'))))+ 
  scale_colour_manual(values = c('FALSE' = '#cbc9e2', 
                                 'TRUE' = '#54278f'),
                      labels  = c(expression(italic('y')['rep']),
                                  expression(italic('y')))) + 
  
  guides(alpha= 'none', 
         colour=guide_legend()) + 
  scale_x_continuous(expression(paste('Predicted Weighted Diffusion Coefficient (',Km**2~year**-1, ')' )),
                     breaks = log1p(c(0, 10^(seq(from = 2, to = 10, by = 4)))),
                     labels = expression(0,  1%*%10^2,  1%*%10^6, 1%*%10^10),
                     #limits = c(-0.01, log1p(10^9)),
                     expand = c(0.02,0.02))+
  scale_y_continuous('Density',expand = c(0,0)) 