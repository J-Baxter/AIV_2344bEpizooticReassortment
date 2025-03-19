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


diffusionmodel_simple_fit %>%
  avg_predictions( by = 'collection_regionname')

diffusionmodel_simple3_fit %>%
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





diffusionmodel1_fit_gamma%>%
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
  
  

