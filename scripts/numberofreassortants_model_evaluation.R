

###################################### MCMC Diagnostics #############################################
t <- get_variables(numbers_model)
# Trace plot
numbers_model %>%
  gather_draws(., !!!syms(t)) %>%
  filter(grepl('^continent|^beta|^theta|^lp', .variable)) %>%
  ggplot(aes(x = .iteration,
             y = .value, 
             col = as.factor(.chain)))+
  geom_line(alpha = 0.8) + 
  facet_wrap(~ .variable,
             ncol = 2,
             scale  = 'free_y',
             strip.position = 'left')+
  scale_colour_brewer(palette = 'GnBu', 'Chains') +
  theme_minimal(base_size = 8) + 
  theme(legend.position = 'bottom')

# Autocorrelation
stan_acf <- posterior::as_draws_array(numbers_model, nchains = 4) %>% 
  mcmc_acf()

stan_acf$data %>% 
  filter(!grepl('y_rep', Parameter)) %>%
  ggplot(aes(y = AC, 
             x = Lag,
             colour = as.factor(Chain))) +
  geom_path() + 
  facet_wrap(~Parameter, labeller = label_parsed, ncol = 4) + 
  theme_minimal()  + 
  scale_colour_brewer('Chain') +
  theme(legend.position = 'bottom',
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 10),
        legend.text = element_text(size = 8))

ggsave('~/Downloads/flu_plots/class_autocorrelation.jpeg',
       dpi = 360,
       height = 29,
       width = 21,
       units = 'cm')

 # Might be worth checking if any of the params are at the 'worry about' threshold neff/n



###################################### Posterior Checks ############################################
y_rep_matrix <- numbers_model$draws('y_rep') %>%
  posterior::as_draws_matrix()

# Global
ppc_bars(y =  data_processed_2 %>% pull(n_reassortants),
         yrep = y_rep_matrix[sample(0:16000, 100),]) %>% 
  .$data %>%
  ggplot() + 
  geom_bar(aes(x = x,
               y = y_obs),
           stat = 'identity',
           fill = '#cbc9e2',
           colour = '#cbc9e2',
           alpha = 0.7) + 
  geom_pointinterval(aes(x = x,
                         y=m,
                         ymin = l,
                         ymax = h), orientation = 'x',
                     colour = '#54278f')

# Grouped by Continent
ppc_bars_grouped(y =  data_processed_2 %>% pull(n_reassortants), 
                 group = data_processed_2 %>% pull(collection_regionname) %>% as.factor(),
                 yrep = y_rep_matrix[sample(0:16000, 250),]) %>%
  .$data %>%
  #mutate(group = case_when(group == 1 ~ 'Africa',
                           #group == 2 ~ 'Asia',
                           #group == 3 ~ 'Northern and Central America',
                           #group == 4 ~ 'Europe') ) %>%
  
  ggplot() + 
  geom_bar(aes(x = x,
               y = y_obs),
           stat = 'identity',
           fill = '#cbc9e2',
           colour = '#cbc9e2',
           alpha = 0.7) + 
  geom_pointinterval(aes(x = x,
                         y=m,
                         ymin = l,
                         ymax = h), orientation = 'x',
                     colour = '#54278f') +
  
  facet_wrap(~group, labeller  = as_labeller(str_to_title)) + 
  scale_x_discrete('Reassortants per Month')+
  scale_y_continuous('Count', expand = c(0,0)) +
  theme_classic() + 
  theme(strip.text = element_text(face = 'bold', size = 10),
        strip.background = element_blank(),
        legend.title = element_blank(),
        legend.position = 'bottom',
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 10),
        legend.text = element_text(size = 8))


ggsave('~/Downloads/flu_plots/class_ppc.jpeg',
       dpi = 360,
       device = 'jpeg' ,
       height = 12,
       width = 17, 
       units = 'cm')


# Identifiability (Prior and Posterior Plots)
test_fit_linear %>%
  gather_draws(., !!!syms(t)) %>%
  mutate(type = 'posterior') %>%
  filter(grepl('^continent_specific_theta', .variable)) %>%
  ggplot() + 
  
  geom_density(aes(x = .value#,
                   #y = after_stat(density)
  ),
  inherit.aes = F, 
  fill = '#1b9e77') +
  
  stat_function(fun = dbeta,
                args = list(shape1 = 2, shape2 =5),
                fill = '#d95f02',
                geom = 'area',
                alpha = 0.5) +
  
  
  facet_wrap(~`.variable`)



###################################### Residuals Checks ############################################

simulated_residuals <- createDHARMa(
  simulatedResponse = t(y_rep_matrix[sample(0:16000, 250),]),
  observedResponse = data_processed_2 %>% pull(n_reassortants)
)


# QQ plot
qq_data <- data.frame(
  sample = sort(simulated_residuals$scaledResiduals),
  theoretical = sort(ppoints(length(simulated_residuals$scaledResiduals)))
)

ggplot(qq_data, aes(sample = sample)) +
  stat_qq(distribution = stats::qunif) +
  stat_qq_line(distribution = stats::qunif) +
  scale_y_continuous('Observed')+
  scale_x_continuous('Expected')+
  theme_minimal() +
  theme(axis.text = element_text(size = 8),
        axis.title = element_text(size = 10),
        legend.text = element_text(size = 8))

ggsave('~/Downloads/flu_plots/numbers_qq.jpeg',
       dpi = 360,
       height = 15,
       width = 15,
       units = 'cm')