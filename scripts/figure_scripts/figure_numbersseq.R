

inv_logit <- function(x){
  return(exp(x)/(1+exp(x)))
}

continent_specific_detection_samples <- numbers_model_2 %>%
  gather_draws(., continent_specific_detection[i]) 

beta_sequences_samples <-  numbers_model_2 %>%
  gather_draws(., beta_sequences) %>%
  ungroup() %>% 
  dplyr::select(.draw,
         .chain, 
         .iteration, 
         beta_sequences = .value)

year_detection_samples <- numbers_model_2 %>%
  gather_draws(., year_detection[i]) %>% 
  dplyr::select(.draw,
         .chain,
         .iteration,
         year_detection = .value)


detection_probabilities <- continent_specific_detection_samples %>%
  rename(continent_specific_detection = .value) %>%
  inner_join(beta_sequences_samples ,
             by = c(".draw", ".chain", ".iteration")) %>%
  inner_join(year_detection_samples, 
             by = c(".draw", ".chain", ".iteration"), 
             relationship = "many-to-many") %>%
  mutate(sequences = list(log1p(1:75))) %>%
  unnest(sequences) %>%
  mutate(
    probability = inv_logit(continent_specific_detection + beta_sequences * sequences + year_detection)
  )



detection_probabilities %>%
  ggplot(aes(x = expm1(sequences), y = probability)) +
  stat_lineribbon(point_interval = "median_hdci", 
                  alpha = 0.7) + 
  
  scale_y_continuous('P(Detection)',
                     breaks = seq(0.2, 1, by = 0.2),
                     labels = seq(0.2, 1, by = 0.2),
                     expand = c(0,0))+ 
  scale_x_continuous('Sequences (n)',
                     expand = c(0,0),
                     breaks = seq(0,75,by=25)) +
  scale_fill_brewer() + 
  global_theme+ 
  theme(strip.placement  = 'inside',
        strip.text = element_text(face = 'bold', size = 10),
        strip.background = element_blank(),
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 9),
        legend.text = element_text(size = 8))


year_abundance_samples <- numbers_model_2 %>%
  gather_draws(., year_abundance[i]) %>% 
  dplyr::select(.draw,
                .chain,
                .iteration,
                year_detection = .value) 
  
year_abundance_samples %>%
  ggplot(aes(x = year_detection, y = as.factor(i))) + 
  stat_slabinterval(point_interval = "median_hdci")
