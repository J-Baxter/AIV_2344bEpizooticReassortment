# interpretation of stan model
diffusionmodel1_fit_gamma_19 %>% 
  epred_draws(newdata = diffusion_data) %>%  
  group_by(.draw, collection_regionname) %>% 
  summarise(avg_epred = mean(.epred), .groups = "drop") %>% 
  group_by(collection_regionname) %>%
  median_hdci(avg_epred)


#1) continent stratified estimates of the observed number of reassortants per month
test_pred <- as_draws_df(numbers_model_2$draws('y_rep') ) %>%
  pivot_longer(cols = starts_with("y_rep["), names_to = "row", values_to = ".epred") %>%
  mutate(row = as.integer(str_extract(row, "\\d+"))) %>%
  left_join(data_processed_2 %>% rowid_to_column('row'),
            by = 'row')

test_pred %>%
  group_by(.draw, collection_regionname) %>% 
  summarise(avg_epred = mean(.epred), .groups = "drop") %>% 
  group_by(collection_regionname) %>%
  median_hdci(avg_epred)


#1) continent stratified estimates of the true (latent) number of reassortants per month
test_pred_2 <- as_draws_df(numbers_model_2$draws('N_rep') ) %>%
  pivot_longer(cols = starts_with("N_rep["), names_to = "row", values_to = ".epred") %>%
  mutate(row = as.integer(str_extract(row, "\\d+"))) %>%
  left_join(data_processed_2 %>% rowid_to_column('row'),
            by = 'row')

test_pred_2 %>%
  group_by(.draw, collection_regionname) %>% 
  summarise(avg_epred = mean(.epred), .groups = "drop") %>% 
  group_by(collection_regionname) %>%
  median_hdci(avg_epred)


rpois(100, exp(1.77 + -0.06697910 * 0.000000 + -0.03163115 ))

#2) ZI 
t <- get_variables(numbers_model_2)
# Trace plot
numbers_model_2 %>%
  gather_draws(., !!!syms(t)) %>%
  filter(grepl('theta', .variable)) %>%
  
  group_by(.draw, .variable) %>% 
  summarise(avg_value = mean(.value), .groups = "drop") %>% 
  group_by(.variable) %>%
  median_hdci(avg_value)


#2) Probability of detection


#3) Effect of number of sequences on p(detect)


#4) Effect of number of cases on N


#5) Grouped average across years



draws <- as_draws_df(numbers_model_2$draws())
newdata <- tibble(x = seq(-2, 2, length.out = 50))

# Add fitted draws
preds <- draws %>%
  select(.draw, alpha, beta) %>%
  crossing(newdata) %>%
  mutate(
    y_hat = alpha + beta * x
  )

avg_preds <- preds %>%
  group_by(x) %>%
  summarise(
    estimate = mean(y_hat),
    lower = quantile(y_hat, 0.025),
    upper = quantile(y_hat, 0.975)
  )
