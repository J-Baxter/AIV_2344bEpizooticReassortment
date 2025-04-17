####################################################################################################
####################################################################################################
## Script name:
##
## Purpose of script:
##
## Date created: 2025-04-15
##
##
########################################## SYSTEM OPTIONS ##########################################
options(scipen = 6, digits = 7) 
memory.limit(30000000) 

  
########################################## DEPENDENCIES ############################################
# Packages
library(tidyverse)
library(magrittr)
library(cmdstanr)
library(posterior)
library(bayesplot)


# User functions


############################################## DATA ################################################



####################################### Test basic model ############################################
# Compile the model
basic_mod <- cmdstan_model('./scripts/stan_models/n_reassortants_basic.stan')

real_data <- list(N = nrow(data_processed),
                  K = data_processed %>% pull(n_reassortants) %>% max(),
                  C = data_processed %>% pull(collection_regionname) %>% n_distinct(),
                  y = data_processed %>% pull(n_reassortants),
                  continent = data_processed %>% pull(collection_regionname) %>% as.factor() %>% as.numeric())


# Run the model
test_fit_new <- basic_mod$sample(
  data = real_data,
  seed = 42,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 500,
  iter_sampling = 2000
)


print(test_fit_new$summary())

# posterior predictive check
y_rep_matrix <- test_fit_new$draws('y_rep') %>%
  posterior::as_draws_matrix()

ppc_bars(y =  data_processed %>% pull(n_reassortants), yrep = y_rep_matrix[sample(500:2000, 100),])


t <- get_variables(test_fit_new)
# Trace plot
test_fit_new %>%
  gather_draws(., !!!syms(t)) %>%
  filter(grepl('^lambda|^p|^theta|^lp', .variable)) %>%
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


# Trace - Rank plot
as_draws_df(test_fit_new) %>% 
  mcmc_rank_overlay(regex_pars = '^lambda|^p|^theta|^lp') %>%
  .$data %>%
  
  ggplot() + 
  geom_step(aes(x = bin_start, 
                y = n,
                colour = chain),
            linewidth = 0.8) + 
  scale_y_continuous(expand = c(0,0)) +
  scale_x_continuous(expand = c(0,0)) + 
  facet_wrap(~parameter,  ncol = 3) +
  theme_minimal(base_size = 8) + 
  
  scale_colour_brewer(palette = 'GnBu', 'Chains') 


# Autocorrelation
stan_acf <- posterior::as_draws_array(test_fit_new, nchains = 4) %>% 
  mcmc_acf()

stan_acf$data %>% 
  filter(!grepl('y_rep', Parameter)) %>%
  ggplot(aes(y = AC, 
             x = Lag,
             colour = as.factor(Chain))) +
  geom_path() + 
  facet_wrap(~Parameter) + 
  theme_classic() + 
  scale_color_brewer()


# Plot continent specific lambda (and priors)
test_fit_new %>%
  gather_draws(., !!!syms(t)) %>%
  mutate(type = 'posterior') %>%
  filter(grepl('^lambda', .variable)) %>%
  ggplot() + 
  geom_histogram(aes(x = .value,
                     y = after_stat(density)),
                 inherit.aes = F, 
                 binwidth = 0.1, 
                 fill = '#1b9e77') +
  
  stat_function(fun = dnorm,
                args = list(mean = 3, sd = 1.5),
                fill = '#d95f02',
                geom = 'area',
                alpha = 0.5) +
  
  facet_grid(rows = vars(`.variable`))


# Plot continent specific detection (and priors)
test_fit_new %>%
  gather_draws(., !!!syms(t)) %>%
  mutate(type = 'posterior') %>%
  filter(grepl('^p', .variable)) %>%
  ggplot() + 
  # geom_histogram(aes(x = .value,
  #   y = after_stat(density)),
  #  inherit.aes = F, 
  #binwidth = 0.1, 
  #fill = '#1b9e77') +
  
  geom_density(aes(x = .value#,
                   #y = after_stat(density)
  ),
  inherit.aes = F, 
  binwidth = 0.1, 
  fill = '#1b9e77') +
  
  stat_function(fun = dbeta,
                args = list(shape1 = 2, shape2 = 2),
                fill = '#d95f02',
                geom = 'area',
                alpha = 0.5) +
  
  facet_wrap(~`.variable`)


################################### Test full linear model #########################################
# Compile the model
linear_mod <- cmdstan_model('./scripts/stan_models/n_reassortants_full_linear.stan')

linear_data <- list(N = nrow(data_processed),
                  K = data_processed %>% pull(n_reassortants) %>% max(),
                  C = data_processed %>% pull(collection_regionname) %>% n_distinct(),
                  y = data_processed %>% pull(n_reassortants),
                  continent = data_processed %>% pull(collection_regionname) %>% as.factor() %>% as.numeric(),
                  cases =  data_processed %>% pull(woah_susceptibles_log1p),
                  sequences =  data_processed %>% pull(n_sequences))


# Run the model
test_fit_linear <- linear_mod$sample(
  data = linear_data,
  seed = 42,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 1000,
  iter_sampling = 5000
)

print(test_fit_linear$summary())

# posterior predictive check
y_rep_matrix <- test_fit_linear$draws('y_rep') %>%
  posterior::as_draws_matrix()

ppc_bars(y =  data_processed %>% pull(n_reassortants), yrep = y_rep_matrix[sample(500:2000, 100),])


t <- get_variables(test_fit_linear)
# Trace plot
test_fit_linear %>%
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


# Trace - Rank plot
as_draws_df(test_fit_linear) %>% 
  mcmc_rank_overlay(regex_pars = '^continent|^beta|^theta|^lp') %>%
  .$data %>%
  
  ggplot() + 
  geom_step(aes(x = bin_start, 
                y = n,
                colour = chain),
            linewidth = 0.8) + 
  scale_y_continuous(expand = c(0,0)) +
  scale_x_continuous(expand = c(0,0)) + 
  facet_wrap(~parameter,  ncol = 3) +
  theme_minimal(base_size = 8) + 
  
  scale_colour_brewer(palette = 'GnBu', 'Chains') 


# Autocorrelation
stan_acf <- posterior::as_draws_array(test_fit_linear, nchains = 4) %>% 
  mcmc_acf()

stan_acf$data %>% 
  filter(!grepl('y_rep', Parameter)) %>%
  ggplot(aes(y = AC, 
             x = Lag,
             colour = as.factor(Chain))) +
  geom_path() + 
  facet_wrap(~Parameter) + 
  theme_classic() + 
  scale_color_brewer()


# Plot continent specific lambda (and priors)
test_fit_linear %>%
  gather_draws(., !!!syms(t)) %>%
  mutate(type = 'posterior') %>%
  filter(grepl('^continent_specific_abundance', .variable)) %>%
  ggplot() + 
  geom_histogram(aes(x = .value,
                     y = after_stat(density)),
                 inherit.aes = F, 
                 binwidth = 0.1, 
                 fill = '#1b9e77') +
  
  stat_function(fun = dnorm,
                args = list(mean = 3, sd = 1.5),
                fill = '#d95f02',
                geom = 'area',
                alpha = 0.5) +
  
  facet_grid(rows = vars(`.variable`))


# Plot continent specific detection (and priors)
test_fit_linear %>%
  gather_draws(., !!!syms(t)) %>%
  mutate(type = 'posterior') %>%
  filter(grepl('^continent_specific_detection', .variable)) %>%
  ggplot() + 

  geom_density(aes(x = .value#,
                   #y = after_stat(density)
  ),
  inherit.aes = F, 
  fill = '#1b9e77') +
  
  #stat_function(fun = dbeta,
                #args = list(shape1 = 1.5, shape2 =1.5),
                #fill = '#d95f02',
                #geom = 'area',
                #alpha = 0.5) +
  stat_function(fun = dunif,
                args = list(min = 0, max= 1),
                fill = '#d95f02',
                geom = 'area',
                alpha = 0.5) +
  
  
  facet_wrap(~`.variable`)



test_fit_linear %>%
  gather_draws(., !!!syms(t)) %>%
  mutate(type = 'posterior') %>%
  filter(grepl('^beta_sequences', .variable)) %>%
  ggplot() + 

  geom_density(aes(x = .value#,
                   #y = after_stat(density)
  ),
  inherit.aes = F, 
  fill = '#1b9e77') +
  
  stat_function(fun = dnorm,
                args = list(mean = 0, sd = 1),
                fill = '#d95f02',
                geom = 'area',
                alpha = 0.5) +

  facet_wrap(~`.variable`)


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

################################### Test full linear + R Eff model #########################################
# Compile the model
raneff_mod <- cmdstan_model('./scripts/stan_models/n_reassortants_single_raneff.stan')

raneff_data <- list(N = nrow(data_processed_2),
                    y = data_processed_2 %>% pull(n_reassortants),
                    K = data_processed_2 %>% pull(n_reassortants) %>% max(),
                    C = data_processed_2 %>% pull(collection_regionname) %>% n_distinct(),
                    Y = data_processed_2 %>% pull(collection_year) %>% n_distinct(),
                    continent_index = data_processed_2 %>% pull(collection_regionname) %>% as.factor() %>% as.numeric(),
                    year_index = data_processed_2 %>% pull(collection_year) %>% as.factor() %>% as.numeric(),
                    cases =  data_processed_2 %>% pull(woah_susceptibles_monthly_log1p),
                    sequences =  data_processed_2 %>% pull(n_sequences_log1p))


# Run the model (00:24)
test_fit_raneff <- raneff_mod$sample(
  data = raneff_data,
  seed = 42,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 500,
  iter_sampling = 2000
)

raneff_parms <- test_fit_raneff$summary()

# posterior predictive check
y_rep_matrix <- test_fit_raneff$draws('y_rep') %>%
  posterior::as_draws_matrix()

ppc_bars(y =  data_processed_2 %>% pull(n_reassortants), yrep = y_rep_matrix[sample(500:8000, 100),])



simulated_residuals <- createDHARMa(
  simulatedResponse = t(y_rep_matrix[sample(500:8000, 100),]),
  observedResponse = data_processed_2 %>% pull(n_reassortants)
)



# QQ plot
qq_data <- data.frame(
  sample = sort(simulated_residuals$scaledResiduals),
  theoretical = sort(ppoints(length(simulated_residuals$scaledResiduals)))
)

ggplot(qq_data, aes(sample = sample)) +
  stat_qq(distribution = stats::qunif) +
  stat_qq_line(distribution = stats::qunif) 

############################################## END #################################################
####################################################################################################
####################################################################################################