stan_code <- 'data {
  int<lower=0> N; // Number of observations
  array[N] int<lower=0> y; // Observed counts
  int<lower=0> K; // Upper bound of population size
  int<lower=1> C; // Number of continents
  array[N] int<lower=1, upper=C> continent; // Continent index for each observation
}

parameters {
  real<lower=0, upper=1> theta; // Zero-inflation probability
  array[C] real<lower=0> lambda; // Poisson rate stratified by continent
  array[C] real<lower=0, upper=1> p; // Detection probability stratified by continent
}

model {
  // Priors
  lambda ~ normal(3, 1.5);
  theta ~ beta(2, 5); // Example prior
  p ~ beta(2, 2); // Example prior
  
  // Loop over number of observations
  for (i in 1:N) {
    vector[K] lp;
    int c = continent[i]; // Current continent

    // Loop over plausible values of K to marginalize out discrete latent variables
    for (j in 1:K) {
      int current_population = y[i] + j - 1;
      
      // Likelihood for y[i] = 0
      if (y[i] == 0)
        lp[j] = log_sum_exp(to_vector({
          bernoulli_lpmf(1 | theta), 
          bernoulli_lpmf(0 | theta) + poisson_lpmf(current_population | lambda[c]),
          bernoulli_lpmf(0 | theta) + poisson_lpmf(current_population | lambda[c]) + binomial_lpmf(y[i] | current_population, p[c])
        }));
      
      // Likelihood for y[i] > 0
      else
        lp[j] = bernoulli_lpmf(0 | theta) + poisson_lpmf(current_population | lambda[c]) + binomial_lpmf(y[i] | current_population, p[c]);
    }
    
    target += log_sum_exp(lp); // Aggregate the probabilities correctly
  }
}'

# Compile the model
mod <- cmdstan_model(write_stan_file(stan_code))

# Define parameters
N <- 300  # Number of observations
K <- 20   # Upper bound of population size
C <- 3    # Number of continents

# True parameters for each continent
true_lambda <- c(2.5, 5.0, 3.5)
true_p <- c(0.1, 0.2, 0.3)
theta <- 0.2  # Zero-inflation probability

# Simulate data
continent <- sample(1:C, N, replace = TRUE)
y <- integer(N)

for (i in 1:N) {
  c <- continent[i]
  if (runif(1) < theta) {
    y[i] <- 0  # Zero-inflated case
  } else {
    # Simulate from Poisson and Binomial
    population_size <- rpois(1, true_lambda[c])
    y[i] <- rbinom(1, population_size, true_p[c])
  }
}

# Print the simulated data
data_list <- list(N = N, K = K, C = C, y = y, continent = continent)
real_data <- list(N = nrow(data_processed),
                  K = data_processed %>% pull(n_reassortants) %>% max(),
                  C = data_processed %>% pull(collection_regionname) %>% n_distinct(),
                  y = data_processed %>% pull(n_reassortants),
                  continent = data_processed %>% pull(collection_regionname) %>% as.factor() %>% as.numeric())


# Run the model
test_fit <- mod$sample(
  data = real_data,
  seed = 42,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 500,
  iter_sampling = 2000
)

print(test_fit$summary())

t <- get_variables(test_fit)

test_fit %>%
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



test_fit %>%
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


ppc_dens_overlay(y = fit$y,
                 yrep = posterior_predict(fit, draws = 50))  

stan_acf <- posterior::as_draws_array(test_fit, nchains = 4) %>% 
  mcmc_acf()

stan_acf$data %>% 
  ggplot(aes(y = AC, 
             x = Lag,
             colour = as.factor(Chain))) +
  geom_path() + 
  facet_wrap(~Parameter) + 
  theme_classic() + 
  scale_color_brewer()
