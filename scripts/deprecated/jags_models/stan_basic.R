stan_code <- 'data {
  int<lower=0> N; // Number of observations
  array[N] int<lower=0> y; // observed counts
  int<lower=0> K; // Upper bound of population size
}

parameters {
  real<lower=0, upper=1> theta; // zero-inflation probability
  real<lower=0> lambda; // Poisson rate
  real<lower=0, upper=1> p; // Detection probability
}

model {
  // Priors
  lambda ~ normal(3, 1);
  theta ~ beta(2, 5); // Example prior
  p ~ beta(2, 5); // Example prior
  
  // Loop over number of observations
  for (i in 1:N) {
    vector[K] lp;
    
    // Loop over plausible values of K to marginalize out discrete latent variables
    for (j in 1:K) {
      int current_population = y[i] + j - 1;
      
      // Likelihood for y[i] = 0
      if (y[i] == 0)
        lp[j] = log_sum_exp(to_vector({
          bernoulli_lpmf(1 | theta), 
          bernoulli_lpmf(0 | theta) + poisson_lpmf(current_population | lambda),
          bernoulli_lpmf(0 | theta) + poisson_lpmf(current_population | lambda) + binomial_lpmf(y[i] | current_population, p)
        }));
      
      // Likelihood for y[i] > 0
      else
        lp[j] = bernoulli_lpmf(0 | theta) + poisson_lpmf(current_population | lambda) + binomial_lpmf(y[i] | current_population, p);
    }
    
    target += log_sum_exp(lp); // Aggregate the probabilities correctly
  }
}
'

# Compile the model
mod <- cmdstan_model(write_stan_file(stan_code))

# Define some example data
true_lambda <- 3.0   # Poisson rate (expected population size)
true_p <- 0.4        # Detection probability
true_theta <- 0.2    # Zero-inflation probability
N <- 150             # Number of observations
K <- 10             # Upper bound on total population size

# Simulate the population size for each observation
true_population <- rpois(N, lambda = true_lambda)

# Simulate observed counts with binomial thinning (detection probability)
y <- rbinom(N, size = true_population, prob = true_p)

# Introduce zero inflation (some zeros due to `theta`)
zero_inflation <- rbinom(N, size = 1, prob = true_theta)
y[zero_inflation == 1] <- 0  # Force some observations to zero

# Create test data list for Stan
test_data <- list(N = N, y = y, K = 10)

# Run the model
fit <- mod$sample(
  data = test_data,
  seed = 42,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 500,
  iter_sampling = 2000
)

print(fit$summary())

true_theta <- 0.2    # Zero-inflation probability 0.153
true_lambda <- 3.0   # Poisson rate (expected population size) 1.6
true_p <- 0.5        # Detection probability 0.928


# Print a summary of the results
print(fit$summary())