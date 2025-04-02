set.seed(42)
N <- 100
lambda <- 2.5
pi <- 0.3
p <- 0.6

# Simulate true event counts with zero inflation
true_counts <- rpois(N, lambda)
inflated_zeros <- rbinom(N, 1, pi)
true_counts[inflated_zeros == 1] <- 0

# Simulate observed counts based on detection
y <- sapply(true_counts, function(n) rbinom(1, n, p))

# Show a portion of the data
head(y)


# Data list
data_list <- list(
  N = N,
  y = y
)

# Initial values
init_list <- list(
  list(lambda = 1.0, pi = 0.2, p = 0.5),
  list(lambda = 1.5, pi = 0.3, p = 0.6),
  list(lambda = 2.0, pi = 0.4, p = 0.7),
  list(lambda = 2.5, pi = 0.5, p = 0.8)
)

# Fit the model
fit <- stan(
  model_code = 'data {
  int<lower=0> N; // number of observations
  int<lower=0> y[N]; // observed counts
}

parameters {
  real<lower=0> lambda; // Poisson rate for true number of events
  real<lower=0, upper=1> pi; // zero-inflation probability
  real<lower=0, upper=1> p; // probability of detection
}

model {

  // Priors
  lambda ~ normal(1, 2);
  pi ~ beta(1.5, 1.5);
  p ~ beta(1.5, 1.5);
  
  for (i in 1:N) {
    real lp_y_i = log(pi * (y[i] == 0)); // Contribution from zero inflation
    for (n in y[i]:min(100, y[i]+10)) { // Choose an upper bound for practical computation
      lp_y_i = log_sum_exp(lp_y_i,
                log(1 - pi) +
                poisson_lpmf(n | lambda) +
                binomial_lpmf(y[i] | n, p));
    }
    target += lp_y_i;
  }
}',
  data = data_list,
  init = init_list,
  chains = 4,
  iter = 2000,
  warmup = 1000,
  seed = 42
)

# Print results
print(fit)