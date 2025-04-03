stan_code <- "data {
  int<lower=0> N; // Number of observations
  array[N] int<lower=0> y; // Observed counts
  int<lower=0> K; // Upper bound of population size
  int<lower=1> C; // Number of continents
  int<lower=1> Y; // Number of years
  array[N] int<lower=1, upper=C> continent; // Continent index for each observation
  array[N] int<lower=1, upper=Y> year_index; // Year index for each observation
}

parameters {
  vector<lower=0>[C] lambda; // Fixed effects for Poisson rate by continent
  vector<lower=0, upper=1>[C] p; // Fixed effects for detection probability by continent
  real<lower=0, upper=1> theta; // Fixed zero-inflation probability

  // Random effects
  matrix[2, Y] z_lambda; // Standard normal for lambda random effects
  matrix[2, Y] z_p; // Standard normal for p random effects
  matrix[2, Y] z_theta; // Standard normal for theta random effects

  cholesky_factor_corr[2] L_Omega_lambda;
  vector<lower=0>[2] sigma_lambda;

  cholesky_factor_corr[2] L_Omega_p;
  vector<lower=0>[2] sigma_p;

  cholesky_factor_corr[2] L_Omega_theta;
  vector<lower=0>[2] sigma_theta;
}

transformed parameters {
  matrix[C, Y] lambda_re;
  matrix[C, Y] p_re;
  matrix[C, Y] theta_re;

  lambda_re = (diag_pre_multiply(sigma_lambda, L_Omega_lambda) * z_lambda)';
  p_re = (diag_pre_multiply(sigma_p, L_Omega_p) * z_p)';
  theta_re = (diag_pre_multiply(sigma_theta, L_Omega_theta) * z_theta)';
}

model {
  // Priors
  lambda ~ normal(3, 1);
  p ~ beta(2, 5);
  theta ~ beta(2, 5);

  L_Omega_lambda ~ lkj_corr_cholesky(1);
  sigma_lambda ~ normal(0, 1);
  to_vector(z_lambda) ~ normal(0, 1);

  L_Omega_p ~ lkj_corr_cholesky(1);
  sigma_p ~ normal(0, 1);
  to_vector(z_p) ~ normal(0, 1);

  L_Omega_theta ~ lkj_corr_cholesky(1);
  sigma_theta ~ normal(0, 1);
  to_vector(z_theta) ~ normal(0, 1);

  // Loop over number of observations
  for (i in 1:N) {
    vector[K] lp;
    int c = continent[i];
    int yr = year_index[i];

    real lambda_obs = exp(log(lambda[c]) + lambda_re[c, yr]);
    real p_obs = inv_logit(logit(p[c]) + p_re[c, yr]);
    real theta_obs = inv_logit(logit(theta) + theta_re[c, yr]);

    // Loop over plausible values of K to marginalize out discrete latent variables
    for (j in 1:K) {
      int current_population = y[i] + j - 1;

      if (y[i] == 0)
        lp[j] = log_sum_exp(to_vector({
          bernoulli_lpmf(1 | theta_obs), 
          bernoulli_lpmf(0 | theta_obs) + poisson_lpmf(current_population | lambda_obs),
          bernoulli_lpmf(0 | theta_obs) + poisson_lpmf(current_population | lambda_obs) + binomial_lpmf(y[i] | current_population, p_obs)
        }));
      else
        lp[j] = bernoulli_lpmf(0 | theta_obs) + poisson_lpmf(current_population | lambda_obs) + binomial_lpmf(y[i] | current_population, p_obs);
    }

    target += log_sum_exp(lp);
  }
}

"

mod <- cmdstan_model(write_stan_file(stan_code))

set.seed(123)

# Parameters
N <- 100  # Number of observations
K <- 20   # Upper bound of population size
C <- 3    # Number of continents
Y <- 5    # Number of years

# True parameters
true_lambda <- c(2.5, 3.0, 3.5)
true_p <- c(0.4, 0.6, 0.5)
true_theta <- 0.2

# Random effects
sigma_lambda <- 0.5
sigma_p <- 0.3
sigma_theta <- 0.2

# Generate data
continent <- sample(1:C, N, replace = TRUE)
year_index <- sample(1:Y, N, replace = TRUE)

lambda_re <- matrix(rnorm(C * Y, 0, sigma_lambda), nrow = C, ncol = Y)
p_re <- matrix(rnorm(C * Y, 0, sigma_p), nrow = C, ncol = Y)
theta_re <- matrix(rnorm(C * Y, 0, sigma_theta), nrow = C, ncol = Y)

y <- integer(N)

for (i in 1:N) {
  c <- continent[i]
  yr <- year_index[i]
  
  lambda <- exp(log(true_lambda[c]) + lambda_re[c, yr])
  p <- plogis(qlogis(true_p[c]) + p_re[c, yr])
  theta <- plogis(qlogis(true_theta) + theta_re[c, yr])
  
  if (runif(1) < theta) {
    y[i] <- 0
  } else {
    population_size <- rpois(1, lambda)
    y[i] <- rbinom(1, population_size, p)
  }
}

# Print the simulated data
data_list <- list(N = N, K = K, C = C, Y = Y, y = y, continent = continent, year_index = year_index)
print(data_list)



# Run the model
fit <- mod$sample(
  data = data_list,
  seed = 42,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 500,
  iter_sampling = 2000
)

print(fit$summary())

