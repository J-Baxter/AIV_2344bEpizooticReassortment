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
  list(lambda = 1.0), 
  list(pi = 0.2),
  list(p = 0.5))

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
  pi ~ beta(1, 1);
  p ~ beta(1, 1);
  
  for (i in 1:N) {
    real lp_y_i = log(pi * (y[i] == 0)); // Contribution from zero inflation
    for (n in y[i]:100) { // Choose an upper bound for practical computation
      lp_y_i = log_sum_exp(lp_y_i,
                log(1 - pi) +
                poisson_lpmf(n | lambda) +
                binomial_lpmf(y[i] | n, p));
    }
    target += lp_y_i;
  }
}',
  data = data_list,
  #init = init_list,
  chains = 4,
  iter = 2000,
  warmup = 1000,
  seed = 42
)

# Print results

mod <- cmdstan_model(file)
fit <- mode$sample(data = data_list,
                   init = init_list,
                   chains = 2,
                   iter = 2000,
                   warmup = 1000,
                   seed = 42)
print(fit)

y <-
  structure(c(0L, 3L, 1L, 1L, 1L, 0L, 1L, 1L, 0L, 0L, 0L, 0L, 1L, 
              0L, 1L, 0L, 1L, 1L, 1L, 0L, 1L, 0L, 0L, 1L, 1L, 0L, 1L, 0L, 1L, 
              1L, 0L, 1L, 3L, 2L, 0L, 2L, 0L, 1L, 0L, 0L, 1L, 1L, 3L, 0L, 3L, 
              0L, 1L, 0L, 0L, 0L, 1L, 0L, 1L, 1L, 0L, 4L, 0L, 1L, 2L, 1L, 1L, 
              0L, 5L, 2L, 0L, 1L, 0L, 0L, 1L, 1L, 1L, 1L, 0L, 0L, 1L, 4L, 1L, 
              2L, 1L, 0L, 1L, 2L, 0L, 1L, 2L, 0L, 0L, 2L, 1L, 0L, 0L, 1L, 0L, 
              0L, 1L, 2L, 1L, 1L, 0L, 0L, 2L, 1L, 2L, 0L, 2L, 2L, 2L, 0L, 1L, 
              0L, 3L, 1L, 2L, 0L, 0L, 1L, 1L, 1L, 2L, 1L, 1L, 2L, 2L, 3L, 2L, 
              1L, 0L, 0L, 1L, 1L, 2L, 1L, 4L, 0L, 0L, 2L, 3L, 0L, 0L, 1L, 1L, 
              2L, 1L, 1L, 0L, 2L, 0L, 2L, 2L, 0L, 0L, 2L, 0L, 0L, 1L, 2L, 1L, 
              0L, 1L, 0L, 0L, 1L, 0L, 0L, 1L, 2L, 3L, 4L, 2L, 0L, 1L, 0L, 1L, 
              0L, 4L, 0L, 0L, 1L, 0L, 1L, 0L, 0L, 0L, 1L, 1L, 0L, 1L, 1L, 1L, 
              0L, 0L, 3L, 1L, 0L, 0L, 1L, 1L, 5L, 1L, 1L, 0L, 0L, 2L, 0L, 1L, 
              0L, 1L, 1L, 0L, 0L, 0L, 2L, 1L, 1L, 0L, 2L, 0L, 1L, 1L, 0L, 0L, 
              1L, 0L, 1L, 1L, 0L, 1L, 0L, 1L, 2L, 0L, 1L, 0L, 1L, 1L, 4L, 3L, 
              2L, 1L, 0L, 1L, 2L, 3L, 1L, 2L, 0L, 1L, 0L, 1L, 1L, 1L, 1L, 2L, 
              0L, 2L, 3L, 2L, 2L, 0L, 1L, 2L, 2L, 3L, 1L, 0L, 0L, 1L, 2L, 0L, 
              1L, 0L, 2L, 0L, 0L, 1L, 3L, 2L, 1L, 0L, 1L, 1L, 2L, 0L, 0L, 3L, 
              0L, 0L, 1L, 1L, 0L, 1L, 0L, 0L, 0L, 1L, 1L, 0L, 1L, 0L, 0L, 0L, 
              1L, 1L, 3L, 2L, 1L, 2L, 1L, 0L, 1L, 2L, 0L, 1L, 0L, 0L, 1L, 1L, 
              1L, 1L, 2L, 1L, 0L, 2L, 3L, 1L, 1L, 2L, 1L, 2L, 1L, 2L, 1L, 2L, 
              0L, 0L, 2L, 1L, 1L, 3L, 1L, 0L, 1L, 1L, 2L, 0L, 1L, 2L, 2L, 1L, 
              0L, 1L, 0L, 0L, 0L, 1L, 2L, 2L, 0L, 2L, 0L, 1L, 2L, 0L, 0L, 1L, 
              1L, 2L, 3L, 0L, 1L, 2L, 0L, 0L, 3L, 2L, 0L, 1L, 0L, 0L, 2L, 1L, 
              1L, 0L, 0L, 1L, 1L, 2L, 0L, 3L, 1L, 0L, 1L, 2L, 0L, 1L, 0L, 1L, 
              4L, 0L, 0L, 0L, 2L, 1L, 0L, 0L, 0L, 0L, 1L, 0L, 0L, 1L, 2L, 0L, 
              0L, 1L, 0L, 2L, 0L, 0L, 0L, 0L, 2L, 0L, 2L, 1L, 0L, 0L, 0L, 1L, 
              2L, 1L, 3L, 0L, 1L, 1L, 2L, 1L, 2L, 1L, 0L, 0L, 1L, 3L, 1L, 0L, 
              0L, 0L, 0L, 0L, 2L, 0L, 0L, 1L, 0L, 2L, 4L, 2L, 1L, 1L, 2L, 2L, 
              0L, 5L, 1L, 1L, 0L, 1L, 0L, 1L, 3L, 1L, 2L, 2L, 0L, 0L, 3L, 1L, 
              2L, 0L, 2L, 1L, 2L, 0L, 0L, 1L, 0L, 0L, 1L, 1L, 0L, 2L, 0L, 0L, 
              0L, 1L, 1L, 2L, 0L, 1L, 1L, 2L, 2L, 1L, 3L, 2L, 2L, 4L, 1L, 0L, 
              2L, 2L, 0L, 3L, 0L, 0L, 0L, 0L, 1L, 2L, 2L, 1L, 0L, 2L, 2L, 1L, 
              1L, 2L, 1L, 1L, 1L, 1L, 0L, 2L, 0L, 0L, 1L, 2L, 1L, 2L, 2L, 1L, 
              3L, 1L, 2L, 0L, 0L, 1L, 2L, 0L, 0L, 1L, 1L, 0L, 0L, 0L, 1L, 0L, 
              0L, 1L, 1L, 1L, 1L, 0L, 0L, 1L, 1L, 3L, 3L, 0L, 0L, 1L, 1L, 0L, 
              0L, 3L, 1L, 0L, 0L, 0L, 1L, 1L, 1L, 0L, 1L, 2L, 0L, 0L, 1L, 1L, 
              1L, 1L, 2L, 1L, 0L, 1L, 1L, 1L, 3L, 0L, 1L), .Dim = c(200L, 3L
              ))
R <-
  200L
T <-
  3L
K <-
  100L
## Parameters monitored
params <- c("lambda", "p")

## MCMC settings
ni <- 1000
nt <- 1
nb <- 500
nc <- 4

## Initial values
inits <- lapply(1:nc, function(i)
  list(p = runif(1, 0, 1),
       lambda = runif(1, 0, 1)))

out <- stan("~/Downloads/binmix.stan",
            data = list(y = y, # Data
                        R = R, # Number of sites
                        T = T, # Temporal replications
                        K = 100), # Upper bound of population size (for marginalisation)
            init = inits, pars = params,
            chains = nc, iter = ni, warmup = nb, thin = nt,
            seed = 1,
            open_progress = FALSE)