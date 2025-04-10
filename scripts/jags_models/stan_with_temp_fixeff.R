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
    
    // Linear predictors - requires updating
    real lambda = exp(alpha_abundance + X_abundance[j] * beta_abundance);
    real p = inv_logit(alpha_detection + X_detection[j] * beta_detection);

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