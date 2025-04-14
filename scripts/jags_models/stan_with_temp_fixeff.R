stan_code <- 'data {
  int<lower=0> N; // Number of observations
  array[N] int<lower=0> y; // Observed counts
  int<lower=0> K; // Upper bound of reassortants
  int<lower=1> C; // Number of continents
  array[N] int<lower=1, upper=C> continent; // Continent index for each observation
  array[N] real cases; // Additional data for cases
  array[N] real sequences; // Additional data for sequences
}

parameters {
  real<lower=0, upper=1> theta; // Zero-inflation probability
  array[C] real<lower=0> continent_specific_abundance; // Poisson rate stratified by continent
  array[C] real<lower=0, upper=1> continent_specific_detection; // Detection probability stratified by continent
  real beta_cases; // Coefficient for additional cases data
  real beta_sequences; // Coefficient for additional sequences data
}

model {
  // Priors
  // Abundance Model
  continent_specific_abundance ~ normal(3, 1.5);
  beta_cases ~ normal(0, 1);
  
  // Detection model
  theta ~ beta(2, 5); // Example prior
  beta_sequences ~ normal(0, 1);
  
  // Zero Inflation Model
  continent_specific_detection ~ beta(2, 2); // Example prior
  
  // Loop over number of observations
  for (i in 1:N) {
    vector[K] lp;
    int c = continent[i]; // Current continent
    
    // Linear predictors 
    real lambda = exp(continent_specific_abundance[c] + beta_cases * cases[i]);
    real p = inv_logit(continent_specific_detection[c] + beta_sequences * sequences[i]);

    // Loop over plausible values of K to marginalise out discrete latent variables
    for (j in 1:K) {
    
      int current_population = y[i] + j - 1;
      
      // Likelihood for y[i] = 0
      if (y[i] == 0){
      
      vector[2] components;
      components[1] = bernoulli_lpmf(1 | theta);
      components[2] = bernoulli_lpmf(0 | theta) +
                      poisson_lpmf(current_population | lambda) +
                      binomial_lpmf(0 | current_population, p); //CHECK THIS LOGIC - YOU HAVE REMOVED SOMETHING
                      
      lp[j] = log_sum_exp(components);
      
      // Likelihood for y[i] > 0
      } else {
      
      lp[j] = bernoulli_lpmf(0 | theta) +
              poisson_lpmf(current_population | lambda) +
              binomial_lpmf(y[i] | current_population, p);
      }
    
     // Aggregate the probabilities 
    target += log_sum_exp(lp);
  }
}'