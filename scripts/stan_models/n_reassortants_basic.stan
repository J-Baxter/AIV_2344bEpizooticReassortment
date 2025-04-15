//
// This Stan program defines a simple model, with a
// vector of values 'y' modeled as normally distributed
// with mean 'mu' and standard deviation 'sigma'.
//
// Learn more about model development with Stan at:
//
//    http://mc-stan.org/users/interfaces/rstan.html
//    https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started
//

data {
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
  // Abundance Model
  lambda ~ normal(3, 1.5);
  theta ~ beta(2, 5); 
  
  // Detection Model
  p ~ beta(2, 2); 
  
  // Loop over number of observations
  for (i in 1:N) {
    vector[K] lp;
    int c = continent[i]; // Current continent

    // Loop over plausible values of K to marginalize out discrete latent variables
    for (j in 1:K) {
      int current_population = y[i] + j - 1;
      
      // Likelihood for y[i] = 0
      
      if (y[i] == 0){
        vector[3] components;
        components[1] = bernoulli_lpmf(1 | theta);
        components[2] =  bernoulli_lpmf(0 | theta) + poisson_lpmf(current_population | lambda[c]);
        components[3] =  bernoulli_lpmf(0 | theta) + poisson_lpmf(current_population | lambda[c]) + binomial_lpmf(0 | current_population, p[c]);
        
        lp[j] = log_sum_exp(components);
        
        // Likelihood for y[i] > 0
        } else {
          lp[j] = bernoulli_lpmf(0 | theta) + poisson_lpmf(current_population | lambda[c]) + binomial_lpmf(y[i] | current_population, p[c]);
      }
      
      // Aggregate the probabilities 

    target += log_sum_exp(lp); 
    }
  }
}

