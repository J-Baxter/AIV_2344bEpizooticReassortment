
// This Stan program defines a simple model, with a
// vector of values 'y' modeled as normally distributed
// with mean 'mu' and standard deviation 'sigma'.
//
// Learn more about model development with Stan at:
//
//    http://mc-stan.org/users/interfaces/rstan.html
//    https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started
//

// The input data 
data {
  int<lower=0> N; // Number of observations
  array[N] int<lower=0> y; // Observed counts
  int<lower=0> K; // Upper bound of reassortants
  int<lower=1> C; // Number of continents
  int<lower=1> Y; // Number of years
  array[N] int<lower=1, upper=C> continent_index; // Continent index for each observation
  array[N] int<lower=1, upper=Y> year_index; // Year index for each observation
  array[N] real cases; // Additional data for cases
  array[N] real sequences; // Additional data for sequences
}


// Declared parameters
parameters {
  array[C] real<lower=0, upper=1> continent_specific_theta; // Zero-inflation probability stratified by continent
  array[C] real<lower=0> continent_specific_abundance; // Poisson rate stratified by continent
  array[C] real<lower=0, upper=1> continent_specific_detection; // Detection probability stratified by continent
  
  // 'fixed' effects
  real beta_cases; // Coefficient for additional cases data
  real beta_sequences; // Coefficient for additional sequences data
  
  // 'random' effects
  matrix[C, Y] z_abundance; // Standard normal for abundance random effects
  matrix[C, Y] z_detection; // Standard normal for detection random effects
  //matrix[2, Y] z_theta; // Standard normal for theta random effects

  cholesky_factor_corr[C] L_Omega_abundance;
  vector<lower=0>[C] sigma_abundance;

  cholesky_factor_corr[C] L_Omega_detection;
  vector<lower=0>[C] sigma_detection;

  //cholesky_factor_corr[2] L_Omega_theta;
 // vector<lower=0>[2] sigma_theta;
}


transformed parameters {
  matrix[C, Y] abundance_re;
  matrix[C, Y] detection_re;
  //matrix[C, Y] theta_re;

  abundance_re = diag_pre_multiply(sigma_abundance, L_Omega_abundance) * z_abundance;
  detection_re = diag_pre_multiply(sigma_detection, L_Omega_detection) * z_detection;
  //theta_re = (diag_pre_multiply(sigma_theta, L_Omega_theta) * z_theta)';
}


// The model to be estimated
model {
  // Priors
  // Abundance Model
  continent_specific_abundance ~ normal(3, 1.5);
  beta_cases ~ normal(0, 1);
  L_Omega_abundance ~ lkj_corr_cholesky(1);
  sigma_abundance ~ normal(0, 1);
  to_vector(z_abundance) ~ normal(0, 1);
 
  
  // Detection model
  continent_specific_detection ~ beta(1,1);//beta(1.5, 1.5);
  beta_sequences ~ normal(0, 1);
  L_Omega_detection ~ lkj_corr_cholesky(1);
  sigma_detection ~ normal(0, 1);
  to_vector(z_detection) ~ normal(0, 1);
  

  // Zero Inflation Model
  continent_specific_theta ~ beta(2, 5); 
  //L_Omega_theta ~ lkj_corr_cholesky(1);
  //sigma_theta ~ normal(0, 1);
  //to_vector(z_theta) ~ normal(0, 1);
  

  // Loop over number of observations
  for (i in 1:N) {
    vector[K] lp;
    int c = continent_index[i]; // Current continent
    int yr = year_index[i]; //Current year

    
    // Linear predictors 
    real lambda = exp(continent_specific_abundance[c] + beta_cases * cases[i] + abundance_re[c, yr] ); 
    real p = inv_logit(continent_specific_detection[c] + beta_sequences * sequences[i]); 

    // Loop over plausible values of K to marginalise out discrete latent variables
    for (j in 1:K) {
    
      int current_population = y[i] + j - 1;
      
      // Likelihood for y[i] = 0
      if (y[i] == 0){
      
      vector[3] components;
      components[1] = bernoulli_lpmf(1 | continent_specific_theta[c]);
      components[2] = bernoulli_lpmf(0 | continent_specific_theta[c]) + poisson_lpmf(current_population | lambda);
      components[3] = bernoulli_lpmf(0 | continent_specific_theta[c]) + poisson_lpmf(current_population | lambda) + binomial_lpmf(0 | current_population, p);
                      
      lp[j] = log_sum_exp(components);
      
      // Likelihood for y[i] > 0
      } else {
      
      lp[j] = bernoulli_lpmf(0 | continent_specific_theta[c]) + poisson_lpmf(current_population | lambda) + binomial_lpmf(y[i] | current_population, p);
      }
    }
    // Aggregate the probabilities 
    target += log_sum_exp(lp);
  }
}


// Replications for the posterior predictive distribution
generated quantities {
  array[N] int y_rep; 

  for (i in 1:N) {
    int c = continent_index[i]; // Current continent
    int yr = year_index[i]; //Current year
    
    // Draw a latent count from the zero-inflated Poisson
    int current_population;
    if (bernoulli_rng(continent_specific_theta[c]) == 1) {
      // Zero inflation: y_rep[i] = 0
      y_rep[i] = 0;
    } else {
      // Simulate from Poisson
      current_population = poisson_rng(exp(continent_specific_abundance[c] + beta_cases * cases[i] + abundance_re[c, yr] ));
      
      // Simulate observed count from a binomial
      y_rep[i] = binomial_rng(current_population, inv_logit(continent_specific_detection[c] + beta_sequences * sequences[i] + detection_re[c, yr] ));
    }
  }
}


