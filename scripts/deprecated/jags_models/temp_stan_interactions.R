'data {
  int<lower=0> N; // Number of observations
  int<lower=0> K; // Upper bound of population size
  int<lower=1> C; // Number of continents
  array[N] int<lower=1, upper=C> continent; // Continent index for each observation
  array[N] int<lower=0> y; // Observed counts
  int<lower=1> R; // Number of collection regions
  array[N] int<lower=1, upper=R> collection_regionname; // Collection region index
  array[N] real<lower=0> n_sequences; // Number of sequences
  array[N] real<lower=0> n_reassortants; // Number of reassortants
  array[N] real<lower=0> woah_susceptibles; // Susceptibles measure
}

parameters {
  real<lower=0, upper=1> theta; // Zero-inflation probability
  array[C] real<lower=0> lambda; // Poisson rate stratified by continent
  array<C] real<lower=0, upper=1> p; // Detection probability stratified by continent
  array[R] real alpha; // Intercept for collection regions
  real beta_woah; // Coefficient for woah_susceptibles
  real beta_nseq; // Coefficient for n_sequences
  real beta_nreas; // Coefficient for n_reassortants
}

model {
  // Priors
  lambda ~ normal(3, 1.5);
  theta ~ beta(2, 5);
  p ~ beta(2, 2);
  alpha ~ normal(0, 1); // Prior for region intercepts
  beta_woah ~ normal(0, 1); // Prior for woah_susceptibles
  beta_nseq ~ normal(0, 1); // Prior for n_sequences
  beta_nreas ~ normal(0, 1); // Prior for n_reassortants
  
  // Loop over observations
  for (i in 1:N) {
    vector[K] lp;
    int c = continent[i]; // Current continent
    int r = collection_regionname[i]; // Current region

    // Linear predictor
    real eta = alpha[r] + beta_woah * woah_susceptibles[i] + 
               beta_nseq * n_sequences[i] + beta_nreas * n_reassortants[i];

    // Loop over plausible values of K
    for (j in 1:K) {
      int current_population = y[i] + j - 1;
      
      if (y[i] == 0)
        lp[j] = log_sum_exp(to_vector({
          bernoulli_lpmf(1 | theta), 
          bernoulli_lpmf(0 | theta) + poisson_lpmf(current_population | lambda[c]),
          bernoulli_lpmf(0 | theta) + poisson_lpmf(current_population | lambda[c]) + binomial_lpmf(y[i] | current_population, p[c])
        })) + eta;
      
      else
        lp[j] = bernoulli_lpmf(0 | theta) + poisson_lpmf(current_population | lambda[c]) + binomial_lpmf(y[i] | current_population, p[c]) + eta;
    }

    target += log_sum_exp(lp); // Aggregate the probabilities correctly
  }
}
'