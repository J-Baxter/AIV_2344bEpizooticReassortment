model <- "
model {
  
  # Likelihood for the True Reassortant Count (Zero-Inflated Poisson)
  for (i in 1:N) {
    # True count of reassortants follows a Poisson distribution
    true_reassortants[i] ~ dpois(lambda[i])
    
    # Lambda (mean) of the Poisson distribution depends on covariates
    log(lambda[i]) <- beta_regionname[collection_regionname[i]] * n_sequences[i] + 
                      beta_reassortants[collection_regionname[i]] * n_reassortants[i] + 
                      beta_collection[collection_regionname[i]] + 
                      season_random[collection_regionname[i], collection_year[i], collection_season[i]] + 
                      year_random[collection_regionname[i], collection_year[i]]
    
    # Zero-Inflation Model for True Reassortants (Bernoulli for zero vs non-zero)
    zi[i] ~ dbern(psi[i])  # Zero-Inflation probability
    logit(psi[i]) <- zi_intercept + 
                     region_random[collection_regionname[i]] + 
                     year_season_random[collection_regionname[i], collection_year[i], collection_season[i]]
  }
  
  # Likelihood for Observed Reassortants (Binomial / Bernoulli)
  for (i in 1:N) {
    # The observed number of reassortants follows a binomial distribution
    observed_reassortants[i] ~ dbin(p[i], true_reassortants[i])
    logit(p[i]) <- beta_detection[collection_regionname[i]] + 
                   region_random[collection_regionname[i]] + 
                   season_random[collection_regionname[i], collection_year[i], collection_season[i]]
  }
  
  # Model for the number of reassortants (n_reassortants)
  for (i in 1:N) {
    n_reassortants[i] ~ dpois(lambda_reassortants[i])
    log(lambda_reassortants[i]) <- beta_regionname_woah[collection_regionname[i]] * woah_susceptibles[i] + 
                                   region_random[collection_regionname[i]] + 
                                   year_season_random[collection_regionname[i], collection_year[i], collection_season[i]]
  }
  
  # Priors for Fixed Effects (Beta Coefficients)
  for (r in 1:R) {
    beta_regionname[r] ~ dnorm(0, 0.001)            # Fixed effect for region and n_sequences
    beta_reassortants[r] ~ dnorm(0, 0.001)           # Fixed effect for region and n_reassortants
    beta_collection[r] ~ dnorm(0, 0.001)            # Fixed effect for collection
    beta_detection[r] ~ dnorm(0, 0.001)             # For the detection probability of observed reassortants
    beta_regionname_woah[r] ~ dnorm(0, 0.001)        # Fixed effect for region and woah_susceptibles
  }
  
  # Priors for Random Effects (Hierarchical Structure)
  for (c in 1:C) {
    region_random[c] ~ dnorm(0, tau_region)  # Random effect for region
  }
  
  for (y in 1:Y) {
    year_random[y] ~ dnorm(0, tau_year)  # Random effect for collection_year within region
  }
  
  for (c in 1:C) {
    for (y in 1:Y) {
      for (s in 1:S) {
        season_random[c, y, s] ~ dnorm(0, tau_season)  # Random effect for collection_season within collection_year and region
      }
    }
  }
  
  # Priors for Variances (Precision of Random Effects)
  tau_region ~ dgamma(0.001, 0.001)
  tau_year ~ dgamma(0.001, 0.001)
  tau_season ~ dgamma(0.001, 0.001)
  
  # Prior for Zero-Inflation Intercept
  zi_intercept ~ dnorm(0, 0.001)
}

"


dummy_data <- list(
  N = 100,
  M = 100,
  K = 100,
  collection_regionname = sample(1:5, 100, replace = TRUE),
  collection_year = sample(2015:2020, 100, replace = TRUE),
  collection_season = sample(1:4, 100, replace = TRUE),
  n_sequences = rpois(100, lambda = 20),
  n_reassortants = rpois(100, lambda = 5),
  woah_susceptibles = sample(1:3, 100, replace = TRUE),
  reassortants_observed = rpois(100, lambda = 10),
  zi = rbinom(100, size = 1, prob = 0.1)
)


inits <- function() {
  list(
    # Fixed Effects Initial Values
    beta_regionname = rep(0, 5),           # Initial values for each region's fixed effect
    beta_reassortants = rep(0, 5),          # Initial values for reassortants' effect
    beta_collection = rep(0, 5),           # Initial values for collection effect
    
    # Random Effects Initial Values
    region_random = rep(0, 5),             # Initial values for regional random effects
    season_random = rep(0, 4),             # Initial values for seasonal random effects
    
    # Hyperparameters (Variances)
    tau_region = 1,                        # Initial value for the variance of region effects
    tau_season = 1,                        # Initial value for the variance of season effects
    
    # Zero-Inflation (zi) Initial Values
    zi_intercept = 0                        # Initial value for zero-inflation intercept
  )
}


te_jags <- jags.model(textConnection(model), 
                        data = dummy_data,
                        n.chains = 2,
                        n.adapt = 5000, 
                        inits = inits)


M <- 150
J <- 2
C <- matrix(NA, nrow = M, ncol = J)

lambda <- 2.5
p <- 0.4
psi <- 0.2

N <- rpois(n = M, lambda = lambda) * rbinom(M, 1, 1-psi)

for (j in 1:J){
  C[,j] <- rbinom(n = M, size = N, prob = p)
}

table(N)
sum(N)
sum(N>0)
mean(N)

table(apply(C, 1, max))
sum(apply(C, 1, max))
sum(apply(C, 1, max)>0)
mean(apply(C, 1, max))

# COMPILE the model 

vote_jags <- jags.model(textConnection(model), 
                        data = list(C = C, M = nrow(C), J = ncol(C)),
                        n.chains = 2,
                        n.adapt = 5000, 
                        inits = function(){
                          list(lambda = runif(1, 0, 10),
                               psi = runif(1, 0, 1), 
                               p = runif(1, 0, 1),
                               N = apply(C, 1, max), 
                               z = rep(1, nrow(C)))
                        })

# SIMULATE the posterior
vote_sim <- coda.samples(model = vote_jags, 
                         variable.names = c("lambda", "p", 'psi'), 
                         n.iter = 25000,
                         thin = 20)

# PLOT the posterior
summary(vote_sim)
