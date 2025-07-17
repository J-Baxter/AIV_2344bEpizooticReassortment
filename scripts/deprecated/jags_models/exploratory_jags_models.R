# constant expected abundance and detection 
# lambda = 2.5
# p = 0.4
# sample size = 150 sites
# repeated measurements = 2

M <- 150
J <- 2
C <- matrix(NA, nrow = M, ncol = J)

lambda <- 2.5
p <- 0.4

N <- rpois(n = M, lambda = lambda)

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

# JAGS model



# DEFINE the model
model <- "model {
  # Priors
  lambda ~ dgamma(0.001, 0.001)   # Prior for mean abundance
  p ~ dunif(0, 1)          # Prior for detection probability
  
  # Loop over sites
  for (i in 1:M) {
    N[i] ~ dpois(lambda)   # True abundance at site i
    
    # Loop over surveys
    for (j in 1:J) {
      C[i, j] ~ dbin(p, N[i])  # Observed count given detection probability
    }
  }
}
"

# COMPILE the model 

vote_jags <- jags.model(textConnection(model), 
                        data = list(C = C, M = nrow(C), J = ncol(C)),
                        n.chains = 2,
                        n.adapt = 5000,
                        inits = function(){list(N = apply(C, 1, max) )})

# SIMULATE the posterior
vote_sim <- coda.samples(model = vote_jags, 
                         variable.names = c("lambda", "p"), 
                         n.iter = 25000,
                         thin = 20)

# PLOT the posterior
summary(vote_sim)
plot(vote_sim, trace = FALSE)


###################################################################################################
# Basic model (complete pooling across sites + time)
# DEFINE the model
# sample size 

data_processed %>%
  filter(collection_regionname != 'central & northern america') %>%
  select(collection_regionname, n_reassortants) %>%
  group_split(collection_regionname, .keep =F) %>%
  lapply(., function(x) as.matrix(x) %>% t()) %>%
  do.call(rbind, .)


# Repeat Measures


C <- data_processed %>%
  filter(collection_regionname != 'central & northern america') %>%
  select(collection_regionname, n_reassortants) %>%
  group_split(collection_regionname, .keep =F) %>%
  lapply(., function(x) as.matrix(x) %>% t()) %>%
  do.call(rbind, .)

n_continents <- nrow(C)
n_observations <- ncol(C)

model <- "model {
  # Priors
  lambda ~ dgamma(0.001, 0.001)   # Prior for mean abundance
  p ~ dunif(0, 1)          # Prior for detection probability
  
  # Loop over sites
  for (i in 1:continents) {
    N[i] ~ dpois(lambda)   # True abundance at site i
    
    # Loop over surveys
    for (j in 1:observations) {
      C[i, j] ~ dbin(p, N[i])  # Observed count given detection probability
    }
  }
}
"

# COMPILE the model 

vote_jags <- jags.model(textConnection(model), 
                        data = list(C = C, continents = nrow(C), observations = ncol(C)),
                        n.chains = 2,
                        n.adapt = 5000,
                        inits = function(){list(N = apply(C, 1, max) )})

# SIMULATE the posterior
vote_sim_2 <- coda.samples(model = vote_jags, 
                         variable.names = c("lambda", "p"), 
                         n.iter = 25000,
                         thin = 20)

# PLOT the posterior
summary(vote_sim_2)

####################################################################################################
# Basic model with fixed effect example
model <- "model {
  # Priors
  alpha ~ dnorm(0, 0.001)  # Intercept

  # Priors for categorical predictor: collection_regionname (factor levels)
  for (r in 1:nRegions) {
    beta_region[r] ~ dnorm(0, 0.1)
  }

  # Continuous predictor: woah_susceptibles_log1p
  beta_susceptibles ~ dnorm(0, 0.1)

  # Priors for categorical predictor: collection_season (factor levels)
  for (s in 1:nSeasons) {
    beta_season[s] ~ dnorm(0, 0.1)
  }

  # Likelihood
  for (i in 1:N) {
    log(lambda[i]) <- alpha + 
                      beta_region[collection_regionname[i]] + 
                      beta_susceptibles * woah_susceptibles_log1p[i] + 
                      beta_season[collection_season[i]]
    
    n_reassortants[i] ~ dpois(lambda[i])  # Poisson likelihood
  }
}
"
# COMPILE the model 

vote_jags <- jags.model(textConnection(model), 
                        data = list(C = C, continents = nrow(C), observations = ncol(C)),
                        n.chains = 2,
                        n.adapt = 5000,
                        inits = function(){list(N = apply(C, 1, max) )})

# SIMULATE the posterior
vote_sim <- coda.samples(model = vote_jags, 
                         variable.names = c("lambda", "p"), 
                         n.iter = 25000,
                         thin = 20)

# PLOT the posterior
summary(vote_sim)
