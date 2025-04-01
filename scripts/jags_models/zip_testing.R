model <- "
model {
  # Specify priors
  # zero-inflation/suitability
  phi ~ dbeta(1, 1)            # Proportion of suitable sites
  theta <-  1-phi                # Zero inflation (ie proportion UNsuitable)
  ltheta <- logit(theta)
  
  # abundance
    lambda ~ dunif(0, 100)       # Mean abundance per site
  
  p ~ dbeta(1, 1)              # Detection probability

  for (i in 1:M) {
    # True abundance (N) follows a Zero-Inflated Poisson (ZIP)
    z[i] ~ dbern(1 - psi)      # Latent indicator: 1 if site is occupied
    N[i] ~ dpois(lambda * z[i])# Zero-inflated Poisson process

    for (j in 1:J) {
      # Observation model: Binomial sampling given true abundance
      C[i, j] ~ dbin(p, N[i])
    }
  }
}
"
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
