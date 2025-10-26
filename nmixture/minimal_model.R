# 1. basic n-mixture model
model_minimal_nmixture <- nimbleCode({
  
  mu ~ dunif(0, 10000)
  p ~ dbeta(1, 1)
  
  for (i in 1:nSites) {
    
    N[i] ~ dpois(mu / p) 
    
    for (j in 1:nVisits) {
      
      # observed data
      y[i, j] ~ dbinom(size = N[i], prob = p)
      
    }
  }
})


# 1. basic n-mixture model - log scale
model_minimal_nmixture_logmu <- nimbleCode({
  
  log_mu ~ dunif(-100, 100)
  p ~ dbeta(1, 1)
  
  for (i in 1:nSites) {
    
    N[i] ~ dpois(exp(log_mu) / p) 
    
    for (j in 1:nVisits) {
      
      # observed data
      y[i, j] ~ dbinom(size = N[i], prob = p)
      
    }
  }
})

# 1. basic n-mixture model - log scale
model_minimal_nmixture_logp <- nimbleCode({
  
  mu ~ dunif(0, 10000)
  log_p ~ dunif(-100, 100)
  
  p <- 1 / (1 + exp(-log_p))
  
  for (i in 1:nSites) {
    
    N[i] ~ dpois(mu / p) 
    
    for (j in 1:nVisits) {
      
      # observed data
      y[i, j] ~ dbinom(size = N[i], prob = p)
      
    }
  }
})

# nimbleecology n-mixture
model_minimal_nmixture_NE <- nimbleCode({
  
  lambda ~ dunif(0, 10000)
  p ~ dunif(0, 1)
  
  for (i in 1:nSites) {
    
    y[i, 1:nVisits] ~ dNmixture_s(lambda = lambda, p = p, 
                                  Nmin = -1, Nmax = -1, len = nVisits)
  }
})



