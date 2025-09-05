# 1. basic n-mixture model
model_minimal_nmixture <- nimbleCode({
  
  lambda ~ dunif(0, 1000)
  p ~ dbeta(1, 1)
  
  for (i in 1:nSites) {
    
    N[i] ~ dpois(lambda) 
    
    for (j in 1:nVisits) {
      
      # observed data
      y[i, j] ~ dbinom(size = N[i], prob = p)
      
    }
  }
})

