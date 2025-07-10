# 1. basic occupancy model
model_minimal <- nimbleCode({
  
  psi ~ dbeta(1, 1)
  p ~ dbeta(1, 1)
  
  for (i in 1:nSites) {
    
    z[i] ~ dbern(psi) 
    
    # observed data
    y[i] ~ dbinom(z[i] * p, nVisits)
    
  }
  
})
