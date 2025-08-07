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

# 2. site-level covariates
model_cov_occ <- nimbleCode({
  
  p ~ dunif(0, 1)
  for (c in 1:ncov) {
    beta[c] ~ dunif(-100, 100)
  }
  
  # logit link - site-level covariates
  psi[1:nSites] <- ilogit(x_site[1:nSites, 1:ncov] %*% beta[1:ncov])
  
  for (i in 1:nSites) {
    
    z[i] ~ dbern(psi[i])
    
    for (j in 1:nVisits) {
      
      # observed data
      y[i, j] ~ dbern(z[i] * p)
      
    }
  }
  
})
