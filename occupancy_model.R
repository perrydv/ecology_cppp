# 1. basic version
# 2-3. covariates on detection and occupancy
# 4. spatial random effect
# 5. multi-species?


# 1. basic occupancy model
model_basic <- nimbleCode({
  
  psi ~ beta(1, 1)
  p ~ beta(1, 1)
  
  for (i in 1:nSites) {
    
    z[i] ~ dbern(psi)
    
    for (j in 1:nVisits) {
      
      # observed data
      y[i, j] ~ dbern(z[i] * p)
      
      # log deviance for observed data
      D_obs[i, j] <- -2 * log(dbern(y[i, j], z[i] * p) + 1e-6)  
      
      # posterior predictive replicated data
      y_rep[i, j] ~ dbern(z[i] * p)
      D_rep[i, j] <- -2 * log(dbern(y_rep[i, j], z[i] * p) + 1e-6)
      
    }
  }
  
  # total discrepancy measures
  D_obs_total <- sum(D_obs[1:nSites, 1:nVisits])
  D_rep_total <- sum(D_rep[1:nSites, 1:nVisits])
  
  
})

# 2. covariates on occupancy
model_cov_occ <- nimbleCode({
  
  p ~ dunif(0, 1)
  for (c in 1:ncov) {
    beta[c] ~ dunif(-10, 10)
  }
  
  for (i in 1:nSites) {
    
    z[i] ~ dbern(psi[i])
    
    # logit link - site-level covariates
    psi[i] <- ilogit(beta[1:ncov] %*% x_site[i, 1:ncov])
    
    for (j in 1:nVisits) {
      
      # observed data
      y[i, j] ~ dbern(z[i] * p)
      
      # log deviance for observed data
      D_obs[i, j] <- -2 * log(dbern(y[i, j], z[i] * p) + 1e-6) 
      
      # posterior predictive replicated data
      y_rep[i, j] ~ dbern(z[i] * p)
      D_rep[i, j] <- -2 * log(dbern(y_rep[i, j], z[i] * p) + 1e-6)
      
    }
  }
  
  # total discrepancy measures
  D_obs_total <- sum(D_obs[1:nSites, 1:nVisits])
  D_rep_total <- sum(D_rep[1:nSites, 1:nVisits])
  
})

# 3. covariates on occupancy and detection
model_cov_occ_det <- nimbleCode({
  
  p ~ dunif(0, 1)
  for (c in 1:ncov_o) {
    beta_o[c] ~ dunif(-10, 10)
  }
  for (c in 1:ncov_p) {
    beta_p[c] ~ dunif(-10, 10)
  }
  
  for (i in 1:nSites) {
    
    z[i] ~ dbern(psi[i])
    
    # logit link - site-level covariates
    psi[i] <- ilogit(beta_o[1:ncov_o] %*% x_site[i, 1:ncov_o])
    
    for (j in 1:nVisits) {
      
      # logit link - detection-level covariates
      p[i, j] <- ilogit(beta_o[1:ncov_p] %*% x_visit[i, j, 1:ncov_p])
      
      # observed data
      y[i, j] ~ dbern(z[i] * p[i, j])
      
      # log deviance for observed data
      D_obs[i, j] <- -2 * log(dbern(y[i, j], z[i] * p[i, j]) + 1e-6) 
      
      # posterior predictive replicated data
      y_rep[i, j] ~ dbern(z[i] * p[i, j])
      D_rep[i, j] <- -2 * log(dbern(y_rep[i, j], z[i] * p[i, j]) + 1e-6)
      
      
    }
  }
  
  # total discrepancy measures
  D_obs_total <- sum(D_obs[1:nSites, 1:nVisits])
  D_rep_total <- sum(D_rep[1:nSites, 1:nVisits])
  
})

# 4. spatial random effect
model_spatial_ranef <- nimbleCode({
  
  psi_mean ~ dunif(-10, 10)
  psi_sd ~ dunif(0, 100)
  p ~ dbeta(1, 1)
  
  for(r in 1:nRegions){
    
    psi_logit[r] ~ dnorm(psi_mean, psi_sd)
    psi[r] <- ilogit(psi_logit[r])
    
  }
  
  for(i in 1:nSites) {
    
    z[i] ~ dbern(psi[region[i]])
    
    for(j in 1:nVisits) {
      
      # observed data
      y[i, j] ~ dbern(z[i] * p)
      
      # log deviance for observed data
      D_obs[i, j] <- -2 * log(dbern(y[i, j], z[i] * p) + 1e-6)  
      
      # posterior predictive replicated data
      y_rep[i, j] ~ dbern(z[i] * p)
      D_rep[i, j] <- -2 * log(dbern(y_rep[i, j], z[i] * p) + 1e-6)
      
    }
  }
  
  # total discrepancy measures
  D_obs_total <- sum(D_obs[1:nSites, 1:nVisits])
  D_rep_total <- sum(D_rep[1:nSites, 1:nVisits])
  
})

# 5. multi-species (https://doi.org/10.1111/2041-210X.12587)
# model_multispec <- nimbleCode({
  
# })

