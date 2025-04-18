
####################
# model variations #
####################


# 1. basic version
# 2-3. covariates on detection and occupancy
# 4. spatial random effect
# 5. basic (binomial detection)


# 1. basic occupancy model
model_basic <- nimbleCode({
  
  psi ~ dbeta(1, 1)
  p ~ dbeta(1, 1)
  
  for (i in 1:nSites) {
    
    z[i] ~ dbern(psi) 
    
    # posterior predictive replicated data
    z_rep[i] ~ dbern(psi)
    
    for (j in 1:nVisits) {
      
      # observed data
      y[i, j] ~ dbern(z[i] * p)
      
      # log deviance - observed data
      D_obs[i, j] <- -2 * log(dbinom(y[i, j], 1, z[i] * p) + 1e-6)
      
      # posterior predictive replicated data - conditioned on latent state
      y_rep[i, j] ~ dbern(z[i] * p)
      
      # posterior predictive replicated data - not conditioned on latent state
      y_rep_latent[i, j] ~ dbern(z_rep[i] * p)
      
      # deviance - pp data - conditioned on latent state
      D_rep[i, j] <- -2 * 
        log(dbinom(y_rep[i, j], 1, z[i] * p) + 1e-6)
      
      # deviance - pp data - not conditioned on latent state
      D_rep_latent[i, j] <- -2 * 
        log(dbinom(y_rep_latent[i, j], 1, z_rep[i] * p) + 1e-6)
      
    }
    
    # expected values
    y_exp[i, 1:nVisits] <- z[i] * p
    y_exp_rep[i, 1:nVisits] <- z_rep[i] * p
    
  }
  
  ##############################
  # total discrepancy measures #
  ##############################
  
  # deviance
  D_obs_total <- sum(D_obs[1:nSites, 1:nVisits])
  D_rep_total <- sum(D_rep[1:nSites, 1:nVisits])
  D_rep_latent_total <- sum(D_rep_latent[1:nSites, 1:nVisits])

  # chi-squared
  chi_obs_total <- calc_chi(y[1:nSites, 1:nVisits],
                            y_exp[1:nSites, 1:nVisits])
  chi_rep_total <- calc_chi(y_rep[1:nSites, 1:nVisits],
                            y_exp[1:nSites, 1:nVisits])
  chi_rep_latent_total <- calc_chi(y_rep_latent[1:nSites, 1:nVisits],
                                   y_exp_rep[1:nSites, 1:nVisits])

  # likelihood ratio
  ratio_obs_total <- calc_ratio(y[1:nSites, 1:nVisits],
                                y_exp[1:nSites, 1:nVisits])
  ratio_rep_total <- calc_ratio(y_rep[1:nSites, 1:nVisits],
                                y_exp[1:nSites, 1:nVisits])
  ratio_rep_latent_total <- calc_ratio(y_rep_latent[1:nSites, 1:nVisits],
                                       y_exp_rep[1:nSites, 1:nVisits])

  # freeman tukey
  tukey_obs_total <- calc_tukey(y[1:nSites, 1:nVisits],
                                y_exp[1:nSites, 1:nVisits])
  tukey_rep_total <- calc_tukey(y_rep[1:nSites, 1:nVisits],
                                y_exp[1:nSites, 1:nVisits])
  tukey_rep_latent_total <- calc_tukey(y_rep_latent[1:nSites, 1:nVisits],
                                       y_exp_rep[1:nSites, 1:nVisits])

  
})

# 2. covariates on occupancy
model_cov_occ <- nimbleCode({
  
  p ~ dunif(0, 1)
  for (c in 1:ncov) {
    beta[c] ~ dunif(-10, 10)
  }
  
  # logit link - site-level covariates
  psi[1:nSites] <- ilogit(x_site[1:nSites, 1:ncov] %*% beta[1:ncov])
  
  
  for (i in 1:nSites) {
    
    z[i] ~ dbern(psi[i])
    
    # posterior predictive replicated data
    z_rep[i] ~ dbern(psi[i])
    
    for (j in 1:nVisits) {
      
      # observed data
      y[i, j] ~ dbern(z[i] * p)
      
      # log deviance for observed data
      D_obs[i, j] <- -2 * log(dbinom(y[i, j], 1, z[i] * p) + 1e-6) 
      
      # posterior predictive replicated data - conditioned on latent state
      y_rep[i, j] ~ dbern(z[i] * p)
      
      # posterior predictive replicated data - not conditioned on latent state
      y_rep_latent[i, j] ~ dbern(z_rep[i] * p)
      
      # deviance - pp data - conditioned on latent state
      D_rep[i, j] <- -2 * 
        log(dbinom(y_rep[i, j], 1, z[i] * p) + 1e-6)
      
      # deviance - pp data - not conditioned on latent state
      D_rep_latent[i, j] <- -2 * 
        log(dbinom(y_rep_latent[i, j], 1, z_rep[i] * p) + 1e-6)
      
    }
    
    # expected values
    y_exp[i, 1:nVisits] <- z[i] * p
    y_exp_rep[i, 1:nVisits] <- z_rep[i] * p
    
  }
  
  ##############################
  # total discrepancy measures #
  ##############################
  
  # deviance
  D_obs_total <- sum(D_obs[1:nSites, 1:nVisits])
  D_rep_total <- sum(D_rep[1:nSites, 1:nVisits])
  D_rep_latent_total <- sum(D_rep_latent[1:nSites, 1:nVisits])
  
  # chi-squared
  chi_obs_total <- calc_chi(y[1:nSites, 1:nVisits],
                            y_exp[1:nSites, 1:nVisits])
  chi_rep_total <- calc_chi(y_rep[1:nSites, 1:nVisits],
                            y_exp[1:nSites, 1:nVisits])
  chi_rep_latent_total <- calc_chi(y_rep_latent[1:nSites, 1:nVisits],
                                   y_exp_rep[1:nSites, 1:nVisits])
  
  # likelihood ratio
  ratio_obs_total <- calc_ratio(y[1:nSites, 1:nVisits],
                                y_exp[1:nSites, 1:nVisits])
  ratio_rep_total <- calc_ratio(y_rep[1:nSites, 1:nVisits],
                                y_exp[1:nSites, 1:nVisits])
  ratio_rep_latent_total <- calc_ratio(y_rep_latent[1:nSites, 1:nVisits],
                                       y_exp_rep[1:nSites, 1:nVisits])
  
  # freeman tukey
  tukey_obs_total <- calc_tukey(y[1:nSites, 1:nVisits],
                                y_exp[1:nSites, 1:nVisits])
  tukey_rep_total <- calc_tukey(y_rep[1:nSites, 1:nVisits],
                                y_exp[1:nSites, 1:nVisits])
  tukey_rep_latent_total <- calc_tukey(y_rep_latent[1:nSites, 1:nVisits],
                                       y_exp_rep[1:nSites, 1:nVisits])
  
})

# 3. covariates on occupancy and detection
model_cov_occ_det <- nimbleCode({
  
  for (c in 1:ncov_o) {
    beta_o[c] ~ dunif(-10, 10)
  }
  for (c in 1:ncov_p) {
    beta_p[c] ~ dunif(-10, 10)
  }
  
  # logit link - site-level covariates
  psi[1:nSites] <- ilogit(x_site[1:nSites, 1:ncov_o] %*% beta_o[1:ncov_o])
  
  for (i in 1:nSites) {
    
    z[i] ~ dbern(psi[i])
    
    # posterior predictive replicated data
    z_rep[i] ~ dbern(psi[i])
    
    # logit link - detection-level covariates
    p[i, 1:nVisits] <- ilogit(x_visit[i, 1:nVisits, 1:ncov_p] %*%
                                beta_o[1:ncov_p])
    
    for (j in 1:nVisits) {
      
      # observed data
      y[i, j] ~ dbern(z[i] * p[i, j])
      
      # log deviance for observed data
      D_obs[i, j] <- -2 * log(dbinom(y[i, j], 1, z[i] * p[i, j]) + 1e-6) 
      
      # posterior predictive replicated data - conditioned on latent state
      y_rep[i, j] ~ dbern(z[i] * p[i, j])
      
      # posterior predictive replicated data - not conditioned on latent state
      y_rep_latent[i, j] ~ dbern(z_rep[i] * p[i, j])
      
      # deviance - pp data - conditioned on latent state
      D_rep[i, j] <- -2 * 
        log(dbinom(y_rep[i, j], 1, z[i] * p[i, j]) + 1e-6)
      
      # deviance - pp data - not conditioned on latent state
      D_rep_latent[i, j] <- -2 * 
        log(dbinom(y_rep_latent[i, j], 1, z_rep[i] * p[i, j]) + 1e-6)
      
      # expected values
      y_exp[i, j] <- z[i] * p[i, j]
      y_exp_rep[i, j] <- z_rep[i] * p[i, j]
      
    }
  }
  
  ##############################
  # total discrepancy measures #
  ##############################
  
  # deviance
  D_obs_total <- sum(D_obs[1:nSites, 1:nVisits])
  D_rep_total <- sum(D_rep[1:nSites, 1:nVisits])
  D_rep_latent_total <- sum(D_rep_latent[1:nSites, 1:nVisits])
  
  # chi-squared
  chi_obs_total <- calc_chi(y[1:nSites, 1:nVisits],
                            y_exp[1:nSites, 1:nVisits])
  chi_rep_total <- calc_chi(y_rep[1:nSites, 1:nVisits],
                            y_exp[1:nSites, 1:nVisits])
  chi_rep_latent_total <- calc_chi(y_rep_latent[1:nSites, 1:nVisits],
                                   y_exp_rep[1:nSites, 1:nVisits])
  
  # likelihood ratio
  ratio_obs_total <- calc_ratio(y[1:nSites, 1:nVisits],
                                y_exp[1:nSites, 1:nVisits])
  ratio_rep_total <- calc_ratio(y_rep[1:nSites, 1:nVisits],
                                y_exp[1:nSites, 1:nVisits])
  ratio_rep_latent_total <- calc_ratio(y_rep_latent[1:nSites, 1:nVisits],
                                       y_exp_rep[1:nSites, 1:nVisits])
  
  # freeman tukey
  tukey_obs_total <- calc_tukey(y[1:nSites, 1:nVisits],
                                y_exp[1:nSites, 1:nVisits])
  tukey_rep_total <- calc_tukey(y_rep[1:nSites, 1:nVisits],
                                y_exp[1:nSites, 1:nVisits])
  tukey_rep_latent_total <- calc_tukey(y_rep_latent[1:nSites, 1:nVisits],
                                       y_exp_rep[1:nSites, 1:nVisits])
  
})

# 4. spatial random effect
model_spatial_ranef <- nimbleCode({
  
  psi_mean ~ dunif(-10, 10)
  psi_sd ~ dunif(0, 100)
  p ~ dbeta(1, 1)
  
  psi[1:nRegions] <- ilogit(psi_logit[1:nRegions])
  
  for(r in 1:nRegions){
    
    psi_logit[r] ~ dnorm(psi_mean, psi_sd)
    
    # posterior predictive replicated data
    psi_logit_rep[r] ~ dnorm(psi_mean, psi_sd)
    
  }
  
  for(i in 1:nSites) {
    
    z[i] ~ dbern(psi[region[i]])
    
    # posterior predictive replicated data
    z_rep[i] ~ dbern(ilogit(psi_logit_rep[region[i]]))
    
    for(j in 1:nVisits) {
      
      # observed data
      y[i, j] ~ dbern(z[i] * p)
      
      # log deviance for observed data
      D_obs[i, j] <- -2 * log(dbinom(y[i, j], 1, z[i] * p) + 1e-6)  
      
      # posterior predictive replicated data - conditioned on latent state
      y_rep[i, j] ~ dbern(z[i] * p)
      
      # posterior predictive replicated data - not conditioned on latent state
      y_rep_latent[i, j] ~ dbern(z_rep[i] * p)
      
      # deviance - pp data - conditioned on latent state
      D_rep[i, j] <- -2 * 
        log(dbinom(y_rep[i, j], 1, z[i] * p) + 1e-6)
      
      # deviance - pp data - not conditioned on latent state
      D_rep_latent[i, j] <- -2 * 
        log(dbinom(y_rep_latent[i, j], 1, z_rep[i] * p) + 1e-6)
      
    }
    
    # expected values
    y_exp[i, 1:nVisits] <- z[i] * p
    y_exp_rep[i, 1:nVisits] <- z_rep[i] * p
    
  }
  
  ##############################
  # total discrepancy measures #
  ##############################
  
  # deviance
  D_obs_total <- sum(D_obs[1:nSites, 1:nVisits])
  D_rep_total <- sum(D_rep[1:nSites, 1:nVisits])
  D_rep_latent_total <- sum(D_rep_latent[1:nSites, 1:nVisits])
  
  # chi-squared
  chi_obs_total <- calc_chi(y[1:nSites, 1:nVisits],
                            y_exp[1:nSites, 1:nVisits])
  chi_rep_total <- calc_chi(y_rep[1:nSites, 1:nVisits],
                            y_exp[1:nSites, 1:nVisits])
  chi_rep_latent_total <- calc_chi(y_rep_latent[1:nSites, 1:nVisits],
                                   y_exp_rep[1:nSites, 1:nVisits])
  
  # likelihood ratio
  ratio_obs_total <- calc_ratio(y[1:nSites, 1:nVisits],
                                y_exp[1:nSites, 1:nVisits])
  ratio_rep_total <- calc_ratio(y_rep[1:nSites, 1:nVisits],
                                y_exp[1:nSites, 1:nVisits])
  ratio_rep_latent_total <- calc_ratio(y_rep_latent[1:nSites, 1:nVisits],
                                       y_exp_rep[1:nSites, 1:nVisits])
  
  # freeman tukey
  tukey_obs_total <- calc_tukey(y[1:nSites, 1:nVisits],
                                y_exp[1:nSites, 1:nVisits])
  tukey_rep_total <- calc_tukey(y_rep[1:nSites, 1:nVisits],
                                y_exp[1:nSites, 1:nVisits])
  tukey_rep_latent_total <- calc_tukey(y_rep_latent[1:nSites, 1:nVisits],
                                       y_exp_rep[1:nSites, 1:nVisits])
  
})

# 5. basic occupancy model (binomial detection)
model_basic_bin <- nimbleCode({
  
  psi ~ dbeta(1, 1)
  p ~ dbeta(1, 1)
  
  for (i in 1:nSites) {
    
    z[i] ~ dbern(psi) 
    
    # posterior predictive replicated data
    z_rep[i] ~ dbern(psi)
    
    # observed data
    y[i] ~ dbinom(z[i] * p, nVisits)
    
    # log deviance - observed data
    D_obs[i] <- -2 * log(dbinom(y[i], nVisits, z[i] * p) + 1e-6)
    
    # posterior predictive replicated data - conditioned on latent state
    y_rep[i] ~ dbinom(z[i] * p, nVisits)
    
    # posterior predictive replicated data - not conditioned on latent state
    y_rep_latent[i] ~ dbinom(z_rep[i] * p, nVisits)
    
    # deviance - pp data - conditioned on latent state
    D_rep[i] <- -2 * 
      log(dbinom(y_rep[i], nVisits, z[i] * p) + 1e-6)
    
    # deviance - pp data - not conditioned on latent state
    D_rep_latent[i] <- -2 * 
      log(dbinom(y_rep_latent[i], nVisits, z_rep[i] * p) + 1e-6)
    
    # expected values
    y_exp[i] <- z[i] * p * nVisits
    y_exp_rep[i] <- z_rep[i] * p * nVisits
    
  }
  
  ##############################
  # total discrepancy measures #
  ##############################
  
  # deviance
  D_obs_total <- sum(D_obs[1:nSites])
  D_rep_total <- sum(D_rep[1:nSites])
  D_rep_latent_total <- sum(D_rep_latent[1:nSites])
  
  # chi-squared
  chi_obs_total <- calc_chi_betabin(y[1:nSites], y_exp[1:nSites])
  chi_rep_total <- calc_chi_betabin(y_rep[1:nSites], y_exp[1:nSites])
  chi_rep_latent_total <- calc_chi_betabin(y_rep_latent[1:nSites],
                                           y_exp_rep[1:nSites])

  # likelihood ratio
  ratio_obs_total <- calc_ratio_betabin(y[1:nSites], y_exp[1:nSites],
                                        nVisits)
  ratio_rep_total <- calc_ratio_betabin(y_rep[1:nSites], y_exp[1:nSites],
                                        nVisits)
  ratio_rep_latent_total <- calc_ratio_betabin(y_rep_latent[1:nSites],
                                               y_exp_rep[1:nSites],
                                               nVisits)

  # freeman tukey
  tukey_obs_total <- calc_tukey_betabin(y[1:nSites], y_exp[1:nSites])
  tukey_rep_total <- calc_tukey_betabin(y_rep[1:nSites], y_exp[1:nSites])
  tukey_rep_latent_total <- calc_tukey_betabin(y_rep_latent[1:nSites],
                                               y_exp_rep[1:nSites])
  
  
})


######################################
# functions for discrepancy measures #
######################################

calc_chi <- nimbleFunction(
  
  run = function(y = double(2), 
                 y_exp = double(2))
  {
    returnType(double(0))
    
    # calculate chi squared discrepancy measure
    chi_out <- (y - y_exp) ^ 2 / (y + 1e-6)
    
    return(sum(chi_out))
  }
)

calc_ratio <- nimbleFunction(
  
  run = function(y = double(2), 
                 y_exp = double(2))
  {
    returnType(double(0))
    
    # calculate likelihood ratio discrepancy measure
    ratio_out <- 2 * (y * log((y + 1e-6) / (y_exp + 1e-6)) + (1 - y) * 
                        log((1 - y + 1e-6)/(1 - y_exp + 1e-6)))
    
    return(sum(ratio_out))
  }
)

calc_tukey <- nimbleFunction(
  
  run = function(y = double(2), 
                 y_exp = double(2))
  {
    returnType(double(0))
    
    # calculate freeman tukey discrepancy measure
    tukey_out <- (sqrt(y) - sqrt(y_exp)) ^ 2
    
    return(sum(tukey_out))
  }
)

calc_chi_betabin <- nimbleFunction(
  
  run = function(y = double(1), 
                 y_exp = double(1))
  {
    returnType(double(0))
    
    # calculate chi squared discrepancy measure
    chi_out <- (y - y_exp) ^ 2 / (y + 1e-6)
    
    return(sum(chi_out))
  }
)

calc_ratio_betabin <- nimbleFunction(
  
  run = function(y = double(1), 
                 y_exp = double(1),
                 nVisits = double(0))
  {
    returnType(double(0))
    
    # calculate likelihood ratio discrepancy measure
    ratio_out <- 2 * (y * log((y + 1e-6) / (y_exp + 1e-6)) + (nVisits - y) * 
                        log((nVisits - y + 1e-6)/(nVisits - y_exp + 1e-6)))
    
    return(sum(ratio_out))
  }
)

calc_tukey_betabin <- nimbleFunction(
  
  run = function(y = double(1), 
                 y_exp = double(1))
  {
    returnType(double(0))
    
    # calculate freeman tukey discrepancy measure
    tukey_out <- (sqrt(y) - sqrt(y_exp)) ^ 2
    
    return(sum(tukey_out))
  }
)

##############################
# functions for running MCMC #
##############################

fit_basic <- function(compiled_model, compiled_mcmc, niter, nburnin, thin) {
  
  # function for generating initial values
  inits <- function(y) list(psi = runif(1, 0, 1), p = runif(1, 0, 1), 
                            z = pmin(rowSums(y), 1),
                            z_rep = rep(NA, nrow(y)),
                            y_rep = matrix(NA, nrow(y), ncol(y)),
                            y_rep_latent = matrix(NA, nrow(y), ncol(y))) 
  
  # add initial values
  compiled_model$setInits(inits(compiled_model$y))
  
  # generate samples
  samples <- runMCMC(compiled_mcmc, niter = niter, 
                     nburnin = nburnin, thin = thin)
  
  return(samples)
}


fit_cov_occ <- function(compiled_model, compiled_mcmc, ncov, 
                        niter, nburnin, thin) {
  
  # function for generating initial values
  inits <- function(y, ncov) list(beta = rnorm(ncov, 0, 1), 
                                  p = runif(1, 0, 1), 
                                  z = pmin(rowSums(y), 1),
                                  z_rep = rep(NA, nrow(y)),
                                  y_rep = matrix(NA, nrow(y), ncol(y)),
                                  y_rep_latent = matrix(NA, nrow(y), ncol(y)))  
  
  # add initial values
  compiled_model$setInits(inits(compiled_model$y, ncov))
  
  # generate samples
  samples <- runMCMC(compiled_mcmc, niter = niter, 
                     nburnin = nburnin, thin = thin)
  
  return(samples)
}


fit_cov_occ_det <- function(compiled_model, compiled_mcmc, ncov_o, 
                            ncov_p, niter, nburnin, thin) {
  
  # function for generating initial values
  inits <- function(y, ncov_o, ncov_p) list(beta_o = rnorm(ncov_o, 0, 1), 
                                            beta_p = rnorm(ncov_p, 0, 1), 
                                            z_rep = rep(NA, nrow(y)),
                                            z = pmin(rowSums(y), 1),
                                            y_rep = matrix(NA, nrow(y), 
                                                           ncol(y)),
                                            y_rep_latent = matrix(NA, nrow(y), 
                                                                  ncol(y)))  
  
  # add initial values
  compiled_model$setInits(inits(compiled_model$y, ncov_o, ncov_p))
  
  # generate samples
  samples <- runMCMC(compiled_mcmc, niter = niter, 
                     nburnin = nburnin, thin = thin)
  
  return(samples)
}

fit_spatial_ranef <- function(compiled_model, compiled_mcmc, region,
                              nRegions, niter, nburnin, thin) {
  
  # function for generating initial values
  inits <- function(y, nRegions) list(
    psi_mean = 0, psi_sd = 1, p = runif(1, 0, 1),
    psi_logit = logit(runif(nRegions, 0, 1)), 
    z = pmin(rowSums(y), 1), z_rep = rep(NA, nrow(y)),
    y_rep = matrix(NA, nrow(y), ncol(y)),
    y_rep_latent = matrix(NA, nrow(y), ncol(y))
  )  
  
  # add initial values
  compiled_model$setInits(inits(compiled_model$y, nRegions))
  
  # generate samples
  samples <- runMCMC(compiled_mcmc, niter = niter, 
                     nburnin = nburnin, thin = thin)
  
  return(samples)
}

fit_basic_betabin <- function(compiled_model, compiled_mcmc, 
                              niter, nburnin, thin) {
  
  # function for generating initial values
  inits <- function(y) list(psi = runif(1, 0, 1), p = runif(1, 0, 1), 
                            z = pmin(y, 1),
                            z_rep = rep(NA, length(y)),
                            y_rep = rep(NA, length(y)),
                            y_rep_latent = rep(NA, length(y))) 
  
  # add initial values
  compiled_model$setInits(inits(compiled_model$y))
  
  # generate samples
  samples <- runMCMC(compiled_mcmc, niter = niter, 
                     nburnin = nburnin, thin = thin)
  
  return(samples)
}

