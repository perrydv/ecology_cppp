
####################
# model variations #
####################


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
      
      # log deviance - observed data
      D_obs[i, j] <- -2 * log(dbern(y[i, j], z[i] * p) + 1e-6)
      
      # posterior predictive replicated data
      y_rep[i, j] ~ dbern(z[i] * p)
      
      # deviance - pp data
      D_rep[i, j] <- -2 * log(dbern(y_rep[i, j], z[i] * p) + 1e-6)
      
    }
    
    # expected values
    y_exp[i, 1:nVisits] <- z[i] * p
    
  }
  
  ##############################
  # total discrepancy measures #
  ##############################
  
  # deviance
  D_obs_total <- sum(D_obs[1:nSites, 1:nVisits])
  D_rep_total <- sum(D_rep[1:nSites, 1:nVisits])
  
  # chi-squared
  chi_obs_total <- calc_chi(y[1:nSites, 1:nVisits], 
                            y_exp[1:nSites, 1:nVisits])
  chi_rep_total <- calc_chi(y_rep[1:nSites, 1:nVisits], 
                            y_exp[1:nSites, 1:nVisits])
  
  # likelihood ratio
  ratio_obs_total <- calc_ratio(y[1:nSites, 1:nVisits], 
                                y_exp[1:nSites, 1:nVisits])
  ratio_rep_total <- calc_ratio(y_rep[1:nSites, 1:nVisits], 
                                y_exp[1:nSites, 1:nVisits])
  
  # freeman tukey
  tukey_obs_total <- calc_tukey(y[1:nSites, 1:nVisits], 
                                y_exp[1:nSites, 1:nVisits])
  tukey_rep_total <- calc_tukey(y_rep[1:nSites, 1:nVisits], 
                                y_exp[1:nSites, 1:nVisits])
  
  
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
      
      # deviance - pp data
      D_rep[i, j] <- -2 * log(dbern(y_rep[i, j], z[i] * p) + 1e-6)
      
    }
    
    # expected values
    y_exp[i, 1:nVisits] <- z[i] * p
    
  }
  
  ##############################
  # total discrepancy measures #
  ##############################
  
  # deviance
  D_obs_total <- sum(D_obs[1:nSites, 1:nVisits])
  D_rep_total <- sum(D_rep[1:nSites, 1:nVisits])
  
  # chi-squared
  chi_obs_total <- calc_chi(y[1:nSites, 1:nVisits], 
                            y_exp[1:nSites, 1:nVisits])
  chi_rep_total <- calc_chi(y_rep[1:nSites, 1:nVisits], 
                            y_exp[1:nSites, 1:nVisits])
  
  # likelihood ratio
  ratio_obs_total <- calc_ratio(y[1:nSites, 1:nVisits], 
                                y_exp[1:nSites, 1:nVisits])
  ratio_rep_total <- calc_ratio(y_rep[1:nSites, 1:nVisits], 
                                y_exp[1:nSites, 1:nVisits])
  
  # freeman tukey
  tukey_obs_total <- calc_tukey(y[1:nSites, 1:nVisits], 
                                y_exp[1:nSites, 1:nVisits])
  tukey_rep_total <- calc_tukey(y_rep[1:nSites, 1:nVisits], 
                                y_exp[1:nSites, 1:nVisits])
  
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
      
      # deviance - pp data
      D_rep[i, j] <- -2 * log(dbern(y_rep[i, j], z[i] * p[i, j]) + 1e-6)
      
      # expected values
      y_exp[i, j] <- z[i] * p[i, j]
      
    }
  }
  
  ##############################
  # total discrepancy measures #
  ##############################
  
  # deviance
  D_obs_total <- sum(D_obs[1:nSites, 1:nVisits])
  D_rep_total <- sum(D_rep[1:nSites, 1:nVisits])
  
  # chi-squared
  chi_obs_total <- calc_chi(y[1:nSites, 1:nVisits], 
                            y_exp[1:nSites, 1:nVisits])
  chi_rep_total <- calc_chi(y_rep[1:nSites, 1:nVisits], 
                            y_exp[1:nSites, 1:nVisits])
  
  # likelihood ratio
  ratio_obs_total <- calc_ratio(y[1:nSites, 1:nVisits], 
                                y_exp[1:nSites, 1:nVisits])
  ratio_rep_total <- calc_ratio(y_rep[1:nSites, 1:nVisits], 
                                y_exp[1:nSites, 1:nVisits])
  
  # freeman tukey
  tukey_obs_total <- calc_tukey(y[1:nSites, 1:nVisits], 
                                y_exp[1:nSites, 1:nVisits])
  tukey_rep_total <- calc_tukey(y_rep[1:nSites, 1:nVisits], 
                                y_exp[1:nSites, 1:nVisits])
  
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
      
      # deviance - pp data
      D_rep[i, j] <- -2 * log(dbern(y_rep[i, j], z[i] * p) + 1e-6)
      
    }
    
    # expected values
    y_exp[i, 1:nVisits] <- z[i] * p
    
  }
  
  ##############################
  # total discrepancy measures #
  ##############################
  
  # deviance
  D_obs_total <- sum(D_obs[1:nSites, 1:nVisits])
  D_rep_total <- sum(D_rep[1:nSites, 1:nVisits])
  
  # chi-squared
  chi_obs_total <- calc_chi(y[1:nSites, 1:nVisits], 
                            y_exp[1:nSites, 1:nVisits])
  chi_rep_total <- calc_chi(y_rep[1:nSites, 1:nVisits], 
                            y_exp[1:nSites, 1:nVisits])
  
  # likelihood ratio
  ratio_obs_total <- calc_ratio(y[1:nSites, 1:nVisits], 
                                y_exp[1:nSites, 1:nVisits])
  ratio_rep_total <- calc_ratio(y_rep[1:nSites, 1:nVisits], 
                                y_exp[1:nSites, 1:nVisits])
  
  # freeman tukey
  tukey_obs_total <- calc_tukey(y[1:nSites, 1:nVisits], 
                                y_exp[1:nSites, 1:nVisits])
  tukey_rep_total <- calc_tukey(y_rep[1:nSites, 1:nVisits], 
                                y_exp[1:nSites, 1:nVisits])
  
})

# 5. multi-species (https://doi.org/10.1111/2041-210X.12587)
# model_multispec <- nimbleCode({
  
# })


######################################
# functions for discrepancy measures #
######################################

calc_chi <- nimbleFunction(
  
  run = function(y = double(2), 
                 y_exp = double(2))
  {
    returnType(double(0))
    
    # calculate chi squared discrepancy measure
    chi_out <- (y - y_exp) ^ 2 / y
    
    return(sum(chi_out))
  }
)
assign("calc_chi", calc_chi, envir = .GlobalEnv)

calc_ratio <- nimbleFunction(
  
  run = function(y = double(2), 
                 y_exp = double(2))
  {
    returnType(double(0))
    
    # calculate likelihood ratio discrepancy measure
    ratio_out <- 2 * y * log(y / y_exp + 1e-6)
    
    return(sum(ratio_out))
  }
)
assign("calc_ratio", calc_ratio, envir = .GlobalEnv)

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
assign("calc_tukey", calc_tukey, envir = .GlobalEnv)

##############################
# functions for running MCMC #
##############################

fit_basic <- function(model_code, y, nSites, nVisits, 
                      niter, nburnin, thin) {
  
  constants <- list(nSites = nSites, nVisits = nVisits)
  
  data <- list(y = y)
  
  inits <- function() list(psi = runif(1, 0, 1), p = runif(1, 0, 1), 
                           z = pmin(rowSums(y), 1))  
  
  model <- nimbleModel(model_code, constants = constants, 
                       data = data, inits = inits)
  
  compiled_model <- compileNimble(model)
  
  mcmc_conf <- configureMCMC(model, 
                             monitors = c("psi", "p", "z",
                                          "D_obs_total", "D_rep_total",
                                          "chi_obs_total", "chi_rep_total",
                                          "ratio_obs_total", "ratio_rep_total",
                                          "tukey_obs_total", "tukey_rep_total"))
  mcmc <- buildMCMC(mcmc_conf)
  compiled_mcmc <- compileNimble(mcmc, project = compiled_model)
  
  samples <- runMCMC(compiled_mcmc, niter = niter, 
                     nburnin = nburnin, thin = thin)
  
  return(samples)
}


fit_cov_occ <- function(model_code, y, nSites, nVisits, ncov, 
                        niter, nburnin, thin) {
  
  constants <- list(nSites = nSites, nVisits = nVisits)
  
  data <- list(y = y)
  
  inits <- function() list(beta = rnorm(ncov, 0, 1), p = runif(1, 0, 1), 
                           z = pmin(rowSums(y), 1))  
  
  model <- nimbleModel(model_code, constants = constants, 
                       data = data, inits = inits)
  
  compiled_model <- compileNimble(model)
  
  mcmc_conf <- configureMCMC(model, 
                             monitors = c("beta", "p", "z",
                                          "D_obs_total", "D_rep_total",
                                          "chi_obs_total", "chi_rep_total",
                                          "ratio_obs_total", "ratio_rep_total",
                                          "tukey_obs_total", "tukey_rep_total"))
  mcmc <- buildMCMC(mcmc_conf)
  compiled_mcmc <- compileNimble(mcmc, project = compiled_model)
  
  samples <- runMCMC(compiled_mcmc, niter = niter, 
                     nburnin = nburnin, thin = thin)
  
  return(samples)
}


fit_cov_occ_det <- function(model_code, y, nSites, nVisits, ncov_o, ncov_p, 
                            niter, nburnin, thin) {
  
  constants <- list(nSites = nSites, nVisits = nVisits)
  
  data <- list(y = y)
  
  inits <- function() list(beta_o = rnorm(ncov_o, 0, 1), 
                           beta_p = rnorm(ncov_p, 0, 1), 
                           z = pmin(rowSums(y), 1))  
  
  model <- nimbleModel(model_code, constants = constants, 
                       data = data, inits = inits)
  
  compiled_model <- compileNimble(model)
  
  mcmc_conf <- configureMCMC(model, 
                             monitors = c("beta", "p", "z",
                                          "D_obs_total", "D_rep_total",
                                          "chi_obs_total", "chi_rep_total",
                                          "ratio_obs_total", "ratio_rep_total",
                                          "tukey_obs_total", "tukey_rep_total"))
  mcmc <- buildMCMC(mcmc_conf)
  compiled_mcmc <- compileNimble(mcmc, project = compiled_model)
  
  samples <- runMCMC(compiled_mcmc, niter = niter, 
                     nburnin = nburnin, thin = thin)
  
  return(samples)
}

fit_spatial_ranef <- function(model_code, y, nSites, nVisits, nRegions, 
                              niter, nburnin, thin) {
  
  constants <- list(nSites = nSites, nVisits = nVisits)
  
  data <- list(y = y)
  
  inits <- function() list(psi_mean = 0, psi_sd = 1, p = runif(1, 0, 1),
                           psi_logit = logit(runif(nRegions, 0, 1)),
                           z = pmin(rowSums(y), 1))  
  
  model <- nimbleModel(model_code, constants = constants, 
                       data = data, inits = inits)
  
  compiled_model <- compileNimble(model)
  
  mcmc_conf <- configureMCMC(model, 
                             monitors = c("beta", "p", "z",
                                          "D_obs_total", "D_rep_total",
                                          "chi_obs_total", "chi_rep_total",
                                          "ratio_obs_total", "ratio_rep_total",
                                          "tukey_obs_total", "tukey_rep_total"))
  mcmc <- buildMCMC(mcmc_conf)
  compiled_mcmc <- compileNimble(mcmc, project = compiled_model)
  
  samples <- runMCMC(compiled_mcmc, niter = niter, 
                     nburnin = nburnin, thin = thin)
  
  return(samples)
}

