library(nimble)
library(parallel)

# load models
source("occupancy_model.R")
source("occupancy_model_latent.R")

# load functions to simulate data
source("simulate_data.R")

#############
# constants #
#############

n_datasets <- 200
niter <- 5000
nburnin <- 1000
thin <- 5

# data size
nRegions <- 10
nSites <- nRegions * 5
nVisits <- 6
nCov_site <- 2
nCov_detect <- 2


#################################
# non-independence across sites #
#################################

# parameter values
p <- 0.3
psi <- 0.6
psi_sd1 <- 2.5
psi_sd2 <- 5
psi_sd3 <- 10

# scenario 1: sd = 2.5
pvalues_scen1 <- matrix(NA, nrow = n_datasets, ncol = 4)
colnames(p_values_basic) <- c("deviance", "chi_sq", "lik_ratio", "ftukey")


##########################################
# generate null distribution of p values #
##########################################

###################################################
# model 1 - basic, conditioned on latent states 

p_values_basic <- matrix(NA, nrow = n_datasets, ncol = 4)
colnames(p_values_basic) <- c("deviance", "chi_sq", "lik_ratio", "ftukey")

# create model instance and compile
constants <- list(nSites = nSites, nVisits = nVisits)
model <- nimbleModel(model_basic, constants = constants, 
                     data = list(y = simulate_basic(
                       params = list(p = p, psi = psi), nSites, nVisits)))
compiled_model <- compileNimble(model)

# configure MCMC
mcmc_conf <- configureMCMC(model, 
                           monitors = c("psi", "p", "z",
                                        "D_obs_total", "D_rep_total",
                                        "chi_obs_total", "chi_rep_total",
                                        "ratio_obs_total", "ratio_rep_total",
                                        "tukey_obs_total", "tukey_rep_total"
                                        ))

# build MCMC
mcmc <- buildMCMC(mcmc_conf)

# compile mcmc
compiled_mcmc <- compileNimble(mcmc, project = compiled_model)

for (i in 1:n_datasets) {
  
  # simulate data
  y <- simulate_basic(
    params = list(p = p, psi = psi), 
    nSites, nVisits
  )
  compiled_model$y <- y
  
  # fit model
  samples <- fit_basic(compiled_model, compiled_mcmc,
                       niter, nburnin, thin)
  
  # calculate p values
  p_values_basic[i, "deviance"] <- mean(
    samples[, "D_rep_total"] > samples[, "D_obs_total"]
  )
  p_values_basic[i, "chi_sq"] <- mean(
    samples[, "chi_rep_total"] > samples[, "chi_obs_total"]
  )
  p_values_basic[i, "lik_ratio"] <- mean(
    samples[, "ratio_rep_total"] > samples[, "ratio_obs_total"]
  )
  p_values_basic[i, "ftukey"] <- mean(
    samples[, "tukey_rep_total"] > samples[, "tukey_obs_total"]
  )
  
}
saveRDS(p_values_basic, "pvalues/null/basic.rds")

###################################################
# model 1 - basic, not conditioned on latent states 

p_values_basic <- matrix(NA, nrow = n_datasets, ncol = 4)
colnames(p_values_basic) <- c("deviance", "chi_sq", "lik_ratio", "ftukey")

# create model instance and compile
constants <- list(nSites = nSites, nVisits = nVisits)
model <- nimbleModel(model_basic_latent, constants = constants, 
                     data = list(y = simulate_basic(
                       params = list(p = p, psi = psi), nSites, nVisits)))
compiled_model <- compileNimble(model)

# configure MCMC
mcmc_conf <- configureMCMC(model, 
                           monitors = c("psi", "p", "z",
                                        "D_obs_total", "D_rep_total",
                                        "chi_obs_total", "chi_rep_total",
                                        "ratio_obs_total", "ratio_rep_total",
                                        "tukey_obs_total", "tukey_rep_total"
                           ))

# build MCMC
mcmc <- buildMCMC(mcmc_conf)

# compile mcmc
compiled_mcmc <- compileNimble(mcmc, project = compiled_model)

for (i in 1:n_datasets) {
  
  # simulate data
  y <- simulate_basic(
    params = list(p = p, psi = psi), 
    nSites, nVisits
  )
  compiled_model$y <- y
  
  # fit model
  samples <- fit_basic_latent(compiled_model, compiled_mcmc,
                              niter, nburnin, thin)
  
  # calculate p values
  p_values_basic[i, "deviance"] <- mean(
    samples[, "D_rep_total"] > samples[, "D_obs_total"]
  )
  p_values_basic[i, "chi_sq"] <- mean(
    samples[, "chi_rep_total"] > samples[, "chi_obs_total"]
  )
  p_values_basic[i, "lik_ratio"] <- mean(
    samples[, "ratio_rep_total"] > samples[, "ratio_obs_total"]
  )
  p_values_basic[i, "ftukey"] <- mean(
    samples[, "tukey_rep_total"] > samples[, "tukey_obs_total"]
  )
  
}
saveRDS(p_values_basic, "pvalues/null/basic_latent.rds")


##############################################################
# model 2 - occupancy covariates, conditioned on latent states 

p_values_cov_occ <- matrix(NA, nrow = n_datasets, ncol = 4)
colnames(p_values_cov_occ) <- c("deviance", "chi_sq", "lik_ratio", "ftukey")

# simulate data to build model
x_site <- matrix(NA, nrow = nSites, ncol = nCov_site + 1) # covariate data
x_site[, 1] <- 1
for (j in 1:nCov_site) {
  x_site[, j + 1] <- rnorm(nSites, 0, 1)
}

y <- simulate_cov_occ(
  params = list(beta = runif(nCov_site + 1, -1, 1), p = p), 
  x_site, nSites, nVisits
)

# create model instance and compile
constants <- list(nSites = nSites, nVisits = nVisits, ncov = nCov_site + 1)
model <- nimbleModel(model_cov_occ, constants = constants,
                     data = list(x_site = x_site, y = y))
compiled_model <- compileNimble(model)

# configure MCMC
mcmc_conf <- configureMCMC(model, 
                           monitors = c("beta", "p", "z",
                                        "D_obs_total", "D_rep_total",
                                        "chi_obs_total", "chi_rep_total",
                                        "ratio_obs_total", "ratio_rep_total",
                                        "tukey_obs_total", "tukey_rep_total"))

# build and compile MCMC
mcmc <- buildMCMC(mcmc_conf)
compiled_mcmc <- compileNimble(mcmc, project = compiled_model)

for (i in 1:n_datasets) {
  
  # simulate data
  x_site <- matrix(NA, nrow = nSites, ncol = nCov_site + 1) # covariate data
  x_site[, 1] <- 1
  for (j in 1:nCov_site) {
    x_site[, j + 1] <- rnorm(nSites, 0, 1)
  }
  
  y <- simulate_cov_occ(
    params = list(beta = runif(nCov_site + 1, -1, 1), p = p), 
    x_site, nSites, nVisits
  )
  
  # add data
  compiled_model$y <- y
  compiled_model$x_site <- x_site
  
  # fit model
  samples <- fit_cov_occ(compiled_model, compiled_mcmc, nCov_site + 1,
                         niter, nburnin, thin)
  
  # calculate p values
  p_values_cov_occ[i, "deviance"] <- mean(
    samples[, "D_rep_total"] > samples[, "D_obs_total"]
  )
  p_values_cov_occ[i, "chi_sq"] <- mean(
    samples[, "chi_rep_total"] > samples[, "chi_obs_total"]
  )
  p_values_cov_occ[i, "lik_ratio"] <- mean(
    samples[, "ratio_rep_total"] > samples[, "ratio_obs_total"]
  )
  p_values_cov_occ[i, "ftukey"] <- mean(
    samples[, "tukey_rep_total"] > samples[, "tukey_obs_total"]
  )
  
}
saveRDS(p_values_cov_occ, "pvalues/null/cov_occ.rds")

##################################################################
# model 2 - occupancy covariates, not conditioned on latent states 

p_values_cov_occ <- matrix(NA, nrow = n_datasets, ncol = 4)
colnames(p_values_cov_occ) <- c("deviance", "chi_sq", "lik_ratio", "ftukey")

# simulate data to build model
x_site <- matrix(NA, nrow = nSites, ncol = nCov_site + 1) # covariate data
x_site[, 1] <- 1
for (j in 1:nCov_site) {
  x_site[, j + 1] <- rnorm(nSites, 0, 1)
}

y <- simulate_cov_occ(
  params = list(beta = runif(nCov_site + 1, -1, 1), p = p), 
  x_site, nSites, nVisits
)

# create model instance and compile
constants <- list(nSites = nSites, nVisits = nVisits, ncov = nCov_site + 1)
model <- nimbleModel(model_cov_occ_latent, constants = constants,
                     data = list(x_site = x_site, y = y))
compiled_model <- compileNimble(model)

# configure MCMC
mcmc_conf <- configureMCMC(model, 
                           monitors = c("beta", "p", "z",
                                        "D_obs_total", "D_rep_total",
                                        "chi_obs_total", "chi_rep_total",
                                        "ratio_obs_total", "ratio_rep_total",
                                        "tukey_obs_total", "tukey_rep_total"))

# build and compile MCMC
mcmc <- buildMCMC(mcmc_conf)
compiled_mcmc <- compileNimble(mcmc, project = compiled_model)

for (i in 1:n_datasets) {
  
  # simulate data
  x_site <- matrix(NA, nrow = nSites, ncol = nCov_site + 1) # covariate data
  x_site[, 1] <- 1
  for (j in 1:nCov_site) {
    x_site[, j + 1] <- rnorm(nSites, 0, 1)
  }
  
  y <- simulate_cov_occ(
    params = list(beta = runif(nCov_site + 1, -1, 1), p = p), 
    x_site, nSites, nVisits
  )
  
  # add data
  compiled_model$y <- y
  compiled_model$x_site <- x_site
  
  # fit model
  samples <- fit_cov_occ_latent(compiled_model, compiled_mcmc, nCov_site + 1,
                                niter, nburnin, thin)
  
  # calculate p values
  p_values_cov_occ[i, "deviance"] <- mean(
    samples[, "D_rep_total"] > samples[, "D_obs_total"]
  )
  p_values_cov_occ[i, "chi_sq"] <- mean(
    samples[, "chi_rep_total"] > samples[, "chi_obs_total"]
  )
  p_values_cov_occ[i, "lik_ratio"] <- mean(
    samples[, "ratio_rep_total"] > samples[, "ratio_obs_total"]
  )
  p_values_cov_occ[i, "ftukey"] <- mean(
    samples[, "tukey_rep_total"] > samples[, "tukey_obs_total"]
  )
  
}
saveRDS(p_values_cov_occ, "pvalues/null/cov_occ_latent.rds")


################################################################################
# model 3 - occupancy and detection covariates, not conditioned on latent states

p_values_cov_occ_det <- matrix(NA, nrow = n_datasets, ncol = 4)
colnames(p_values_cov_occ_det) <- c("deviance", "chi_sq", "lik_ratio", "ftukey")

# simulate data to build model
x_site <- matrix(NA, nrow = nSites, ncol = nCov_site + 1) 
x_site[, 1] <- 1
for (j in 1:nCov_site) {
  x_site[, j + 1] <- rnorm(nSites, 0, 1)
}

# simulate data - visit-level covariates
x_visit <- array(NA, dim = c(nSites, nVisits, nCov_detect + 1))
x_visit[, , 1] <- 1
for (j in 1:nCov_detect) {
  x_visit[, , j + 1] <- matrix(rnorm(nSites * nVisits, 0, 1), nrow = nSites)
}

# simulate y
y <- simulate_cov_occ_det(
  params = list(
    beta_o = runif(nCov_site + 1, -1, 1),
    beta_p = runif(nCov_detect + 1, -1, 1)
  ), 
  x_site, x_visit, nSites, nVisits
)

# create model instance and compile
constants <- list(nSites = nSites, nVisits = nVisits,
                  ncov_o = nCov_site + 1, ncov_p = nCov_detect + 1)
model <- nimbleModel(model_cov_occ_det, constants = constants,
                     data = list(y = y, x_site = x_site, x_visit = x_visit))
compiled_model <- compileNimble(model)

# configure MCMC
mcmc_conf <- configureMCMC(model, 
                           monitors = c("beta_o", "beta_p", "p", "z",
                                        "D_obs_total", "D_rep_total",
                                        "chi_obs_total", "chi_rep_total",
                                        "ratio_obs_total", "ratio_rep_total",
                                        "tukey_obs_total", "tukey_rep_total"))
# build and compile MCMC
mcmc <- buildMCMC(mcmc_conf)
compiled_mcmc <- compileNimble(mcmc, project = compiled_model)

for (i in 1:n_datasets) {
  
  # simulate data - site-level covariates
  x_site <- matrix(NA, nrow = nSites, ncol = nCov_site + 1) 
  x_site[, 1] <- 1
  for (j in 1:nCov_site) {
    x_site[, j + 1] <- rnorm(nSites, 0, 1)
  }
  
  # simulate data - visit-level covariates
  x_visit <- array(NA, dim = c(nSites, nVisits, nCov_detect + 1))
  x_visit[, , 1] <- 1
  for (j in 1:nCov_detect) {
    x_visit[, , j + 1] <- matrix(rnorm(nSites * nVisits, 0, 1), nrow = nSites)
  }
  
  # simulate y
  y <- simulate_cov_occ_det(
    params = list(
      beta_o = runif(nCov_site + 1, -1, 1),
      beta_p = runif(nCov_detect + 1, -1, 1)
    ), 
    x_site, x_visit, nSites, nVisits
  )
  
  # add data
  compiled_model$y <- y
  compiled_model$x_site <- x_site
  compiled_model$x_visit <- x_visit
  
  
  # fit model
  samples <- fit_cov_occ_det(compiled_model, compiled_mcmc,
                             nCov_site + 1, nCov_detect + 1, 
                             niter, nburnin, thin)
  
  # calculate p values
  p_values_cov_occ_det[i, "deviance"] <- mean(
    samples[, "D_rep_total"] > samples[, "D_obs_total"]
  )
  p_values_cov_occ_det[i, "chi_sq"] <- mean(
    samples[, "chi_rep_total"] > samples[, "chi_obs_total"]
  )
  p_values_cov_occ_det[i, "lik_ratio"] <- mean(
    samples[, "ratio_rep_total"] > samples[, "ratio_obs_total"]
  )
  p_values_cov_occ_det[i, "ftukey"] <- mean(
    samples[, "tukey_rep_total"] > samples[, "tukey_obs_total"]
  )
  
}
saveRDS(p_values_cov_occ_det, "pvalues/null/cov_occ_det.rds")

############################################################################
# model 3 - occupancy and detection covariates, conditioned on latent states

p_values_cov_occ_det <- matrix(NA, nrow = n_datasets, ncol = 4)
colnames(p_values_cov_occ_det) <- c("deviance", "chi_sq", "lik_ratio", "ftukey")

# simulate data to build model
x_site <- matrix(NA, nrow = nSites, ncol = nCov_site + 1) 
x_site[, 1] <- 1
for (j in 1:nCov_site) {
  x_site[, j + 1] <- rnorm(nSites, 0, 1)
}

# simulate data - visit-level covariates
x_visit <- array(NA, dim = c(nSites, nVisits, nCov_detect + 1))
x_visit[, , 1] <- 1
for (j in 1:nCov_detect) {
  x_visit[, , j + 1] <- matrix(rnorm(nSites * nVisits, 0, 1), nrow = nSites)
}

# simulate y
y <- simulate_cov_occ_det(
  params = list(
    beta_o = runif(nCov_site + 1, -1, 1),
    beta_p = runif(nCov_detect + 1, -1, 1)
  ), 
  x_site, x_visit, nSites, nVisits
)

# create model instance and compile
constants <- list(nSites = nSites, nVisits = nVisits,
                  ncov_o = nCov_site + 1, ncov_p = nCov_detect + 1)
model <- nimbleModel(model_cov_occ_det_latent, constants = constants,
                     data = list(y = y, x_site = x_site, x_visit = x_visit))
compiled_model <- compileNimble(model)

# configure MCMC
mcmc_conf <- configureMCMC(model, 
                           monitors = c("beta_o", "beta_p", "p", "z",
                                        "D_obs_total", "D_rep_total",
                                        "chi_obs_total", "chi_rep_total",
                                        "ratio_obs_total", "ratio_rep_total",
                                        "tukey_obs_total", "tukey_rep_total"))
# build and compile MCMC
mcmc <- buildMCMC(mcmc_conf)
compiled_mcmc <- compileNimble(mcmc, project = compiled_model)

for (i in 1:n_datasets) {
  
  # simulate data - site-level covariates
  x_site <- matrix(NA, nrow = nSites, ncol = nCov_site + 1) 
  x_site[, 1] <- 1
  for (j in 1:nCov_site) {
    x_site[, j + 1] <- rnorm(nSites, 0, 1)
  }
  
  # simulate data - visit-level covariates
  x_visit <- array(NA, dim = c(nSites, nVisits, nCov_detect + 1))
  x_visit[, , 1] <- 1
  for (j in 1:nCov_detect) {
    x_visit[, , j + 1] <- matrix(rnorm(nSites * nVisits, 0, 1), nrow = nSites)
  }
  
  # simulate y
  y <- simulate_cov_occ_det(
    params = list(
      beta_o = runif(nCov_site + 1, -1, 1),
      beta_p = runif(nCov_detect + 1, -1, 1)
    ), 
    x_site, x_visit, nSites, nVisits
  )
  
  # add data
  compiled_model$y <- y
  compiled_model$x_site <- x_site
  compiled_model$x_visit <- x_visit
  
  
  # fit model
  samples <- fit_cov_occ_det_latent(compiled_model, compiled_mcmc,
                                    nCov_site + 1, nCov_detect + 1, 
                                    niter, nburnin, thin)
  
  # calculate p values
  p_values_cov_occ_det[i, "deviance"] <- mean(
    samples[, "D_rep_total"] > samples[, "D_obs_total"]
  )
  p_values_cov_occ_det[i, "chi_sq"] <- mean(
    samples[, "chi_rep_total"] > samples[, "chi_obs_total"]
  )
  p_values_cov_occ_det[i, "lik_ratio"] <- mean(
    samples[, "ratio_rep_total"] > samples[, "ratio_obs_total"]
  )
  p_values_cov_occ_det[i, "ftukey"] <- mean(
    samples[, "tukey_rep_total"] > samples[, "tukey_obs_total"]
  )
  
}
saveRDS(p_values_cov_occ_det, "pvalues/null/cov_occ_det_latent.rds")



#######################################################
# model 4 - spatial ranef, conditioned on latent states

p_values_sp_ranef <- matrix(NA, nrow = n_datasets, ncol = 4)
colnames(p_values_sp_ranef) <- c("deviance", "chi_sq", "lik_ratio", "ftukey")

# simulate data to build model
region <- rep(1:nRegions, each = nSites / nRegions)
y <- simulate_spatial_ranef(
  params = list(
    psi_mean = psi, #(logit scale)
    psi_sd = psi_sd, #(logit scale)
    p = p
  ), 
  nRegions, region,
  nSites, nVisits
)

# create model instance and compile
constants <- list(nSites = nSites, nVisits = nVisits, 
                  nRegions = nRegions, region = region)
model <- nimbleModel(model_spatial_ranef, constants = constants,
                     data = list(y = y))
compiled_model <- compileNimble(model)

# configure MCMC
mcmc_conf <- configureMCMC(model, 
                           monitors = c("psi_mean", "psi_sd", "p", "z",
                                        "D_obs_total", "D_rep_total",
                                        "chi_obs_total", "chi_rep_total",
                                        "ratio_obs_total", "ratio_rep_total",
                                        "tukey_obs_total", "tukey_rep_total"))

# build and compile MCMC
mcmc <- buildMCMC(mcmc_conf)
compiled_mcmc <- compileNimble(mcmc, project = compiled_model)

for (i in 1:n_datasets) {
  
  # simulate data
  y <- simulate_spatial_ranef(
    params = list(
      psi_mean = psi, #(logit scale)
      psi_sd = psi_sd, #(logit scale)
      p = p
    ), 
    nRegions, region,
    nSites, nVisits
  )
  
  # add data
  compiled_model$y <- y
  
  # fit model
  samples <- fit_spatial_ranef(compiled_model, compiled_mcmc, region, nRegions,
                               niter, nburnin, thin)
  
  # calculate p values
  p_values_sp_ranef[i, "deviance"] <- mean(
    samples[, "D_rep_total"] > samples[, "D_obs_total"]
  )
  p_values_sp_ranef[i, "chi_sq"] <- mean(
    samples[, "chi_rep_total"] > samples[, "chi_obs_total"]
  )
  p_values_sp_ranef[i, "lik_ratio"] <- mean(
    samples[, "ratio_rep_total"] > samples[, "ratio_obs_total"]
  )
  p_values_sp_ranef[i, "ftukey"] <- mean(
    samples[, "tukey_rep_total"] > samples[, "tukey_obs_total"]
  )
  
}
saveRDS(p_values_sp_ranef, "pvalues/null/sp_ranef.rds")


###########################################################
# model 4 - spatial ranef, not conditioned on latent states

p_values_sp_ranef <- matrix(NA, nrow = n_datasets, ncol = 4)
colnames(p_values_sp_ranef) <- c("deviance", "chi_sq", "lik_ratio", "ftukey")

# simulate data to build model
region <- rep(1:nRegions, each = nSites / nRegions)
y <- simulate_spatial_ranef(
  params = list(
    psi_mean = psi, #(logit scale)
    psi_sd = psi_sd, #(logit scale)
    p = p
  ), 
  nRegions, region,
  nSites, nVisits
)

# create model instance and compile
constants <- list(nSites = nSites, nVisits = nVisits, 
                  nRegions = nRegions, region = region)
model <- nimbleModel(model_spatial_ranef_latent, constants = constants,
                     data = list(y = y))
compiled_model <- compileNimble(model)

# configure MCMC
mcmc_conf <- configureMCMC(model, 
                           monitors = c("psi_mean", "psi_sd", "p", "z",
                                        "D_obs_total", "D_rep_total",
                                        "chi_obs_total", "chi_rep_total",
                                        "ratio_obs_total", "ratio_rep_total",
                                        "tukey_obs_total", "tukey_rep_total"))

# build and compile MCMC
mcmc <- buildMCMC(mcmc_conf)
compiled_mcmc <- compileNimble(mcmc, project = compiled_model)

for (i in 1:n_datasets) {
  
  # simulate data
  y <- simulate_spatial_ranef(
    params = list(
      psi_mean = psi, #(logit scale)
      psi_sd = psi_sd, #(logit scale)
      p = p
    ), 
    nRegions, region,
    nSites, nVisits
  )
  
  # add data
  compiled_model$y <- y
  
  # fit model
  samples <- fit_spatial_ranef_latent(compiled_model, compiled_mcmc, region, 
                                      nRegions, niter, nburnin, thin)
  
  # calculate p values
  p_values_sp_ranef[i, "deviance"] <- mean(
    samples[, "D_rep_total"] > samples[, "D_obs_total"]
  )
  p_values_sp_ranef[i, "chi_sq"] <- mean(
    samples[, "chi_rep_total"] > samples[, "chi_obs_total"]
  )
  p_values_sp_ranef[i, "lik_ratio"] <- mean(
    samples[, "ratio_rep_total"] > samples[, "ratio_obs_total"]
  )
  p_values_sp_ranef[i, "ftukey"] <- mean(
    samples[, "tukey_rep_total"] > samples[, "tukey_obs_total"]
  )
  
}
saveRDS(p_values_sp_ranef, "pvalues/null/sp_ranef_latent.rds")

