library(nimble)

# load models
source("occupancy_model.R")

# load functions to simulate data
source("simulate_data.R")

# load auxiliary functions
source("utils.R")

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


#################################
# non-independence across sites #
#################################

# parameter values
p <- 0.3
psi <- 0.6
psi_sd <- c(2.5, 5, 10)

# create df to store samples
samples_p <- matrix(NA, nrow = n_datasets * length(psi_sd), 
                    ncol = (niter - nburnin) / thin + 2)
colnames(samples_p) <- c("dataset", "sd", 1:((niter - nburnin) / thin))
samples_psi <- matrix(NA, nrow = n_datasets * length(psi_sd), 
                      ncol = (niter - nburnin) / thin + 2)
colnames(samples_psi) <- c("dataset", "sd", 1:((niter - nburnin) / thin))

# create df to store p values
pvalues_nonind <- matrix(NA, nrow = n_datasets * length(psi_sd), ncol = 9)
colnames(pvalues_nonind) <- c("sd", "deviance", "chi_sq", "lik_ratio", "ftukey",
                              "deviance_latent", "chi_sq_latent", 
                              "lik_ratio_latent", "ftukey_latent")

# create model instance and compile
constants <- list(nSites = nSites, nVisits = nVisits)
model <- nimbleModel(model_basic, constants = constants, 
                     data = list(y = simulate_basic(
                       params = list(p = p, psi = psi), nSites, nVisits)))
compiled_model <- compileNimble(model)

# configure MCMC
mcmc_conf <- configureMCMC(
  model, 
  monitors = c("psi", "p", "z",
               "D_obs_total", "D_rep_total", "D_rep_latent_total",
               "chi_obs_total", "chi_rep_total", "chi_rep_latent_total",
               "ratio_obs_total", "ratio_rep_total", "ratio_rep_latent_total",
               "tukey_obs_total", "tukey_rep_total", "tukey_rep_latent_total")
)

# build MCMC
mcmc <- buildMCMC(mcmc_conf)

# compile mcmc
compiled_mcmc <- compileNimble(mcmc, project = compiled_model)

for (j in 1:length(psi_sd)) {
  for (i in 1:n_datasets) {
    
    # simulate data
    region <- rep(1:nRegions, each = nSites / nRegions)
    y <- simulate_spatial_ranef(
      params = list(
        psi_mean = psi, #(logit scale)
        psi_sd = psi_sd[j], #(logit scale)
        p = p
      ), 
      nRegions, region,
      nSites, nVisits
    )
    compiled_model$y <- y
    
    # fit model
    samples <- fit_basic(compiled_model, compiled_mcmc,
                         niter, nburnin, thin)
    
    # get p-value df index
    index <- (j - 1) * n_datasets + i
    
    # save samples
    samples_psi[index, ] <- c(i, psi_sd[j], samples[, "psi"])
    samples_p[index, ] <- c(i, psi_sd[j], samples[, "p"])
    
    # save sd
    pvalues_nonind[index, "sd"] <- psi_sd[j]
    
    # calculate p values
    pvalues_nonind[index, c("deviance", "chi_sq", "lik_ratio", 
                             "ftukey", "deviance_latent", 
                             "chi_sq_latent", 
                             "lik_ratio_latent", 
                             "ftukey_latent")] <- calc_pvalues(n_measures = 8,
                                                               samples)
    
  }
}

saveRDS(pvalues_nonind, "pvalues/alt/nonind.rds")
saveRDS(samples_p, "posterior/alt/p_nonind.rds")
saveRDS(samples_psi, "posterior/alt/psi_nonind.rds")


###########################
# beta-binomial detection #
###########################

# parameter values
p <- 0.3
rho <- c(0.1, 0.5, 0.9)
psi <- 0.6

# create df to store p values
pvalues_betabin <- matrix(NA, nrow = n_datasets * length(rho), ncol = 9)
colnames(pvalues_betabin) <- c("rho", "deviance", "chi_sq", "lik_ratio", 
                               "ftukey", "deviance_latent", "chi_sq_latent", 
                               "lik_ratio_latent", "ftukey_latent")

# create df to store samples
samples_p <- matrix(NA, nrow = n_datasets * length(rho), 
                    ncol = (niter - nburnin) / thin + 2)
colnames(samples_p) <- c("dataset", "rho", 1:((niter - nburnin) / thin))
samples_psi <- matrix(NA, nrow = n_datasets * length(rho), 
                      ncol = (niter - nburnin) / thin + 2)
colnames(samples_psi) <- c("dataset", "rho", 1:((niter - nburnin) / thin))

# create model instance and compile
constants <- list(nSites = nSites, nVisits = nVisits)
model <- nimbleModel(model_basic_bin, constants = constants, 
                     data = list(y = simulate_betabinomial(
                       params = list(p = p, psi = psi), nSites, nVisits,
                       rho = 0.01)))
compiled_model <- compileNimble(model)

# configure MCMC
mcmc_conf <- configureMCMC(
  model, 
  monitors = c("psi", "p", "z",
               "D_obs_total", "D_rep_total", "D_rep_latent_total",
               "chi_obs_total", "chi_rep_total", "chi_rep_latent_total",
               "ratio_obs_total", "ratio_rep_total", "ratio_rep_latent_total",
               "tukey_obs_total", "tukey_rep_total", "tukey_rep_latent_total")
)

# build MCMC
mcmc <- buildMCMC(mcmc_conf)

# compile mcmc
compiled_mcmc <- compileNimble(mcmc, project = compiled_model)

for (j in 1:length(rho)) {
  for (i in 1:n_datasets) {
    
    # simulate data
    y <- simulate_betabinomial(
      params = list(
        psi = psi,
        p = p
      ), 
      nSites, nVisits, rho[j]
    )
    compiled_model$y <- y
    
    # fit model
    samples <- fit_basic_betabin(compiled_model, compiled_mcmc,
                                 niter, nburnin, thin)
    
    # get p-value df index
    index <- (j - 1) * n_datasets + i
    
    # save samples
    samples_psi[index, ] <- c(i, rho[j], samples[, "psi"])
    samples_p[index, ] <- c(i, rho[j], samples[, "p"])
    
    # save rho
    pvalues_betabin[index, "rho"] <- rho[j]
    
    # calculate p values
    pvalues_betabin[index, c("deviance", "chi_sq", "lik_ratio", 
                             "ftukey", "deviance_latent", "chi_sq_latent", 
                             "lik_ratio_latent", 
                             "ftukey_latent")] <- calc_pvalues(n_measures = 8,
                                                               samples)
    
  }
}

saveRDS(pvalues_betabin, "pvalues/alt/betabinom.rds")
saveRDS(samples_p, "posterior/alt/p_betabinom.rds")
saveRDS(samples_psi, "posterior/alt/psi_betabinom.rds")


#####################
# mixture detection #
#####################

# parameter values
p1 <- 0.2
p2 <- 0.8
pMix <- c(0.1, 0.2, 0.5)
psi <- 0.6

# create df to store p values
pvalues_detMix <- matrix(NA, nrow = n_datasets * length(pMix), ncol = 9)
colnames(pvalues_detMix) <- c("pMix", "deviance", "chi_sq", "lik_ratio", 
                              "ftukey", "deviance_latent", "chi_sq_latent", 
                              "lik_ratio_latent", "ftukey_latent")

# create df to store samples
samples_p <- matrix(NA, nrow = n_datasets * length(pMix), 
                    ncol = (niter - nburnin) / thin + 2)
colnames(samples_p) <- c("dataset", "pMix", 1:((niter - nburnin) / thin))
samples_psi <- matrix(NA, nrow = n_datasets * length(pMix), 
                      ncol = (niter - nburnin) / thin + 2)
colnames(samples_psi) <- c("dataset", "pMix", 1:((niter - nburnin) / thin))

# create model instance and compile
constants <- list(nSites = nSites, nVisits = nVisits)
model <- nimbleModel(model_basic, constants = constants, 
                     data = list(y = simulate_det_pMix(
                       params = list(p1 = p1, p2 = p2, psi = psi), 
                       nSites, nVisits, pMix = 0.1)))
compiled_model <- compileNimble(model)

# configure MCMC
mcmc_conf <- configureMCMC(
  model, 
  monitors = c("psi", "p", "z",
               "D_obs_total", "D_rep_total", "D_rep_latent_total",
               "chi_obs_total", "chi_rep_total", "chi_rep_latent_total",
               "ratio_obs_total", "ratio_rep_total", "ratio_rep_latent_total",
               "tukey_obs_total", "tukey_rep_total", "tukey_rep_latent_total")
)

# build MCMC
mcmc <- buildMCMC(mcmc_conf)

# compile mcmc
compiled_mcmc <- compileNimble(mcmc, project = compiled_model)

for (j in 1:length(pMix)) {
  for (i in 1:n_datasets) {
    
    # simulate data
    y <- simulate_det_pMix(
      params = list(p1 = p1, p2 = p2, psi = psi), 
      nSites, nVisits, pMix[j]
    )
    compiled_model$y <- y
    
    # fit model
    samples <- fit_basic(compiled_model, compiled_mcmc,
                                 niter, nburnin, thin)
    
    # get p-value df index
    index <- (j - 1) * n_datasets + i
    
    # save samples
    samples_psi[index, ] <- c(i, pMix[j], samples[, "psi"])
    samples_p[index, ] <- c(i, pMix[j], samples[, "p"])
    
    # save pMix
    pvalues_detMix[index, "pMix"] <- pMix[j]
    
    # calculate p values 
    pvalues_detMix[index, c("deviance", "chi_sq", "lik_ratio", 
                            "ftukey", "deviance_latent", "chi_sq_latent", 
                            "lik_ratio_latent", 
                            "ftukey_latent")] <- calc_pvalues(n_measures = 8,
                                                              samples)
  }
}

saveRDS(pvalues_detMix, "pvalues/alt/detectionMix.rds")
saveRDS(samples_p, "posterior/alt/p_detectionMix.rds")
saveRDS(samples_psi, "posterior/alt/psi_detectionMix.rds")


