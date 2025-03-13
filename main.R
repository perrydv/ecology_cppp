library(nimble)

# load models
source("occupancy_model.R")

# load functions to simulate data
source("simulate_data.R")

#################
# simulate data #
#################

nSites <- 12
nVisits <- 6
nRegions <- 3
nCov_site <- 2
nCov_detect <- 2
p <- 0.4
psi <- 0.6
psi_sd <- 2.5

# model 1 - basic
y_model1 <- simulate_basic(
  params = list(
    p = p, 
    psi = psi
  ), 
  nSites, nVisits
)


# model 2 - occupancy covariates
x_site <- matrix(NA, nrow = nSites, ncol = nCov_site + 1) # covariate data
x_site[, 1] <- 1
for (i in 1:nCov_site) {
  x_site[, i + 1] <- rnorm(nSites, 0, 1)
}

y_model2 <- simulate_cov_occ(
  params = list(
    beta = runif(nCov_site + 1, -1, 1),
    p = 0.4
  ), 
  x_site, nSites, nVisits
)


# model 3 - occupancy and detection covariates
x_visit <- array(NA, dim = c(nSites, nVisits, nCov_detect + 1)) # covariate data
x_visit[, , 1] <- 1
for (i in 1:nCov_detect) {
  x_visit[, , i + 1] <- matrix(rnorm(nSites * nVisits, 0, 1), nrow = nSites)
}

y_model3 <- simulate_cov_occ_det(
  params = list(
    beta_o = runif(nCov_site + 1, -1, 1),
    beta_d = runif(nCov_detect + 1, -1, 1)
  ), 
  x_site, x_visit,
  nSites, nVisits
)


# model 4 - spatial random effect
region <- rep(1:nRegions, each = nSites / nRegions)

y_model4 <- simulate_spatial_ranef(
  params = list(
    psi_mean = psi, #(logit scale)
    psi_sd = psi_sd, #(logit scale)
    p = p
  ), 
  nRegions, region,
  nSites, nVisits
)



