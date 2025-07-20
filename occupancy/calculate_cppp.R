library(nimble)
library(tidyverse)
library(patchwork)

# load CPPP functions
source("cppp/sally_code/calculateCPPP_original.R")

# load functions to simulate data
source("occupancy/simulate_data.R")

# load minimal model
source("occupancy/minimal_model.R")

# load discrepancy functions
source("occupancy/occupancy_discrepancy_functions.R")

# load auxiliary functions
source("utils.R")


#################
# beta-binomial #
#################

# data size
nSites <- 50
nVisits <- 6
nDatasets <- 100

# parameter values
p <- 0.3
rho <- c(0.01, 0.1, 0.25, 0.5)
psi <- 0.6

# MCMC
niter <- 5000
nburnin <- 1000
thin <- 5
nCalibrationReplicates <- 100

# simulate data
simulated_y <- array(NA, dim = c(nDatasets, length(rho), nSites))
for (n in seq_along(1:nDatasets)) {
  for (i in seq_along(rho)) {
    # simulate data
    simulated_y[n, i, ] <- simulate_betabinomial(
      params = list(psi = psi, p = p),
      nSites, nVisits, rho[i]
    )
  }
}

##
# cppp function inputs
##
# model constants
constants <- list(nSites = nSites, nVisits = nVisits)

# uncompiled model - temporarily add data so that the MCMC samplers get set up
model_uncompiled <- nimbleModel(model_minimal, constants = constants,
                                data = list(y = simulated_y[1, 1, ]))

# param names to monitor in MCMC
mcmc_monitors <- c("psi", "p", "z")

# axis of model breakage
breakage_axis <- rho

# lists of data names, param names, and param indices
data_name_list <- list(
  "y", # conditioned on latent state
  c("y", "z") # not conditioned on latent state
)
param_name_list <- list(
  c("p", "psi", "z"), # conditioned on latent state
  c("p", "psi") # not conditioned on latent state
)
param_indices_list <- list(
  1:(nSites + 2), # conditioned on latent state
  1:2 # not conditioned on latent state
)

# discrepancy functions and arguments
discrepancyFunctions <- list(ratioDiscFunction, devianceDiscFunction,
                             chisqDiscFunction, tukeyDiscFunction)
discrepancyNames <- c("Ratio", "Deviance", "Chi-Square", "Freeman-Tukey")
args <- list(nVisits = nVisits, dataNames = "y",
             latent_occ = "z", prob_detection = "p", prob_occupancy = "psi")
discrepancyFunctionsArgs <- list(args, args, args, args)

# named list of param calculating coverage: values = true param values
coverage_params <- c(psi, p)
names(coverage_params) <- c("psi", "p")

# function to generate initial values for MCMC
init_function_betabin <- function(simulated_data, args) {
  list(psi = runif(1, 0, 1), p = runif(1, 0, 1),
       z = pmin(simulated_data[[1]], 1))
}

# run cppp simulations
betabin_out <- run_cppp_simulations(
  constants, simulated_data = list(simulated_y),
  sim_data_names = "y", model_uncompiled, mcmc_monitors,
  breakage_axis, data_name_list,  param_name_list,  param_indices_list,
  discrepancyFunctions, discrepancyNames, discrepancyFunctionsArgs,
  coverage_params, init_function = init_function_betabin,
  init_args = list(),
  nDatasets, niter, nburnin, thin,
  nCalibrationReplicates,
  condition_on_latent_states = c(TRUE, FALSE)
)

# add new column to output
all_data <- betabin_out %>%
  mutate(all_param = ifelse(psi & p, TRUE, FALSE))
saveRDS(all_data, "occupancy/saved_outputs/output_betabin.rds")


########################################
# site occupancy covariate interaction #
########################################

# data size
nSites <- 50
nVisits <- 6
nDatasets <- 100

# parameter values
p <- 0.1
beta <- c(0, 0.5)
beta2 <- c(0, 1, 3, 5)

# MCMC 
niter <- 5000
nburnin <- 1000
thin <- 5
nCalibrationReplicates <- 100

# simulate data
simulated_y <- array(NA, dim = c(nDatasets, length(beta2), nSites, nVisits))
simulated_x <- array(NA, dim = c(nDatasets, length(beta2), nSites, 
                                 length(beta) + 1))

for (n in seq_along(1:nDatasets)) {
  for (i in seq_along(beta2)) {
    
    # simulate X
    simulated_x[n, i, , ] <- cbind(rep(1, nSites), rnorm(nSites), rnorm(nSites))
    
    # simulate y
    simulated_y[n, i, , ] <- simulate_cov_occ_inter(
      params = list(p = p, beta = c(beta, beta2[i])), 
      simulated_x[n, i, , ], nSites, nVisits
    )
  }
}

##
# cppp function inputs
## 
# model constants
constants <- list(nSites = nSites, nVisits = nVisits, ncov = length(beta))

# uncompiled model - temporarily add data so that the MCMC samplers get set up
model_uncompiled <- nimbleModel(model_cov_occ, constants = constants,
                                data = list(y = simulated_y[1, 1, , ],
                                            x_site = simulated_x[1, 1, , 1:2]))

# param names to monitor in MCMC
mcmc_monitors <- c("psi", "p", "z", "beta")

# axis of model breakage
breakage_axis <- beta2

# lists of data names, param names, and param indices
data_name_list <- list(
  c("y", "z") # not conditioned on latent state
)
param_name_list <- list(
  c("p", "beta") # not conditioned on latent state
)
param_indices_list <- list(
  1:(length(beta) + 1) # not conditioned on latent state
)

# discrepancy functions and arguments
discrepancyFunctions <- list(chisqDiscFunction_z_cov, tukeyDiscFunction_z_cov)
discrepancyNames <- c("Chi-Square", "Freeman-Tukey")
args <- list(y = "y", x_site = "x_site", nVisits = nVisits,
             beta = "beta", latent_occ = "z")
discrepancyFunctionsArgs <- list(args, args)

# named list of param calculating coverage: values = true param values
coverage_params <- c(beta, p)
names(coverage_params) <- c("beta[1]", "beta[2]", "p")

# function to generate initial values for MCMC
init_function_sitecov <- function(simulated_data, args) {
  list(beta = rnorm(args$ncov, 0, 1), p = runif(1, 0, 1), 
       z = pmin(rowSums(simulated_data[[1]]), 1))  
}

# run cppp simulations
sitecov_out <- run_cppp_simulations(
  constants, simulated_data = list(simulated_y, simulated_x[, , , 1:2]), 
  sim_data_names = c("y", "x_site"), model_uncompiled, mcmc_monitors,
  breakage_axis, data_name_list,  param_name_list,  param_indices_list, 
  discrepancyFunctions, discrepancyNames, discrepancyFunctionsArgs,
  coverage_params, init_function = init_function_sitecov,
  init_args = list(ncov = constants$ncov),
  nDatasets, niter, nburnin, thin,
  nCalibrationReplicates,
  condition_on_latent_states = FALSE
) 

# add new column to output
# all_data <- betabin_out %>% 
#   mutate(all_param = ifelse(psi_coverage & p_coverage, TRUE, FALSE))
saveRDS(sitecov_out, "occupancy/saved_outputs/output_sitecov.rds")


