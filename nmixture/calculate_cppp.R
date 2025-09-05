library(nimble)
library(tidyverse)
library(patchwork)

# load CPPP functions
source("cppp/sally_code/calculateCPPP_original.R")

# load functions to simulate data
source("nmixture/simulate_data.R")

# load minimal model
source("nmixture/minimal_model.R")

# load discrepancy functions
source("nmixture/nmixture_discrepancy_functions.R")

# load auxiliary functions
source("utils.R")


###########################
# beta-binomial detection #
###########################

# data size
nSites <- 200
nVisits <- 6
nDatasets <- 100

# parameter values
p <- 0.3
rho <- c(0.01, 0.1, 0.25, 0.5)
lambda <- 100

# MCMC
niter <- 5000
nburnin <- 1000
thin <- 5
nCalibrationReplicates <- 100

# simulate data
simulated_y <- array(NA, dim = c(nDatasets, length(rho), nSites, nVisits))
for (n in seq_along(1:nDatasets)) {
  for (i in seq_along(rho)) {
    # simulate data
    simulated_y[n, i, , ] <- simulate_nmixture_betabin(
      params = list(lambda = lambda, p = p),
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
model_uncompiled <- nimbleModel(model_minimal_nmixture, constants = constants,
                                data = list(y = simulated_y[1, 1, , ]))

# param names to monitor in MCMC
mcmc_monitors <- c("lambda", "p", "N")

# axis of model breakage
breakage_axis <- rho

# lists of data names, param names, and param indices
data_name_list <- list(
  "y", # conditioned on latent state
  c("y", "N") # not conditioned on latent state
)
param_name_list <- list(
  c("p", "lambda", "N"), # conditioned on latent state
  c("p", "lambda") # not conditioned on latent state
)
# param_indices_list <- list(
#   1:(nSites + 2), # conditioned on latent state
#   1:2 # not conditioned on latent state
# )

# discrepancy functions and arguments
discrepancyFunctions <- list(ratioDiscFunction_nmix, devianceDiscFunction_nmix,
                             chisqDiscFunction_nmix, tukeyDiscFunction_nmix)
discrepancyNames <- c("Ratio", "Deviance", "Chi-Square", "Freeman-Tukey")
args <- list(nVisits = nVisits, dataNames = "y",
             latent_N = "N", prob_detection = "p")
discrepancyFunctionsArgs <- list(args, args, args, args)

# named list of param calculating coverage: values = true param values
coverage_params <- c(lambda, p)
names(coverage_params) <- c("lambda", "p")

# function to generate initial values for MCMC
init_function_nmix_betabin <- function(simulated_data, args) {
  list(lambda = runif(1, 10, 1000), p = runif(1, 0, 1),
       N = round(rowMeans(simulated_data[[1]])) * 2)
}

# run cppp simulations
betabin_out <- run_cppp_simulations(
  constants, simulated_data = list(simulated_y),
  sim_data_names = "y", model_uncompiled, mcmc_monitors,
  breakage_axis, data_name_list,  param_name_list,  #param_indices_list,
  discrepancyFunctions, discrepancyNames, discrepancyFunctionsArgs,
  coverage_params, init_function = init_function_nmix_betabin,
  init_args = list(),
  nDatasets, niter, nburnin, thin,
  nCalibrationReplicates,
  condition_on_latent_states = c(TRUE, FALSE)
)

# add new column to output
all_data <- betabin_out %>%
  mutate(all_param = ifelse(lambda & p, TRUE, FALSE))
saveRDS(all_data, "nmixture/saved_outputs/output_betabin.rds")


