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

# # data size
# nSites <- 100
# nVisits <- 6
# nDatasets <- 100
# 
# # parameter values
# p <- 0.3
# rho <- c(0.00001, 0.0005, 0.001, 0.005)
# lambda <- 100
# mu <- lambda * p
# 
# # MCMC
# niter <- 50000
# nburnin <- 2000
# thin <- 10
# nCalibrationReplicates <- 100
# 
# # simulate data
# simulated_y <- array(NA, dim = c(nDatasets, length(rho), nSites, nVisits))
# for (n in seq_along(1:nDatasets)) {
#   for (i in seq_along(rho)) {
#     # simulate data
#     simulated_y[n, i, , ] <- simulate_nmixture_betabin(
#       params = list(lambda = lambda, p = p),
#       nSites, nVisits, rho[i]
#     )
#   }
# }
# 
# ##
# # cppp function inputs
# ##
# # model constants
# constants <- list(nSites = nSites, nVisits = nVisits)
# 
# # uncompiled model - temporarily add data so that the MCMC samplers get set up
# model_uncompiled <- nimbleModel(model_minimal_nmixture, constants = constants,
#                                 data = list(y = simulated_y[1, 1, , ]))
# 
# # param names to monitor in MCMC
# mcmc_monitors <- c("mu", "p", "N")
# 
# # axis of model breakage
# breakage_axis <- rho
# 
# # lists of data names, param names, and param indices
# data_name_list <- list(
#   "y", # conditioned on latent state
#   c("y", "N") # not conditioned on latent state
# )
# param_name_list <- list(
#   c("p", "mu", "N"), # conditioned on latent state
#   c("p", "mu") # not conditioned on latent state
# )
# 
# # discrepancy functions and arguments
# discrepancyFunctions <- list(linkDiscFunction_nmix)
# discrepancyNames <- c("Link")
# args <- list(nVisits = nVisits, dataNames = "y",
#              latent_N = "N", prob_detection = "p")
# # discrepancyFunctionsArgs <- list(args, args, args, args)
# discrepancyFunctionsArgs <- list(args)
# 
# # named list of param calculating coverage: values = true param values
# coverage_params <- c(mu, p)
# names(coverage_params) <- c("mu", "p")
# 
# # function to generate initial values for MCMC
# init_function_nmix <- function(simulated_data, args) {
#   list(mu = runif(1, 10, 1000), p = runif(1, 0, 1),
#        N = round(apply(simulated_data[[1]], 1, median)) * 15)
# }
# 
# # run cppp simulations
# betabin_out <- run_cppp_simulations(
#   constants, simulated_data = list(simulated_y),
#   sim_data_names = "y", model_uncompiled, mcmc_monitors,
#   breakage_axis, data_name_list,  param_name_list,  #param_indices_list,
#   discrepancyFunctions, discrepancyNames, discrepancyFunctionsArgs,
#   coverage_params, init_function = init_function_nmix,
#   init_args = list(),
#   nDatasets, niter, nburnin, thin, nchain,
#   nCalibrationReplicates,
#   condition_on_latent_states = c(TRUE, FALSE)
# )
# 
# # add new column to output
# all_data <- betabin_out %>%
#   mutate(all_param = ifelse(mu & p, TRUE, FALSE))
# saveRDS(all_data, "nmixture/saved_outputs/output_betabin_20250930.rds")


#########################
# unmodeled N variation #
#########################

# # data size
# nSites <- 100
# nVisits <- 6
# nDatasets <- 100
# 
# # parameter values
# p <- 0.3
# rho <- c(1, 0.7, 0.4, 0.1)
# lambda <- 100
# mu <- lambda * p
# 
# # MCMC
# niter <- 50000
# nburnin <- 2000
# thin <- 10
# nCalibrationReplicates <- 100
# 
# # simulate data
# simulated_y <- array(NA, dim = c(nDatasets, length(rho), nSites, nVisits))
# for (n in seq_along(1:nDatasets)) {
#   for (i in seq_along(rho)) {
#     # simulate data
#     simulated_y[n, i, , ] <- simulate_nmixture_Nvar(
#       params = list(lambda = lambda, p = p),
#       nSites, nVisits, rho[i]
#     )
#   }
# }
# 
# ##
# # cppp function inputs
# ##
# # model constants
# constants <- list(nSites = nSites, nVisits = nVisits)
# 
# # uncompiled model - temporarily add data so that the MCMC samplers get set up
# model_uncompiled <- nimbleModel(model_minimal_nmixture, constants = constants,
#                                 data = list(y = simulated_y[1, 1, , ]))
# 
# # param names to monitor in MCMC
# mcmc_monitors <- c("mu", "p", "N")
# 
# # axis of model breakage
# breakage_axis <- rho
# 
# # lists of data names, param names, and param indices
# data_name_list <- list(
#   "y", # conditioned on latent state
#   c("y", "N") # not conditioned on latent state
# )
# param_name_list <- list(
#   c("p", "mu", "N"), # conditioned on latent state
#   c("p", "mu") # not conditioned on latent state
# )
# 
# # discrepancy functions and arguments
# discrepancyFunctions <- list(linkDiscFunction_nmix)
# discrepancyNames <- c("Link")
# args <- list(nVisits = nVisits, dataNames = "y",
#              latent_N = "N", prob_detection = "p")
# # discrepancyFunctionsArgs <- list(args, args, args, args)
# discrepancyFunctionsArgs <- list(args)
# 
# # named list of param calculating coverage: values = true param values
# coverage_params <- c(mu, p)
# names(coverage_params) <- c("mu", "p")
# 
# # function to generate initial values for MCMC
# init_function_nmix <- function(simulated_data, args) {
#   list(mu = runif(1, 10, 1000), p = runif(1, 0, 1),
#        N = round(apply(simulated_data[[1]], 1, median)) * 15)
# }
# 
# # run cppp simulations
# unmodN_out <- run_cppp_simulations(
#   constants, simulated_data = list(simulated_y),
#   sim_data_names = "y", model_uncompiled, mcmc_monitors,
#   breakage_axis, data_name_list,  param_name_list,  #param_indices_list,
#   discrepancyFunctions, discrepancyNames, discrepancyFunctionsArgs,
#   coverage_params, init_function = init_function_nmix,
#   init_args = list(),
#   nDatasets, niter, nburnin, thin, nchain
#   nCalibrationReplicates,
#   condition_on_latent_states = c(TRUE, FALSE)
# )
# 
# # add new column to output
# all_data <- unmodN_out %>%
#   mutate(all_param = ifelse(mu & p, TRUE, FALSE))
# saveRDS(all_data, "nmixture/saved_outputs/output_unmodN.rds")


###################
# double counting #
###################

# data size
nSites <- 100
nVisits <- 6
nDatasets <- 20

# parameter values
p <- 0.3
# gamma <- c(0, 0.05, 0.1, 0.15)
gamma <- 0
lambda <- 100
mu <- lambda * p

# MCMC
nchain <- 4
niter <- 50000
nburnin <- 2000
thin <- 10
nCalibrationReplicates <- 500

# simulate data
simulated_y <- array(NA, dim = c(nDatasets, length(gamma), nSites, nVisits))
for (n in seq_along(1:nDatasets)) {
  for (i in seq_along(gamma)) {
    # simulate data
    simulated_y[n, i, , ] <- simulate_nmixture_dcount(
      params = list(lambda = lambda, p = p),
      nSites, nVisits, gamma[i]
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
mcmc_monitors <- c("mu", "p", "N")

# axis of model breakage
breakage_axis <- gamma

# lists of data names, param names, and param indices
data_name_list <- list(
  "y", # conditioned on latent state
  c("y", "N") # not conditioned on latent state
)
param_name_list <- list(
  c("p", "mu", "N"), # conditioned on latent state
  c("p", "mu") # not conditioned on latent state
)

# discrepancy functions and arguments
discrepancyFunctions <- list(linkDiscFunction_nmix)
discrepancyNames <- c("Link")
args <- list(nVisits = nVisits, dataNames = "y",
             latent_N = "N", prob_detection = "p")
# discrepancyFunctionsArgs <- list(args, args, args, args)
discrepancyFunctionsArgs <- list(args)

# named list of param calculating coverage: values = true param values
coverage_params <- c(mu, p)
names(coverage_params) <- c("mu", "p")

# function to generate initial values for MCMC
init_function_nmix <- function(simulated_data, args) {
  list(mu = runif(1, 10, 1000), p = runif(1, 0, 1),
       N = round(apply(simulated_data[[1]], 1, median)) * 15)
}

# run cppp simulations
dcount_out <- run_cppp_simulations(
  constants, simulated_data = list(simulated_y),
  sim_data_names = "y", model_uncompiled, mcmc_monitors,
  breakage_axis, data_name_list,  param_name_list,  #param_indices_list,
  discrepancyFunctions, discrepancyNames, discrepancyFunctionsArgs,
  coverage_params, init_function = init_function_nmix,
  init_args = list(),
  nDatasets, niter, nburnin, thin, nchain,
  nCalibrationReplicates,
  # condition_on_latent_states = c(TRUE, FALSE),
  condition_on_latent_states = TRUE,
  MCMCcontrol = list(niter = niter,
                     thin = thin,
                     nburnin = nburnin)
)

# add new column to output
all_data <- dcount_out %>%
  mutate(all_param = ifelse(mu & p, TRUE, FALSE))
saveRDS(all_data, "nmixture/saved_outputs/output_dcount_20251006_5.rds")


