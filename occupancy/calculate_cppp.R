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
simulated_data <- array(NA, dim = c(nDatasets, length(rho), nSites))
for (n in seq_along(1:nDatasets)) {
  for (i in seq_along(rho)) {
    # simulate data
    simulated_data[n, i, ] <- simulate_betabinomial(
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
                                data = list(y = simulated_data[1, 1, ]))

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
init_function <- function(y) {
  list(psi = runif(1, 0, 1), p = runif(1, 0, 1), z = pmin(y, 1)) 
}

# run cppp simulations
betabin_out <- run_cppp_simulations(
  constants, simulated_data, model_uncompiled, mcmc_monitors,
  breakage_axis, data_name_list,  param_name_list,  param_indices_list, 
  discrepancyFunctions, discrepancyNames, discrepancyFunctionsArgs,
  coverage_params, init_function, nDatasets, niter, nburnin, thin,
  nCalibrationReplicates
) 

# add new column to output
all_data <- betabin_out %>% 
  mutate(all_param = ifelse(psi & p, TRUE, FALSE))


########
# plot #
########

# print all plots
get_cppp_plot(plot_type = c("density", "dot", "power"), all_data, 
              param = c("p", "psi", "all_param"), 
              breakage_axis_name = "rho", cdtn = c(TRUE, FALSE), 
              print = TRUE)

# save all plots
get_cppp_plot(plot_type = c("density", "dot", "power"), all_data, 
              param = c("p", "psi", "all_param"), 
              breakage_axis_name = "rho", cdtn = c(TRUE, FALSE), 
              print = FALSE, save = TRUE, 
              filepath = "figures/occupancy/betabin")

# print one plot
get_cppp_plot(plot_type = "dot", all_data, param = "all_param", 
              breakage_axis_name = "rho", cdtn = TRUE, print = TRUE)
