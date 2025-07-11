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
for (n in 1:nDatasets) {
  for (i in 1:length(rho)) {
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
model_uncompiled <- model_uncompiled <- nimbleModel(
  model_minimal, constants = constants,
  data = list(y = simulate_betabinomial(params = list(psi = psi, p = p),
                                        nSites, nVisits, rho[1])))

# param names to monitor in MCMC
MCMC_monitors <- c("psi", "p", "z")

# axis of model breakage
breakage_axis <- rho

# lists of data names, param names, and param indices
data_name_list <- list("y", # conditioned on latent state
                       c("y", "z") # not conditioned on latent state
                  )
param_name_list <- list(c("p", "psi", "z"), # conditioned on latent state
                        c("p", "psi") # not conditioned on latent state
)
param_indices_list <- list(1:(nSites + 2), # conditioned on latent state
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
init_function <- function(y) list(psi = runif(1, 0, 1), p = runif(1, 0, 1), 
                                  z = pmin(y, 1)) 

# run cppp simulations
betabin_out <- run_cppp_simulations(
  constants, simulated_data, model_uncompiled, MCMC_monitors,
  breakage_axis, data_name_list,  param_name_list,  param_indices_list, 
  discrepancyFunctions, discrepancyNames, discrepancyFunctionsArgs,
  coverage_params, init_function, nDatasets, niter, nburnin, thin,
  nCalibrationReplicates
) 

# add new column to output
all_data <- betabin_out %>% 
  mutate(all_param = ifelse(psi & p, TRUE, FALSE))

# get density plots
density_condition <- get_cppp_density_plot(all_data, cdtn = T)
ggsave("figures/occupancy/betabin/density_true.png",
       density_condition, dpi = 400, height = 6, width = 6)

density_nocondition <- get_cppp_density_plot(all_data, cdtn = F)
ggsave("figures/occupancy/betabin/density_false.png",
       density_nocondition, dpi = 400, height = 6, width = 6)

# get dot plots
coverpsi_condition <- get_cppp_dot_plot(all_data, param = "psi", 
                                        breakage_axis_name = "rho", cdtn = T)
ggsave("figures/occupancy/betabin/cover_psi_true.png",
       coverpsi_condition, dpi = 400, height = 6, width = 6)

coverpsi_nocondition <- get_cppp_dot_plot(all_data, param = "psi", 
                                          breakage_axis_name = "rho", cdtn = F)
ggsave("figures/occupancy/betabin/cover_psi_false.png",
       coverpsi_nocondition, dpi = 400, height = 6, width = 6)

coverp_condition <- get_cppp_dot_plot(all_data, param = "p", 
                                      breakage_axis_name = "rho", cdtn = T)
ggsave("figures/occupancy/betabin/cover_p_true.png",
       coverp_condition, dpi = 400, height = 6, width = 6)


coverp_nocondition <- get_cppp_dot_plot(all_data, param = "p", 
                                        breakage_axis_name = "rho", cdtn = F)
ggsave("figures/occupancy/betabin/cover_p_false.png",
       coverp_nocondition, dpi = 400, height = 6, width = 6)

coverall_condition <- get_cppp_dot_plot(all_data, param = "all_param", 
                                        breakage_axis_name = "rho", cdtn = T)
ggsave("figures/occupancy/betabin/cover_all_true.png",
       coverall_condition, dpi = 400, height = 6, width = 6)

coverall_nocondition <- get_cppp_dot_plot(all_data, param = "all_param", 
                                          breakage_axis_name = "rho", cdtn = F)
ggsave("figures/occupancy/betabin/cover_all_false.png",
       coverall_nocondition, dpi = 400, height = 6, width = 6)

# get power plots
power_plot_psi <- get_cppp_power_plot(all_data, param = "psi", 
                                      breakage_axis_name = "rho")
ggsave("figures/occupancy/betabin/power_psi.png",
       power_plot_psi, dpi = 400, height = 6, width = 7)    

power_plot_p <- get_cppp_power_plot(all_data, param = "p", 
                                    breakage_axis_name = "rho")
ggsave("figures/occupancy/betabin/power_p.png",
       power_plot_p, dpi = 400, height = 6, width = 7)   


power_plot_all <- get_cppp_power_plot(all_data, param = "all_param", 
                                      breakage_axis_name = "rho")
ggsave("figures/occupancy/betabin/power_all.png",
       power_plot_all, dpi = 400, height = 6, width = 7)   
