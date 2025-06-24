library(nimble)

# load CPPP functions
source("cppp/sally_code/calculateCPPP.R")

# load functions to simulate data
source("simulate_data.R")

# load minimal model
source("cppp/minimal_example/cppp_minimal_model.R")

# load discrepancy functions
source("cppp/minimal_example/cppp_discrepancy_functions.R")

# simulation constants
# data size
nSites <- 50
nVisits <- 6

# parameter values
p <- 0.3
rho <- c(0.01, 0.1, 0.25, 0.5, 0.75, 0.9)
psi <- 0.6

# MCMC 
niter <- 5000
nburnin <- 1000
thin <- 5
nCalibrationReplicates <- 100

# discrepancy function list
discrepancyFunctions <- list(ratioDiscFunction, tukeyDiscFunction, 
                             devianceDiscFunction)
discrepancyFunctionsArgs <- list(
  list(nVisits = nVisits),
  list(nVisits = nVisits),
  list(nVisits = nVisits))

############
# get cppp #
############

# create empty list of lists
cppp_out <- lapply(1:length(rho), function(x) list())

for (i in 1:length(rho)) {
  
  # simulate data
  set.seed(123)
  y <- simulate_betabinomial(params = list(psi = psi, p = p),
                             nSites, nVisits, rho[i])
  
  # create model instance and compile
  constants <- list(nSites = nSites, nVisits = nVisits)
  model_uncompiled <- nimbleModel(model_minimal, constants = constants,
                                  data = list(y = y))
  
  inits <- function(y) list(psi = runif(1, 0, 1), p = runif(1, 0, 1), 
                            z = pmin(y, 1)) 
  
  # add initial values
  model_uncompiled$setInits(inits(model_uncompiled$y))
  
  model <- compileNimble(model_uncompiled)
  
  # configure MCMC
  mcmc_conf <- configureMCMC(model_uncompiled, monitors = c("psi", "p", "z"))
  
  # build MCMC
  mcmc <- buildMCMC(mcmc_conf)
  
  # compile mcmc
  compiled_mcmc <- compileNimble(mcmc, project = model)
  
  # generate samples
  MCMCOutput <- runMCMC(compiled_mcmc, niter = niter, 
                        nburnin = nburnin, thin = thin)
  
  ###############
  # cppp inputs #
  ###############
  
  dataNames <- "y"
  
  condition_options <- c(TRUE, FALSE)
  
  for (j in 1:length(condition_options)) {
    if (condition_options[j]) {
      # if conditioning on latent state
      paramNames <- c("p", "psi", "z")
      simNodes <- "y"
    } else {
      # if *not* conditioning on latent state
      paramNames <- c("p", "psi")
      simNodes <- c("z", "y")
    }
    
    # run calibration
    cppp_out[[i]][[j]] <- runCalibration(
      model = model, dataNames = dataNames, paramNames = paramNames, 
      origMCMCSamples = MCMCOutput, discrepancyFunctions = discrepancyFunctions, 
      discrepancyFunctionsArgs = discrepancyFunctionsArgs,
      nCalibrationReplicates = nCalibrationReplicates,
      returnSamples = FALSE, returnDiscrepancies = FALSE) 
  }
  names(cppp_out[[i]]) <- c("true", "false")
}


