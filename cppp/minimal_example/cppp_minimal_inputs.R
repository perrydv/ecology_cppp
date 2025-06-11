library(nimble)

# load functions to simulate data
source("simulate_data.R")

# load minimal model
source("cppp/minimal_example/cppp_minimal_model.R")

################################################
# simulate data and generate posterior samples #
################################################

# data size
nSites <- 50
nVisits <- 6

# parameter values
p <- 0.3
rho <- 0.5
psi <- 0.6

# MCMC 
niter <- 5000
nburnin <- 1000
thin <- 5

# simulate data
set.seed(123)
y <- simulate_betabinomial(params = list(psi = psi, p = p),
                           nSites, nVisits, rho)

# create model instance and compile
constants <- list(nSites = nSites, nVisits = nVisits)
model_uncompiled <- nimbleModel(model_minimal, constants = constants,
                                data = list(y = y))
model <- compileNimble(model_uncompiled)

# configure MCMC
mcmc_conf <- configureMCMC(model_uncompiled, monitors = c("psi", "p", "z"))

# build MCMC
mcmc <- buildMCMC(mcmc_conf)

# compile mcmc
compiled_mcmc <- compileNimble(mcmc, project = model)

inits <- function(y) list(psi = runif(1, 0, 1), p = runif(1, 0, 1), 
                          z = pmin(y, 1)) 

# add initial values
model$setInits(inits(model$y))

# generate samples
MCMCOutput <- runMCMC(compiled_mcmc, niter = niter, 
                      nburnin = nburnin, thin = thin)


# save data, samples, and compiled model
# saveRDS(model, "cppp/minimal_example/saved_inputs/model.rds")
# saveRDS(MCMCOutput, "cppp/minimal_example/saved_inputs/samples.rds")
# saveRDS(y, "cppp/minimal_example/saved_inputs/y.rds")


#########################
# discrepancy functions #
#########################

## discrepancy base class
discrepancyFunction_BASE <- nimbleFunctionVirtual(
  run = function() returnType(double())
)

## chi-squared discrepancy function
chisqDiscFunction <- nimbleFunction(
  contains = discrepancyFunction_BASE,
  setup = function(discrepancyFunctionsArgs){
    nVisits <- discrepancyFunctionsArgs[["nVisits"]]
  },
  run = function(compiled_model){
    
    # get y_exp
    y_exp <- compiled_model$z * compiled_model$p * nVisits
    
    # calculate chi squared discrepancy measure
    chi_out <- (compiled_model$y - y_exp) ^ 2 / (compiled_model$y + 1e-6)
    
    returnType(double(0)) 
    return(sum(chi_out))
  }
)


###################
# function inputs #
###################

# if *not* conditioning on latent state
paramNames <- c("p", "psi")
simNodes <- c("z", "y")

# if conditioning on latent state
paramNames <- c("p", "psi", "z")
simNodes <- "y"

discrepancyFunctions <- list(chisqDiscFunction)
discrepancyFunctionsArgs <- list(list(nVisits = nVisits))

