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
model <- nimbleModel(model_minimal, constants = constants,
                     data = list(y = y))
compiled_model <- compileNimble(model)

# configure MCMC
mcmc_conf <- configureMCMC(model, monitors = c("psi", "p", "z"))

# build MCMC
mcmc <- buildMCMC(mcmc_conf)

# compile mcmc
compiled_mcmc <- compileNimble(mcmc, project = compiled_model)

inits <- function(y) list(psi = runif(1, 0, 1), p = runif(1, 0, 1), 
                          z = pmin(y, 1)) 

# add initial values
compiled_model$setInits(inits(compiled_model$y))

# generate samples
samples <- runMCMC(compiled_mcmc, niter = niter, 
                   nburnin = nburnin, thin = thin)


# save samples and compiled model
saveRDS(compiled_model, "cppp/minimal_example/saved_inputs/compiled_model.rds")
saveRDS(samples, "cppp/minimal_example/saved_inputs/samples.rds")


#########################
# discrepancy functions #
#########################

## discrepancy base class
discrepancyFunction_BASE <- nimbleFunctionVirtual(
  run = function() returnType(double())
)

## Data log-likelihood discrepancy function
logLikDiscFunction <- nimbleFunction(
  contains = discrepancyFunction_BASE,
  setup = function(model, discrepancyFunctionsArgs){
    dataNames <- discrepancyFunctionsArgs[['dataNames']]
  },
  run = function(){
    disc <- model$getLogProb(dataNames)
    returnType(double(0)) 
    return(disc)
  }
)

calc_chi <- nimbleFunction(
  
  run = function(y = double(2), 
                 y_exp = double(2))
  {
    returnType(double(0))
    
    # expected values
    y_exp[i] <- z[i] * p * nVisits
    y_exp_rep[i] <- z_rep[i] * p * nVisits
    
    # calculate chi squared discrepancy measure
    chi_out <- (y - y_exp) ^ 2 / (y + 1e-6)
    
    return(sum(chi_out))
  }
)

calc_ratio <- nimbleFunction(
  
  run = function(y = double(2), 
                 y_exp = double(2))
  {
    returnType(double(0))
    
    # calculate likelihood ratio discrepancy measure
    ratio_out <- 2 * (y * log((y + 1e-6) / (y_exp + 1e-6)) + (1 - y) * 
                        log((1 - y + 1e-6)/(1 - y_exp + 1e-6)))
    
    return(sum(ratio_out))
  }
)

calc_tukey <- nimbleFunction(
  
  run = function(y = double(2), 
                 y_exp = double(2))
  {
    returnType(double(0))
    
    # calculate freeman tukey discrepancy measure
    tukey_out <- (sqrt(y) - sqrt(y_exp)) ^ 2
    
    return(sum(tukey_out))
  }
)

calc_chi_betabin <- nimbleFunction(
  
  run = function(y = double(1), 
                 y_exp = double(1))
  {
    returnType(double(0))
    
    # calculate chi squared discrepancy measure
    chi_out <- (y - y_exp) ^ 2 / (y + 1e-6)
    
    return(sum(chi_out))
  }
)

calc_ratio_betabin <- nimbleFunction(
  
  run = function(y = double(1), 
                 y_exp = double(1),
                 nVisits = double(0))
  {
    returnType(double(0))
    
    # calculate likelihood ratio discrepancy measure
    ratio_out <- 2 * (y * log((y + 1e-6) / (y_exp + 1e-6)) + (nVisits - y) * 
                        log((nVisits - y + 1e-6)/(nVisits - y_exp + 1e-6)))
    
    return(sum(ratio_out))
  }
)

calc_tukey_betabin <- nimbleFunction(
  
  run = function(y = double(1), 
                 y_exp = double(1))
  {
    returnType(double(0))
    
    # calculate freeman tukey discrepancy measure
    tukey_out <- (sqrt(y) - sqrt(y_exp)) ^ 2
    
    return(sum(tukey_out))
  }
)
