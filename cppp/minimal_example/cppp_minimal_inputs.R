library(nimble)

# load CPPP functions
source("cppp/sally_code/calculateCPPP.R")

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

# inits <- function(y) list(psi = runif(1, 0, 1), p = runif(1, 0, 1), 
#                           z = pmin(y, 1)) 
# 
# # add initial values
# model$setInits(inits(model$y))

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
  setup = function(model, discrepancyFunctionsArgs){
    nVisits <- discrepancyFunctionsArgs[["nVisits"]]
  },
  run = function() {
    
    # get y_exp
    y_exp <- model$z * model$p * nVisits
    
    # calculate chi squared discrepancy measure
    chi_out <- (model$y - y_exp) ^ 2 / (model$y + 1e-6)
    
    returnType(double(0)) 
    return(sum(chi_out))
  }
)

## ratio discrepancy function
ratioDiscFunction <- nimbleFunction(
  contains = discrepancyFunction_BASE,
  setup = function(model, discrepancyFunctionsArgs){
    nVisits <- discrepancyFunctionsArgs[["nVisits"]]
  },
  run = function() {
    
    # get y_exp
    y_exp <- model$z * model$p * nVisits
    
    # calculate likelihood ratio discrepancy measure
    ratio_out <- 2 * (model$y * log((model$y + 1e-6) / 
                                      (y_exp + 1e-6)) + (nVisits - model$y) * 
                        log((nVisits - model$y + 1e-6) / 
                              (nVisits - y_exp + 1e-6)))
    
    returnType(double(0)) 
    return(sum(ratio_out))
  }
)

## tukey discrepancy function
tukeyDiscFunction <- nimbleFunction(
  contains = discrepancyFunction_BASE,
  setup = function(model, discrepancyFunctionsArgs){
    nVisits <- discrepancyFunctionsArgs[["nVisits"]]
  },
  run = function() {
    
    # get y_exp
    y_exp <- model$z * model$p * nVisits
    
    # calculate freeman tukey discrepancy measure
    tukey_out <- (sqrt(model$y) - sqrt(y_exp)) ^ 2
    
    returnType(double(0)) 
    return(sum(tukey_out))
  }
)

## deviance discrepancy function
devianceDiscFunction <- nimbleFunction(
  contains = discrepancyFunction_BASE,
  setup = function(model, discrepancyFunctionsArgs){
    nVisits <- discrepancyFunctionsArgs[["nVisits"]]
  },
  run = function() {
    
    dev_out <- -2 * 
      log(dbinom(model$y, nVisits, model$z * model$p) + 1e-6)
    
    returnType(double(0)) 
    return(sum(dev_out))
  }
)


###################
# function inputs #
###################

dataNames <- "y"

condition <- TRUE

if (condition) {
  # if conditioning on latent state
  paramNames <- c("p", "psi", "z")
  simNodes <- "y"
} else {
  # if *not* conditioning on latent state
  paramNames <- c("p", "psi")
  simNodes <- c("z", "y")
}

discrepancyFunctions <- list(#chisqDiscFunction, 
                             ratioDiscFunction,
                             tukeyDiscFunction, devianceDiscFunction)
discrepancyFunctionsArgs <- list(#list(nVisits = nVisits),
                                 list(nVisits = nVisits),
                                 list(nVisits = nVisits),
                                 list(nVisits = nVisits))


# calculate discrepancies
calcDiscrepanciesFun <- calcDiscrepancies(model,
                                          dataNames,
                                          paramNames,
                                          simNodes,
                                          discrepancyFunctions, 
                                          discrepancyFunctionsArgs)

out_disc <- calcDiscrepanciesFun$run(MCMCOutput[1:5, ])

# calculate PPP
out_ppp <- calculatePPP(MCMCOutput[1:5, ],
                        calcDiscrepanciesFun, 
                        returnDiscrepancies = TRUE)

# run calibration
origMCMCSamples <- MCMCOutput[1:5, ]
mcmcConfFun <- NULL
nCalibrationReplicates <- NULL
MCMCcontrol = list(niter = 500,
                   thin = 1,
                   nburnin = 0)               
returnSamples <- TRUE                            
returnDiscrepancies <- TRUE
calcDisc <- TRUE
parallel <- FALSE 
nCores <- 1

out_cal <- runCalibration(model, dataNames, paramNames, 
                          origMCMCSamples, mcmcConfFun,
                          discrepancyFunctions, discrepancyFunctionsArgs,
                          nCalibrationReplicates, MCMCcontrol,     
                          returnSamples, returnDiscrepancies, calcDisc,
                          parallel, nCores)
