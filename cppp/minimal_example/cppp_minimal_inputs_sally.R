library(nimble)

# load CPPP functions
source("cppp/sally_code/calculateCPPP_original.R")

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
rho <- 0.01
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

inits <- function(y) list(psi = runif(1, 0, 1), 
                          p = runif(1, 0, 1), 
                          z = pmin(y, 1)) 

# add initial values
model$setInits(inits(model$y))
model_uncompiled$setInits(inits(model_uncompiled$y))

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

#########################
## chi-squared discrepancy function
#########################
## SP: I rewrote the function because I had some errors
## once compiled. For example 
## model$p * model$z 
## makes the compiler throwing an error because 
## z is a vector and p a scalar

chisqDiscFunction <- nimbleFunction(
  contains = discrepancyFunction_BASE,
  setup = function(model, discrepancyFunctionsArgs){
    dataNames <- discrepancyFunctionsArgs[['dataNames']]

    nVisits <- discrepancyFunctionsArgs[["nVisits"]]
    latent_occ <- discrepancyFunctionsArgs[["latent_occ"]]
    prob_detection <- discrepancyFunctionsArgs[["prob_detection"]]

  },
  run = function() {
    z <- values(model, latent_occ)
    ## values 
    p <- values(model, prob_detection)[1]
    obs_y <- values(model, dataNames)
    
    # get number of individuals
    nind <- length(obs_y)

    chi_out <- 0
    for(i in 1:nind){
      # get y_exp
      y_exp <- z[i] * p * nVisits
      stat <- (obs_y[i] - y_exp)^2/(obs_y[i] + 1e-6)
      chi_out <- chi_out + stat
    }

    returnType(double(0)) 
    return(chi_out)
  }
)

##################################
## SP: Test discrepancy function
##################################
chisqList = list(nVisits = nVisits, 
                  dataNames = "y",
                  latent_occ = "z", 
                  prob_detection = "p")

# test <- chisqDiscFunction(model = model_uncompiled,
#                           discrepancyFunctionsArgs = chisqList)


# ## use uncompiled model for testing
# ## make sure is initialized - use the same inits

# init_vals <- inits(model$y)
# model_uncompiled$setInits(init_vals)
# test$run()

# ## Since the discrepancy is a nimble function 
# ## this needs to be compiled too. Once compiled, 
# ## it will use the compiled model. The runCalibration() function 
# ## compiles the discrepancies internally if they are not. 

# ctest <- compileNimble(test, project = model)
# 
# model$setInits(init_vals)
# ctest$run()

###################################################
## ratio discrepancy function
###################################################
ratioDiscFunction <- nimbleFunction(
  contains = discrepancyFunction_BASE,
  setup = function(model, discrepancyFunctionsArgs){
    
    dataNames <- discrepancyFunctionsArgs[['dataNames']]

    nVisits <- discrepancyFunctionsArgs[["nVisits"]]
    latent_occ <- discrepancyFunctionsArgs[["latent_occ"]]
    prob_detection <- discrepancyFunctionsArgs[["prob_detection"]]

  },
  run = function() {
    z <- values(model, latent_occ)
    p <- values(model, prob_detection)[1]
    obs_y <- values(model, dataNames)
    
    # get number of individuals
    nind <- length(obs_y)
   
    ratio_out <- 0
    
    for(i in 1:nind){ 
      # get y_exp
      y_exp <- z[i] * p * nVisits
      ratio_out <- ratio_out + 
                2 * (obs_y[i] * log((obs_y[i] + 1e-6) / 
                                      (y_exp + 1e-6)) + (nVisits - obs_y[i]) * 
                        log((nVisits - obs_y[i] + 1e-6) / 
                              (nVisits - y_exp + 1e-6)))
    }

    returnType(double(0)) 
    return(ratio_out)

  }
)

ratioList = list(nVisits = nVisits, 
                  dataNames = "y",
                  latent_occ = "z", 
                  prob_detection = "p")

################################# 
### Testing discrepancy
################################# 

ratioTest <- ratioDiscFunction(model = model,
                               discrepancyFunctionsArgs = ratioList)
# init_vals <- inits(model$y)
# model_uncompiled$setInits(init_vals)
# 
ratioTest$run()


# cratioTest <- compileNimble(ratioTest, project = model)
# model$setInits(init_vals)

# cratioTest$run()
####################


###################
# function inputs #
###################

condition_on_latent_states <- TRUE

if (condition_on_latent_states) {
  # if conditioning on latent state
  ## SP: you can provide z in paramNames
  ## this assumes that you use values from the original MCMC
  paramNames <- c("p", "psi", "z")
  dataNames <- "y"
} else {
  # if *not* conditioning on latent state
  ## SP: just omit z from paramNames - you can check which notes
  ## the runCalibration() function simulates just 
  ## checking which are the simNodes 

  ## SP - a note: the runCalibration function assumed that 
  ## the matrix of MCMC samples contains only values of 
  ## sampled for parameters in paramNames
  paramNames <- c("p", "psi")
  dataNames <- "y"
}

## Check which are simulated notes
simNodes <- unique(c(model$expandNodeNames(dataNames), 
  model$getDependencies(paramNames, includeData = FALSE, self=FALSE)))


discrepancyFunctions <- list(#chisqDiscFunction, 
                             ratioDiscFunction,
                             tukeyDiscFunction, 
                             devianceDiscFunction)
discrepancyFunctionsArgs <- list(#chisqList, 
                                 ratioList,chisqList, 
                                 ratioList)


# run calibration
if (condition_on_latent_states) {
  origMCMCSamples <- MCMCOutput
  } else {

  ## paramNames are matched the column of the MCMc samples matrix
  ## SP: I added an internal line of code that makes sure of this
  origMCMCSamples <- MCMCOutput[, c("p", "psi") ]
}

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

## Pass the name of the uncompiled model 

out_cal <- runCalibration(model = model_uncompiled, 
                          dataNames, 
                          paramNames, 
                          origMCMCSamples, 
                          mcmcConfFun,
                          discrepancyFunctions, 
                          discrepancyFunctionsArgs,
                          nCalibrationReplicates,
                          MCMCcontrol,     
                          returnSamples, 
                          returnDiscrepancies, 
                          calcDisc,
                          parallel, 
                          nCores)

ppp <- out_cal$obsPPP
cppp <- rep(NA, length(out_cal$obsPPP))
for (i in 1:length(cppp)) {
  cppp[i] <- mean(out_cal$repPPP[i, ] <= ppp[i])
}
