library(nimble)

# load CPPP functions
source("cppp/sally_code/calculateCPPP.R")

## Read newcomb data
data <- list(y = read.table("cppp/newcomb_example/light.txt")$V1)

##------------------------------------##
## Newcomb - speed of light data 
##------------------------------------##
modelCode <- nimbleCode({
	for(i in 1:n){
		y[i] ~ dnorm(mu, sd = sigma)
	}

	## noninformative priors
	mu ~ dflat()
	log(sigma) ~ dflat()
})

##------------------------------------##
inits <- list(mu = 0, log_sigma = 2)
constants <- list(n = length(data$y))

## model

model <- nimbleModel(code 		= modelCode, 
					 data 		= data, 
					 inits 		= inits, 
					 constants 	= constants)

cModel 	<- compileNimble(model)
mcmc    <- buildMCMC(model, monitors = c("mu", "log_sigma"))
cMcmc   <- compileNimble(mcmc, project = model)

samples <- runMCMC(cMcmc, niter = 5000, nburnin = 1000)

## save and re-read
##	saveRDS(samples, file = paste0(dirExample, "/MCMCSamples.rds"))
## origMCMCSamples <- readRDS(paste0(dirExample, "/MCMCSamples.rds"))



##------------------------------------##
## base virtual class
##------------------------------------##

## discrepancy base class
discrepancyFunction_BASE <- nimbleFunctionVirtual(
  run = function() returnType(double())
)

##------------------------------------##
## Discrepancies for newcomb example
##------------------------------------##

## 1. minimum observation
minObservation <- nimbleFunction(
  contains = discrepancyFunction_BASE,
  setup = function(model, discrepancyFunctionsArgs){
    dataNames <- discrepancyFunctionsArgs[['dataNames']]
  },
  run = function(){
    dataVals <- values(model, dataNames)
    stat <- min(dataVals)
    returnType(double(0))  # must return a double(0)
    return(stat)
  }
)


test <- minObservation(model = model, 
                       discrepancyFunctionsArgs = list(dataNames = names(data)))

test$run()
ctest <- compileNimble(test, project = model)
ctest$run()
##############

## make a function to call R sort function in nimble
sortR <- nimbleRcall(prototype 	= function(x = double(1)){}, 
					 Rfun 		= "sort", 
					 returnType = double(1))


asymmetryDisc <- nimbleFunction(
  contains = discrepancyFunction_BASE,
  setup = function(model, discrepancyFunctionsArgs){
    dataNames <- discrepancyFunctionsArgs[['dataNames']]
   	meanParam <- discrepancyFunctionsArgs[['meanParam']]
  },
  run = function(){
    dataVals    <- values(model, dataNames)
    meanVal     <- values(model, meanParam)[1]
    sortedVals  <- sortR(dataVals)
    stat 	      <- abs(sortedVals[61] - meanVal) - abs(sortedVals[6] - meanVal)
    returnType(double(0))  # must return a double(0)
    return(stat)
})


test2 <- asymmetryDisc(model = model, 
                       discrepancyFunctionsArgs = 
                       list(dataNames = names(data), 
                            meanParam = "mu"))

test2$run()
ctest2 <- compileNimble(test2, project = model)
ctest2$run()

##################
## Set values to run CPPP
origMCMCSamples <- samples
paramNames    <- colnames(origMCMCSamples)
dataNames     <- names(data)
mcmcConfFun   <- NULL


## discrepancy functions
discrepancyFunctions <- list(minObservation, asymmetryDisc)
discrepancyFunctionsArgs <- list(list(dataNames = names(data)), 
                 list(dataNames = names(data), 
                    meanParam = "mu"))


nCalibrationReplicates <- 1000 
nIterMCMC <- 200

## Settings for MCMC
MCMCcontrol <- list(niter   = nIterMCMC,
                    nburnin = 0, 
                    thin  = 1)

## Flags
returnSamples <- TRUE                            
returnDiscrepancies <- TRUE
calcDisc <- TRUE
parallel <- FALSE 
nCores <- 1



## run calibration
time <- system.time(
  out <- runCalibration(## this comes from model.R file
                      model = model,
                      dataNames = dataNames,
                      paramNames = paramNames,
                      origMCMCSamples = origMCMCSamples,         
                      mcmcConfFun = mcmcConfFun,
                      discrepancyFunctions = discrepancyFunctions, 
                      discrepancyFunctionsArgs = discrepancyFunctionsArgs,
                      ## this comes from args or defaults
                      nCalibrationReplicates = nCalibrationReplicates,
                      MCMCcontrol = MCMCcontrol,                  
                      returnSamples = returnSamples,
                      returnDiscrepancies = returnDiscrepancies,
                      calcDisc = calcDisc, 
                      parallel = parallel, 
                      nCores   = nCores) 
)
out$time <- time
## save results


