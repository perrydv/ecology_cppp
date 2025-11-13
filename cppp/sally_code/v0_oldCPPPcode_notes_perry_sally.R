###########################################################################
## This file contains CPPP code ## Sally
###########################################################################
## Overall organization is as follows:
## runCPPP is the master function for the whole procedure.
##    - get an MCMC sample for the real data if not provided
##    - run the PPP for the MCMC sample for the real data, including standard deviation on the "P-value" using moving block bootstrap (MBB)
##    - create a CPPTask function for parallelization.  Each CPPTask does some iterations of:
##            - simulate fake data from the prior (to be discussed)
##            - run a (short) MCMC for the fake data
##            - run the PPP for the MCMC sample for the fake data, including uncertainty by MBB.
##    - calculate expected value of number of simulated PPPs that are less than real-data PPP (we could also come up with a standard error for this).
##
## Each PPP is run as follows:
## simulatePPP is the master function for calculating a CPPP and its MBB std. dev.
##     - calculatePPP (note similar names) is a nimbleFunction that does the main work.
##
## simulatePPP works as follows:
##     - for each row of the posterior sample,
##           - put the values back in the model,
##           - calculate the discrepancy function for the real data
##           - simulate new data
##           - calculate the discrepancy function for the new data
###########################################################################
## Some high-level notes.
## 1. We want to support doing the calibration from a choice other than the prior, notably the real-data posterior.
## 2. Some variables could use better names.  This seems to be common as software is developed.  An original notion of a name becomes out of date with how a variable is used, or we see a collection of names that need better disambiguation from each other.
## 3. An alternative to MBB might be to use a single estimate of ESS and determine Monte Carlo error from that.
##    MBB has some nice aspects.  We would need to come up with block sizes.  A draw back to MBB is that it limits somewhat the shortness of runs that we can use.  This issue -- how to estimate uncertainties for many short runs -- is perhaps the most technical nugget of the project.
## 4. I think calculatePPP has some (potentially serious) redundancies.
## 5. We will want to give some thought to flexibility required to simulate new data conditioned on latent states or new latent states & new data.
## 6. We should replace syntax like "calculate(model, nodes)" with "model$calculate(nodes)".  Both work, but the latter is our more modern and sensible syntax.

## I put some detailed comments in code.
###########################################################################
## Sally
## LOOK FOR ## REVISE ## FOR CODE TO CERTAINLY ADJUST
###########################################################################
## High level design ideas (Perry and Sally discussion 1/14/21)
##
## calculatePPP: Return a matrix with observed and simulated discrepancies
##               Possibly: option if observed does not need to be calculated?
##               Iterate over all rows of the MCMCoutput (not random rows, not random blocks)
##
## Would it make sense to break apart the steps in "simulatePPP"?
## Somehow I think of one function to run the MCMC calibration replicate, one function to get discrepancies, 
## and one function to calculate PPP from discrepances.
## (Maybe that leads to something like simulatePPP, but with a less confusing (to me!) name?
##
## Design goals for us to keep in mind:
##   - Multiple discrepancy measures
##   - Discrepancy measures applied to saved MCMC samples without re-running everything.
##   - Need to insert steps to estimate uncertainty of final result.
###########################################################################


## virtual functions
pppFuncVirtual <- nimbleFunctionVirtual(
    run = function(N = integer(0)) returnType(double(0))
)

discrepancyFunction_BASE <- nimbleFunctionVirtual(
    run = function() returnType(double(0))
)
###########################################################################
## Nimble function to compute one PPP
###########################################################################

calculatePPP <- nimbleFunction(
    setup = function(model,
                     dataNames,
                     paramNames,
                     MCMCIter,
                     thin,
                     burnIn,
                     discrepancyFunction,
                     discrepancyFunctionArgs){

        ## Get parameter and data dependencies
        paramDependencies <- model$getDependencies(paramNames)
        dataDependencies  <- model$getDependencies(dataNames)
        
        ## REVISE ## How to support multiple discrepancies?? 
        ## SP create a nimbleFunctionList (to allow multiple discrepancy functions)
        discrepancyFunctionList <- nimbleFunctionList(discrepancyFunction_BASE)
        discrepancyFunctionList[[1]] <- discrepancyFunction(model, discrepancyFunctionArgs)

    },
    run = function(MCMCOutput = double(2), nMCMCSamples = integer(0)){
        countDiscrepancy   <- numeric(nMCMCSamples)
        origDataValues     <- values(model, dataNames)
        
        for(i in 1:nMCMCSamples){
            randNum <- ceiling(runif(1, 0, (MCMCIter)/thin - burnIn - 1 ))
            
            ## put MCMC values from sampled iteration in the model
            values(model, paramNames) <<- MCMCOutput[burnIn + randNum, ]
            
            ## calculate deviance from real data
            model$calculate(paramDependencies)

            ## REVISE ##
            ## PdV: This seems redundant.  We should be able to do this once for each element 
            ## in the posterior.  A possible exception is if the discrepancy function itself does any simulation.
            obsDiscrepancy <- discrepancyFunctionList[[1]]$run() 

            ## simulate data from model, then calculate deviance
            model$simulate(dataNames, includeData = TRUE)
            model$calculate(dataDependencies)
            simDiscrepancy <- discrepancyFunctionList[[1]]$run()

            ## count the number of simulated data sets where the
            ## simulated data has a higher deviance than the observed
            if(simDiscrepancy >= obsDiscrepancy) countDiscrepancy[i] <- 1
            values(model, dataNames) <<- origDataValues
            model$calculate(dataDependencies)
        }
        
        ## REVISE ## is this ok as a mean or need to correct for nMCMCSamples + 1
        ppp <- mean(countDiscrepancy)
        returnType(double(0))
        return(ppp)
    }
)


###########################################################################
## Nimble function to compute CPPP
###########################################################################

simulatePPP <- function(MCMCIter,
                         burnIn,
                         nMCMCSamples,
                         cCalcPPP,
                         cMCMC = NULL,
                         samples = NULL,
                         returnSamples){ ## PdV: This seems to need a default or an "if(missing(returnSamples))"
    ## calculate proportion of deviances that are greater or equal to the
    ## observed value
    if(is.null(samples)){
        cMCMC$run(MCMCIter, reset = TRUE) ## PdV: Should we use runMCMC?
        samples <- as.matrix(cMCMC$mvSamples)
    }
    
    modelPPP <- cCalcPPP$run(samples, nMCMCSamples)
    
    if(!is.finite(modelPPP))    modelPPP <- NA
    return(list(modelPPP = modelPPP,
                samples = if(returnSamples == TRUE) samples else NULL))
}

#' Calculates the calibrated posterior predictive p-value (CPPP) for a given NIMBLE model.
#'
#' EXPERIMENTAL.  Takes a NIMBLE model and calculates the CPPP for it.  The CPPP measures the degree of  "surprise" in the observed data given the specified model.
#'
#' @param model an uncompiled nimble model
#' @param dataNames character vector of the names of the data nodes of the model that should be re-simulated to calculate the CPPP.  
#'  It is allowed to leave out some names of data elements in the model, in which case only the data names provided here will be re-simulated for CPPP calculation.
#' @param mcmcCreator (optional) a function that takes the \code{model} as its sole argument, and returns an uncompiled nimble MCMC algorithm created from a call to \code{buildMCMC()}.
#' This argument can be used to provide a custom MCMC algorithm for use in CPPP calculation (e.g. an algorithm with custom samplers for different model nodes).
#' If this argument not provided, NIMBLE's default MCMC samplers will be used.
#' @param discrepancyFunction a \code{nimbleFunction} specifying the discrepancy function to be used.  See details.
#' @param discrepancyFunctionArgs an R \code{list} containing informaton to be used in the \code{discrepancyFunction}.  See details.
#' @param cpppControl	(optional) an R \code{list} with parameters governing the CPPP calculation.  See details for specific parameters.
#' @param mcmcControl (optional) an R \code{list} with parameters governing the MCMC  algorithm,  See details for specific parameters.
#' @param origMCMCOutput (optional) a matrix with samples from the posterior distribution of the model for all model parameters.
#' Each matrix column should contain samples from a different model parameter, and should be named as that parameter.  This argument should be supplied if
#' a user has already run an MCMC and obtained posterior samples from the model, and wishes to use these samples in the CPPP algorithm.
#' @param returnSamples a logical argument indicating whether to return all posterior samples from all simulated data sets.  This can result in a very large
#' returned object,  Defaults to \code{FALSE}.
#' @param nCores an integer argument indicating the number of cpu cores to use for CPPP calculation.  In order to use >1 core, the \code{parallel} R package must be
#' installed.  Only Mac and Linux operating systems support multiple cores at this time.
#'
#' @author Nicholas Michaud and Lauren Ponisio
#' @export
#' @details The CPPP algorithm used is based on that introduced by Hjoort et al. (2006), with some modifications.  Specfically, for each simulated data set and posterior predictive p-value
#' calculated from that data set, a bootstrap estimate of the Monte-Carlo error in the p-value is calculated as well.  These errors are used to produce an error-adjusted calibrated posterior predictive p-value.
#'
#'  @section The \code{discrepancyFunction} Argument:
#' The \code{discrepancyFunction} argument is a \code{nimbleFunction} that will calculate the discrepancy between the model and the observations.  The choice of discrepancy function is highly model-specific,
#' and as such no default is provided for this argument.  The \code{discrepancyFunction} provided must follow these rules:
#' \enumerate{
#'   \item The \code{discrepancyFunction} must be a \code{nimbleFunction} with \code{setup} code and {run} code.  See the User's Manual for more information on writing \code{nimbleFunctions}.
#'   \item The \code{setup} code to \code{discrepancyFunction} must accept two arguments: \code{model},
#'   which is the model the CPPP is being calculated for, and \code{discrepancyFunctionArgs},
#'   which is an R list with any needed information for calculating the discrepancy (e.g. node names, tuning parameters, etc.).
#'   \item The \code{run} code of the discrepancy function must take no arguments and return a scalar double, and must \code{contain} the \code{discrepancyFunction_BASE} virtual \code{nimbleFunction}.
#' }
#' The \code{discrepancyFunctionArgs} argument to \code{runCPPP} is passed to the \code{discrepancyFunction}.  See the examples section below for an example of a discrepancy function.
#'
#'
#' @section The \code{cpppControl} Argument:
#' The \code{cpppControl} argument is a list with the following elements:
#'	\itemize{
#'	  \item{\code{nCalibrationReplicates}}	{
#'      an integer argument determining how many simulated posterior predictive p-values should be calculated to estimate the distribution of p-values.  Defaults to 100.
#'	  }
#'	  \item{\code{nMCMCSamples}}{
#'	    an integer argument determining how many posterior samples should be used to calculate the posterior predictive p-value for each simulated data set.  Defaults to 1000.
#'	  }
#'  }
#' @section The \code{mcmcControl} Argument:
#' The \code{mcmcControl} argument is a list with the following elements:
#'	\itemize{
#'	  \item{\code{nMCMCIters}}	{
#'      an integer argument determining how many MCMC iterations should be run for each posterior predictive p-value calculation.  
#'      Defaults to 10000, but should probably be manually set.
#'	  }
#'	\item{\code{burnInProp}}{
#'	   the proportion of the MCMC chain to discard as burn-in for each posterior predictive p-value calculation.  Must be between 0 and 1.  Defaults to 0.
#'	}
#'  }
#'
#' @return
#' an R list with four elements:
#' \itemize{
#' \item{\code{CPPP}} {The CPPP value, measuring the amount of "surprise" the data has given the model.  Lower values indicate greater "surprise".  Will be between 0 and 1.}
#' \item{\code{observed.PPP}} {The posterior predictive p-value for the observed data, along with an estimate of the Monte Carlo error for this p-value.}
#' \item{\code{simulated.PPPs}} {Posterior predictive p-values for each of the simulated data sets, along with error estimates for each p-value.  The set of simulated posterior predictive p-values provide an empirical distribution of p-values, which is used to calculate the CPPP.}
#' \item{\code{samples}} {An R list, only returned when \code{returnSamples = TRUE}.  Each element of this list will be a matrix of posterior samples from the model given a set of simulated data.  There will be \code{nCalibrationReplicates} sets of samples.}
#' }
#'
#' @references
#' Hjort, N. L., Dahl, F. A., & Steinbakk, G. H. (2006). Post-processing posterior predictive p values. \emph{Journal of the American Statistical Association}, 101(475), 1157-1174.
#' @examples
#' \dontrun{
#'
#'
#' pumpCode <- nimbleCode({
#'   for (i in 1:N){
#'     theta[i] ~ dgamma(alpha, beta)
#'     lambda[i] <- theta[i]*t[i]
#'     x[i] ~ dpois(lambda[i])
#'   }
#'   alpha ~ dexp(1.0)
#'   beta ~ dgamma(0.1,1.0)
#' })
#'
#' pumpConsts <- list(N = 10,
#'                    t = c(94.3, 15.7, 62.9, 126, 5.24,
#'                          31.4, 1.05, 1.05, 2.1, 10.5))
#'
#' pumpData <- list(x = c(5, 1, 5, 14, 3, 19, 1, 1, 4, 22))
#'
#' pumpModel <- nimbleModel(code=pumpCode,
#'                          constants=pumpConsts,
#'                          data=pumpData)
#'
#' # Define the discrepancy function.
#' # Note that the discrepancy function must be a nimbleFunction that takes
#' # model and discFunctionArgs as its only two arguments.
#'
#' pumpDiscFunction <- nimbleFunction(
#'   setup = function(model, discFunctionArgs){
#'     dataNames <- discFunctionArgs[['dataNames']]
#'   },
#'   run = function(){
#'     dataVals <- values(model, dataNames)
#'     lambdaVals <- values(model, 'lambda')
#'     disc <- 0
#'     for(i in 1:dim(dataVals)[1]){
#'       disc <- disc + ((dataVals[i] - lambdaVals[i])/sqrt(lambdaVals[i])) # subtract the expected value, divide by the sd.
#'                                                                          # Note here the expected value and sd are the same,
#'                                                                          # as the data comes from a poisson distribution.
#'     }
#'     returnType(double(0))  # must return a double(0)
#'     return(disc)
#'   },
#'   contains = discrepancyFunction_BASE
#' )
#'
#' pumpCPPP <- runCPPP(pumpModel,
#'                     dataNames = 'x',
#'                     discrepancyFunction = pumpDiscFunction,
#'                     discrepancyFunctionArgs = list(dataNames = 'x'),
#'                     mcmcControl = list(nMCMCIters = 20000))
#' }

runCPPP <-  function(model,
                     dataNames, ## names of the data ## PdV: Maybe this should default to all data nodes.  We can check if(missing(dataNames)) {...}
                     mcmcCreator = NULL,
                     discrepancyFunction = NULL,
                     discrepancyFunctionArgs = list(),
                     cpppControl = list(), ## rename "control"?
                     mcmcControl = list(), ## 
                     origMCMCOutput = NULL,
                     returnSamples = FALSE,
                     nCores = 1){
    
    ## number of calibration replicates (to estimate ppp null distribution)
    nCalibrationReplicates <- cpppControl[['nCalibrationReplicates']] 
    if(is.null(nCalibrationReplicates)){
        print("Defaulting to 100 simulated PPP values.")
        nCalibrationReplicates <- 100
    }
    
    ## number of samples from posterior to use for ppp calculation    
    nMCMCSamples <- cpppControl[['nMCMCSamples']] ## Our "M"
    if(is.null(nMCMCSamples)){
        nMCMCSamples <- 1000
    }
    
    ## REVISE ##
    ## If mcmcCreator is null, use a "standard" MCMC creator
    if(is.null(mcmcCreator)){ ## PdV: Questions: Should this return an MCMC configuration instead of the built MCMC?
        ## Name as "mcmcBuilder" or "mcmcConfFun"??
        ## Do we also want to support an existing MCMC configuration and/or built MCMC?
        ## (These are maybe secondary details, to work out after we get the basic PPP/CPPP steps in place.)
        ## A common workflow might be to use the *same* built and compiled MCMC for the calibration replicates.
        
    	## SP: check for thinning parameter 
    	## that can be also specified in the mcmcCreator function
        ## PdV: Possible behavior is that one over-rides the other.  E.g. If mcmcControl$thin is provided, maybe it should over-ride the thin set by mcmcCreator (or vice-versa?).
        if(is.null(mcmcControl$thin)) thin <- 1 else thin <- mcmcControl$thin

        ## PdV: It is not even clear if we want to do these steps here, as they might only be necessary inside of each thread if using parallel option and the original MCMC sample is provided.
        newMCMCSpec <- configureMCMC(model,
                                        print=FALSE,
                                        thin=thin)
        newMCMC <- buildMCMC(newMCMCSpec)
    }
    else{
        newMCMC <- mcmcCreator(model)
    }

    
    ##-------------------------------------##
    ## REVISE ##
    ## SP: need to get the thin value from an MCMC object;
    ## not sure what's the best way to do it 
    ## this is a temp solution
    ## PdV: Not sure, but maybe there is a processing flow where newMCMC comes from mcmcCreator, and that sets the thin, and we need to determine that.
    ## I think you pulled out the right value, but we probably want to add a simple "get" method to the built MCMC, like get_thin() or something like that.
    ## Curious when we *really* need to know this.  Maybe if thinning of "real" MCMC vs. calibration replicates is different.
    thin <- newMCMC$thinFromConfVec[1]
  	##-------------------------------------##
  
    ##-------------------------------------##
    ## REVISE ##
    ## following controls to check
    ## NOTICE: the code is set up such that one can use a different number
    ## of MCMC samples to compute the PPP than the number of iterations
    ## for the MCMC for each calibration replicate

    nMCMCIters <- mcmcControl[['nMCMCIters']]  ## Or is this our M?  We need to clarify nMCMCSamples vs nMCMCIters.  Maybe one was for MBB and is now unnecessary?

    ## PdV Sanity: Following is for the *original MCMC output* if provided.
    if(!is.null(origMCMCOutput)){
        origMCMCIter <- nrow(origMCMCOutput)*thin ## PdV: The use of `thin` here looks hazardous.  If we need to know the thin used for the original MCMC sample, we might need that from the user.
        ## If the user provided the MCMC used go generate "original" samples, then we can read the thin value from that.
        if(origMCMCIter <= 1){
            stop("origMCMCOutput has no rows!")
        }
        if(is.null(nMCMCIters)){
            nMCMCIters <- origMCMCIter
        }
        if(nMCMCIters < 2){
            stop("nMCMCIters must be an integer > 1")
        }
    }
    else if(is.null(nMCMCIters)){ ## origMCMCoutput was not provided
        nMCMCIters <- 10000       ## PdV: 10000 seems large, given our whole motivation! 1000? 500?
        ## PdV: Another idea: If the goal is say ESS = 500, and we have an original sample, then we could use it to estimate number of iterations needed to achieve ESS = 500 (or whatever).
        warning("Defaulting to 10,000 MCMC iterations for each MCMC run.")
    }
    else{ ## origMCMCoutput was not provided but nMCMCIters was provided.
        origMCMCIter <- nMCMCIters ## PdV & SP: Not sure it makes sense to set default "original MCMC" size be as small as the replicate sizes.  Also, we might want to run the original MCMC to determine necessary length of replicates!
    }

    ## Set up burnin
    ## PdV: Let's give this fresh thought!
    burnInProp <-  mcmcControl[['burnInProp']]
    if(is.null(burnInProp)){
        burnInProp <- 0
    }
    else if(burnInProp >= 1 | burnInProp < 0){
        stop("burnInProp needs to be between 0 and 1")
    }
    ## PdV & SP: We may not need long burn-in because we should have good inits
    ## PdV: Another issue is about adaptive tuning parameters.  We may want to turn off adaptation in the "short" runs and simply copy "converged" tuning parameters from the long run.  We have some code for that but not nicely built in the package.  Daniel dealt with this previously for Dao's paper, so let's ask Daniel.

    ## Some basic checks:
    if(is.null(discrepancyFunction)){
        stop("No discrepancy function specified.")
    }
    if(!inherits(model, "RmodelBaseClass")){
        stop("model is not an Rmodel")
    }
    if(is(model$CobjectInterface, "uninitializedField")){ ## The user provided n uncompiled a model that has not been compiled.   (We might need to check first if the user provided a compiled model!)
        compileNimble(model)
    }


    ## PdV: This looks like it is trying to check is dataNames is valid.
    ## I am not sure the first check will work correctly, because dataNames may be a character vector.
    testDataNames <- try(model[[dataNames]], silent=TRUE)
    if(inherits(testDataNames, "try-error")){
        stop(paste("dataNames", dataNames,
                   "is not the name of the data in model"))
    } else{ 
        ## PdV: This test makes sense to check if all dataNames are really data nodes in the model.
        ## However, shouldn't it be valid to simulate non-data nodes (random effects) as part of CPPP?
        ## Maybe we should just check that all nodes are in the model and are stochastic?
        test2DataNames <- all(model$expandNodeNames(dataNames) %in%
                              model$getNodeNames(dataOnly=TRUE))
        if(test2DataNames == FALSE){
            stop(paste("dataNames", dataNames,
                       "is not the name of the data in model"))
        }
    }

    ## PdV: We could actually calculate the size of all objects necessary for returnSamples == TRUE.
    if(returnSamples == TRUE){
        warning("returnSamples = TRUE, so all posterior samples from all simulations will be returned.  This may produce
              a very large return object.")
    }

    ## keep track of the real data
    ## PdV: Let's keep an eye out for inefficiency if dataNames is large. values() is inefficient in uncompiled execution.
    origData <- values(model, dataNames)

    ## PdV: Make original MCMC if needed.
    ## Should we consider the case that the use provided a compiled MCMC?
    if(is.null(origMCMCOutput)){
        CnewMCMC <- compileNimble(newMCMC, project = model, resetFunctions = TRUE)
        CnewMCMC$run(nMCMCIters)
        origMCMCOutput <- as.matrix(CnewMCMC$mvSamples)
    }
    burnIn <- ceiling(burnInProp*(nMCMCIters/thin))

    ## PdV: I am not clear what the purpose of this error-trapping is.  Maybe it catches in the newMCMC is for the wrong model?
    ## We could potentially do it more efficiently or skip it?
    paramNames <- colnames(origMCMCOutput)
    testParamNames <- lapply(paramNames, function(x){
        test.this.param <- try(model[[x]], silent=TRUE)
        if(inherits(test.this.param, "try-error")){
            # stop(paste("first stop paramNames", x,
            #            "are not parameters in model"))
        }
    })
    ## PdV: This seems like a more effective error-trapping step.  I'm not sure if the above one would catch anything different?
    test2ParamNames <- all(model$expandNodeNames(paramNames) %in%
                           model$getNodeNames(includeData=FALSE,
                                              stochOnly=TRUE))
    if(test2ParamNames == FALSE){
        browser()
        ## stop(paste("second stop paramNames", paramNames,
        ##          "are not parameters in model"))
    }


    ## Compute original PPP from the model 
    ## Create a calculatePPP function associated to the model
    ## than will be compuled and used simulatePPP
    paramDependencies <- model$getDependencies(paramNames)
    ## This will calculate the observed PPP.
    modelCalculatePPP <- calculatePPP(model                     = model,
                                      dataNames                 = dataNames,
                                      paramNames                = paramNames,
                                      MCMCIter                  = origMCMCIter,
                                      thin                      = thin,
                                      burnIn                    = burnIn,
                                      discrepancyFunction       = discrepancyFunction,
                                      discrepancyFunctionArgs   = discrepancyFunctionArgs)
    cModelCalculatePPP <- compileNimble(modelCalculatePPP, project = model)

    ## calculate observed PPP
    observedPPP <- simulatePPP(MCMCIter         = origMCMCIter,
                               burnIn           = burnIn,
                               nMCMCSamples     = nMCMCSamples,
                               cCalcPPP         = cModelCalculatePPP,
                               samples          = origMCMCOutput,
                               returnSamples    = returnSamples)
################################################################################

## SP: this can be outside the runCPPP function.
    CPPPTask <- function(coreNum, nCalibrationsPerCore) {
        ## REVISE ##
       	## From Perry: We might need replicate=TRUE
        ## replicate = TRUE replicate all current values and data flags from 
        ## given model in the new model.  If TRUE, then the data and 
        ## inits arguments are not used.  
        repModel <- model$newModel(replicate = TRUE) 
        ## REVISE ## EFFICENCY
        ## From Perry: If nCores is 1, it would be nice to not have to rebuild and recompile even this one time,
        ## since we have the objects used outside CPPTask.  This is about efficiency, not correctness.
        cRepModel <- compileNimble(repModel)
        if(is.null(mcmcCreator)){
            repMCMCSpec <-   configureMCMC(repModel,
                                            print=FALSE,
                                            nthin=thin)
            repMCMC <- buildMCMC(repMCMCSpec)
        }
        else{
            repMCMC <- mcmcCreator(repModel)
        }

        repModelCalculatePPP <- calculatePPP(model                     = repModel,
                                             dataNames                 = dataNames,
                                             paramNames                = paramNames,
                                             MCMCIter                  = nMCMCIters,
                                             thin                      = thin,
                                             burnIn                    = burnIn,
                                             discrepancyFunction       = discrepancyFunction,
                                             discrepancyFunctionArgs   = discrepancyFunctionArgs)

        cRepMcmc <- compileNimble(repMCMC, project = repModel)
        cRepModelCalculatePPP <- compileNimble(repModelCalculatePPP, project = repModel)
	
		## Get dependencies of nodes to simulate 
        ## nodes that are monitored during the origin MCMC will not be overwritten
        ## CHECK ## Random effects?
		simNodes <- repModel$getDependencies(paramNames, self=FALSE) 

		out <- list()
        for(iteration in 1:nCalibrationsPerCore){
            message(paste("refitting data iteration", iteration, "for core number", coreNum))
            
            ################################################################################
			## DOUBLE CHECK ## SIMULATION OF NEW DATA
			## To sample form the posterior predictive we use samples values form the MCMC
			## Here I sample indices to use in cRepModel before simulate 
			indexDataGen <- sample(1:(origMCMCIter/thin - burnIn), 1)
 			## put posterior samples for the parameters in the model
			values(cRepModel,paramNames ) <- origMCMCOutput[indexDataGen[], ]
			## simulate data for the calibration replicate
			cRepModel$simulate(simNodes, includeData = TRUE) 
			
			################################################################################
			##--------------------------------------------##
            out[[iteration]] <- simulatePPP(MCMCIter        = nMCMCIters,
                                            burnIn          = burnIn,
                                            nMCMCSamples    = nMCMCSamples,
                                            cCalcPPP        = cRepModelCalculatePPP,
                                            cMCMC           = cRepMcmc,
                                            returnSamples   = returnSamples)
        }
        return(out)
    }
    ################################################################################

    libError <-  try(library('parallel'), silent = TRUE)
    if(inherits(libError, 'try-error') && nCores > 1){
        warning("The 'parallel' package must be installed to use multiple cores for CPPP calculation.  Defaulting to one core.")
        nCores <- 1
    }

    if(nCores > 1){
        ## simulate PPP values
        calibrationOutput <- mclapply(X = 1:nCores, FUN = CPPPTask,
                                      nCalibrationsPerCore =
                                       ceiling(nCalibrationReplicates/nCores))
    }
    else{
        calibrationOutput <- lapply(X = 1, FUN = CPPPTask,
                                    nCalibrationsPerCore =
                                     ceiling(nCalibrationReplicates/nCores))
    }

    ## extract PPP replicates
    replicatedPPP <- sapply(sapply(calibrationOutput,
                             function(x) x), function(x)  x$modelPPP)
    
    ## calculate CPPP
    CPPP <- mean(replicatedPPP <= observedPPP$modelPPP)


    ## extract chains and diagnostics
    if(returnSamples){
        calibrationSamples <- lapply(sapply(calibrationOutput, function(x) x),
                              function(x) x$samples)
    }
    else calibrationSamples <- NULL

    ## output real data and model
    values(model, dataNames) <- origData

    out <- list(CPPP            = CPPP,
                observedPPP     = c(estimate=observedPPP$modelPPP),
                replicatedPPPs  = c(estimate=replicatedPPP),
                samples         = calibrationSamples)
    return(out)
}
