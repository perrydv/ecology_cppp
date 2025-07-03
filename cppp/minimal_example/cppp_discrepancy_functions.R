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
    
    dataNames <- discrepancyFunctionsArgs[["dataNames"]]
    nVisits <- discrepancyFunctionsArgs[["nVisits"]]
    latent_occ <- discrepancyFunctionsArgs[["latent_occ"]]
    prob_detection <- discrepancyFunctionsArgs[["prob_detection"]]
    
  },
  run = function() {
    z <- values(model, latent_occ)
    ## values 
    p <- values(model, prob_detection)[1]
    obs_y <- values(model, dataNames)
    
    # get number of sites
    nsites <- length(obs_y)
    
    chi_out <- 0
    for(i in 1:nsites){
      # get y_exp
      y_exp <- z[i] * p * nVisits
      stat <- (obs_y[i] - y_exp) ^ 2 / (obs_y[i] + 1e-6)
      chi_out <- chi_out + stat
    }
    
    returnType(double(0)) 
    return(chi_out)
  }
)

## ratio discrepancy function
ratioDiscFunction <- nimbleFunction(
  contains = discrepancyFunction_BASE,
  setup = function(model, discrepancyFunctionsArgs){
    
    dataNames <- discrepancyFunctionsArgs[["dataNames"]]
    nVisits <- discrepancyFunctionsArgs[["nVisits"]]
    latent_occ <- discrepancyFunctionsArgs[["latent_occ"]]
    prob_detection <- discrepancyFunctionsArgs[["prob_detection"]]
    
  },
  run = function() {
    z <- values(model, latent_occ)
    p <- values(model, prob_detection)[1]
    obs_y <- values(model, dataNames)
    
    # get number of sites
    nsites <- length(obs_y)
    
    ratio_out <- 0
    
    for(i in 1:nsites){ 
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

## tukey discrepancy function
tukeyDiscFunction <- nimbleFunction(
  contains = discrepancyFunction_BASE,
  setup = function(model, discrepancyFunctionsArgs){
    
    dataNames <- discrepancyFunctionsArgs[["dataNames"]]
    nVisits <- discrepancyFunctionsArgs[["nVisits"]]
    latent_occ <- discrepancyFunctionsArgs[["latent_occ"]]
    prob_detection <- discrepancyFunctionsArgs[["prob_detection"]]
    
  },
  run = function() {
    z <- values(model, latent_occ)
    p <- values(model, prob_detection)[1]
    obs_y <- values(model, dataNames)
    
    # get number of sites
    nsites <- length(obs_y)
    
    tukey_out <- 0
    
    for(i in 1:nsites){ 
      # get y_exp
      y_exp <- z[i] * p * nVisits
      tukey_out <- tukey_out + (sqrt(obs_y[i]) - sqrt(y_exp)) ^ 2
    }
    
    returnType(double(0)) 
    return(tukey_out)
    
  }
)

## deviance discrepancy function
devianceDiscFunction <- nimbleFunction(
  contains = discrepancyFunction_BASE,
  setup = function(model, discrepancyFunctionsArgs){
    
    dataNames <- discrepancyFunctionsArgs[["dataNames"]]
    nVisits <- discrepancyFunctionsArgs[["nVisits"]]
    latent_occ <- discrepancyFunctionsArgs[["latent_occ"]]
    prob_detection <- discrepancyFunctionsArgs[["prob_detection"]]
    
  },
  run = function() {
    z <- values(model, latent_occ)
    p <- values(model, prob_detection)[1]
    obs_y <- values(model, dataNames)
    
    # get number of sites
    nsites <- length(obs_y)
    
    dev_out <- 0
    
    for(i in 1:nsites){ 
      dev_out <- dev_out + -2 * 
        log(dbinom(obs_y[i], nVisits, z[i] * p) + 1e-6)
    }
    
    returnType(double(0)) 
    return(dev_out)
    
  }
)


##################################
# discrepancy function arguments #
##################################

# chisqList <- list(nVisits = nVisits, 
#                   dataNames = "y",
#                   latent_occ = "z", 
#                   prob_detection = "p")
# 
# ratioList <- list(nVisits = nVisits, 
#                   dataNames = "y",
#                   latent_occ = "z", 
#                   prob_detection = "p")
# 
# tukeyList <- list(nVisits = nVisits, 
#                   dataNames = "y",
#                   latent_occ = "z", 
#                   prob_detection = "p")
# 
# devList <- list(nVisits = nVisits, 
#                 dataNames = "y",
#                 latent_occ = "z", 
#                 prob_detection = "p")