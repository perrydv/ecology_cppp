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
      stat <- (obs_y[i] - y_exp) ^ 2 / (y_exp + 1e-6)
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

########################################################
# discrepancy function not conditioned on latent state #
########################################################

## chi-squared discrepancy function
chisqDiscFunction_NoLatent <- nimbleFunction(
  contains = discrepancyFunction_BASE,
  setup = function(model, discrepancyFunctionsArgs){
    
    dataNames <- discrepancyFunctionsArgs[["dataNames"]]
    nVisits <- discrepancyFunctionsArgs[["nVisits"]]
    prob_occupancy <- discrepancyFunctionsArgs[["prob_occupancy"]]
    prob_detection <- discrepancyFunctionsArgs[["prob_detection"]]
  },
  run = function() {
    psi <- values(model, prob_occupancy)[1]
    ## values 
    p <- values(model, prob_detection)[1]
    obs_y <- values(model, dataNames)
    # get y_exp
    y_exp <- psi * p * nVisits
    
    # get number of sites
    nsites <- length(obs_y)
    
    chi_out <- 0
    for(i in 1:nsites){
      stat <- (obs_y[i] - y_exp) ^ 2 / (y_exp + 1e-6)
      chi_out <- chi_out + stat
    }
    
    returnType(double(0)) 
    return(chi_out)
  }
)

## ratio discrepancy function
ratioDiscFunction_NoLatent <- nimbleFunction(
  contains = discrepancyFunction_BASE,
  setup = function(model, discrepancyFunctionsArgs){
    
    dataNames <- discrepancyFunctionsArgs[["dataNames"]]
    nVisits <- discrepancyFunctionsArgs[["nVisits"]]
    prob_occupancy <- discrepancyFunctionsArgs[["prob_occupancy"]]
    prob_detection <- discrepancyFunctionsArgs[["prob_detection"]]
  },
  run = function() {
    psi <- values(model, prob_occupancy)[1]
    ## values 
    p <- values(model, prob_detection)[1]
    obs_y <- values(model, dataNames)
    # get y_exp
    y_exp <- psi * p * nVisits
    
    # get number of sites
    nsites <- length(obs_y)
    
    ratio_out <- 0
    
    for(i in 1:nsites){ 
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
tukeyDiscFunction_NoLatent <- nimbleFunction(
  contains = discrepancyFunction_BASE,
  setup = function(model, discrepancyFunctionsArgs){
    
    dataNames <- discrepancyFunctionsArgs[["dataNames"]]
    nVisits <- discrepancyFunctionsArgs[["nVisits"]]
    prob_occupancy <- discrepancyFunctionsArgs[["prob_occupancy"]]
    prob_detection <- discrepancyFunctionsArgs[["prob_detection"]]
  },
  run = function() {
    psi <- values(model, prob_occupancy)[1]
    ## values 
    p <- values(model, prob_detection)[1]
    obs_y <- values(model, dataNames)
    # get y_exp
    y_exp <- psi * p * nVisits
    
    # get number of sites
    nsites <- length(obs_y)
    
    tukey_out <- 0
    
    for(i in 1:nsites){ 
      tukey_out <- tukey_out + (sqrt(obs_y[i]) - sqrt(y_exp)) ^ 2
    }
    
    returnType(double(0)) 
    return(tukey_out)
    
  }
)

## deviance discrepancy function
devianceDiscFunction_NoLatent <- nimbleFunction(
  contains = discrepancyFunction_BASE,
  setup = function(model, discrepancyFunctionsArgs){
    
    dataNames <- discrepancyFunctionsArgs[["dataNames"]]
    nVisits <- discrepancyFunctionsArgs[["nVisits"]]
    prob_occupancy <- discrepancyFunctionsArgs[["prob_occupancy"]]
    prob_detection <- discrepancyFunctionsArgs[["prob_detection"]]
  },
  run = function() {
    psi <- values(model, prob_occupancy)[1]
    ## values 
    p <- values(model, prob_detection)[1]
    obs_y <- values(model, dataNames)
    
    # get number of sites
    nsites <- length(obs_y)    
    
    dev_out <- 0
    
    for(i in 1:nsites){ 
      dev_out <- dev_out + -2 * 
        log(dbinom(obs_y[i], nVisits, psi * p) + 1e-6)
    }
    
    returnType(double(0)) 
    return(dev_out)
    
  }
)

#########################################
# discrepancy function on occupancy (z) #
#########################################



## chi-squared discrepancy function
chisqDiscFunction_z <- nimbleFunction(
  contains = discrepancyFunction_BASE,
  setup = function(model, discrepancyFunctionsArgs){
    
    dataNames <- discrepancyFunctionsArgs[["dataNames"]]
    nVisits <- discrepancyFunctionsArgs[["nVisits"]]
    prob_occupancy <- discrepancyFunctionsArgs[["prob_occupancy"]]
    latent_occ <- discrepancyFunctionsArgs[["latent_occ"]]
  },
  run = function() {
    psi <- values(model, prob_occupancy)[1]
    ## values 
    obs_z <- values(model, latent_occ)
    # get number of sites
    nsites <- length(obs_z)
    
    chi_out <- 0
    for(i in 1:nsites){
      stat <- (obs_z[i] - psi) ^ 2 / (psi + 1e-6)
      chi_out <- chi_out + stat
    }
    
    returnType(double(0)) 
    return(chi_out)
  }
)

## tukey discrepancy function
tukeyDiscFunction_z <- nimbleFunction(
  contains = discrepancyFunction_BASE,
  setup = function(model, discrepancyFunctionsArgs){
    
    dataNames <- discrepancyFunctionsArgs[["dataNames"]]
    nVisits <- discrepancyFunctionsArgs[["nVisits"]]
    prob_occupancy <- discrepancyFunctionsArgs[["prob_occupancy"]]
    latent_occ <- discrepancyFunctionsArgs[["latent_occ"]]
  },
  run = function() {
    psi <- values(model, prob_occupancy)[1]
    ## values 
    obs_z <- values(model, latent_occ)
    
    # get number of sites
    nsites <- length(obs_z)
    
    tukey_out <- 0
    
    for(i in 1:nsites){ 
      tukey_out <- tukey_out + (sqrt(obs_z[i]) - sqrt(psi)) ^ 2
    }
    
    returnType(double(0)) 
    return(tukey_out)
    
  }
)