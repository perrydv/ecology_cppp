#########################
# discrepancy functions #
#########################

## discrepancy base class
discrepancyFunction_BASE <- nimbleFunctionVirtual(
  run = function() returnType(double())
)

## chi-squared discrepancy function
chisqDiscFunction_nmix <- nimbleFunction(
  contains = discrepancyFunction_BASE,
  setup = function(model, discrepancyFunctionsArgs) {
    
    dataNames <- discrepancyFunctionsArgs[["dataNames"]]
    nVisits <- discrepancyFunctionsArgs[["nVisits"]]
    latent_N <- discrepancyFunctionsArgs[["latent_N"]]
    prob_detection <- discrepancyFunctionsArgs[["prob_detection"]]
    
  },
  run = function() {
    
    ## values
    N <- values(model, latent_N)
    p <- values(model, prob_detection)[1]
    obs_y <- matrix(values(model, dataNames), ncol = nVisits)
    
    # get number of sites
    nsites <- length(N)
    
    chi_out <- 0
    for (i in 1:nsites) {
      # get y_exp
      y_exp <- N[i] * p
      for (j in 1:nVisits) {
        stat <- (obs_y[i, j] - y_exp) ^ 2 / (y_exp + 1e-6)
        chi_out <- chi_out + stat
      }
    }
    
    returnType(double(0)) 
    return(chi_out)
  }
)

## ratio discrepancy function
ratioDiscFunction_nmix <- nimbleFunction(
  contains = discrepancyFunction_BASE,
  setup = function(model, discrepancyFunctionsArgs) {
    
    dataNames <- discrepancyFunctionsArgs[["dataNames"]]
    nVisits <- discrepancyFunctionsArgs[["nVisits"]]
    latent_N <- discrepancyFunctionsArgs[["latent_N"]]
    prob_detection <- discrepancyFunctionsArgs[["prob_detection"]]
    
  },
  run = function() {
    
    ## values
    N <- values(model, latent_N)
    p <- values(model, prob_detection)[1]
    obs_y <- matrix(values(model, dataNames), ncol = nVisits)
    
    # get number of sites
    nsites <- length(N)
    
    ratio_out <- 0
    
    for (i in 1:nsites) { 
      # get y_exp
      y_exp <- N[i] * p
      for (j in 1:nVisits) {
        ratio_out <- ratio_out + 
          2 * (obs_y[i, j] * log((obs_y[i, j] + 1e-6) / (y_exp + 1e-6)))
      }
    }
    
    returnType(double(0)) 
    return(ratio_out)
    
  }
)

## tukey discrepancy function
tukeyDiscFunction_nmix <- nimbleFunction(
  contains = discrepancyFunction_BASE,
  setup = function(model, discrepancyFunctionsArgs) {
    
    dataNames <- discrepancyFunctionsArgs[["dataNames"]]
    nVisits <- discrepancyFunctionsArgs[["nVisits"]]
    latent_N <- discrepancyFunctionsArgs[["latent_N"]]
    prob_detection <- discrepancyFunctionsArgs[["prob_detection"]]
    
  },
  run = function() {
    
    ## values
    N <- values(model, latent_N)
    p <- values(model, prob_detection)[1]
    obs_y <- matrix(values(model, dataNames), ncol = nVisits)
    
    # get number of sites
    nsites <- length(N)
    
    tukey_out <- 0
    
    for (i in 1:nsites) { 
      # get y_exp
      y_exp <- N[i] * p
      for (j in 1:nVisits) {
        tukey_out <- tukey_out + (sqrt(obs_y[i, j]) - sqrt(y_exp)) ^ 2
      }
    }
    
    returnType(double(0)) 
    return(tukey_out)
    
  }
)

## deviance discrepancy function
devianceDiscFunction_nmix <- nimbleFunction(
  contains = discrepancyFunction_BASE,
  setup = function(model, discrepancyFunctionsArgs) {
    
    dataNames <- discrepancyFunctionsArgs[["dataNames"]]
    nVisits <- discrepancyFunctionsArgs[["nVisits"]]
    latent_N <- discrepancyFunctionsArgs[["latent_N"]]
    prob_detection <- discrepancyFunctionsArgs[["prob_detection"]]
    
  },
  run = function() {
    
    ## values
    N <- values(model, latent_N)
    p <- values(model, prob_detection)[1]
    obs_y <- matrix(values(model, dataNames), ncol = nVisits)
    
    # get number of sites
    nsites <- length(N)
    
    dev_out <- 0
    
    for (i in 1:nsites) { 
      for (j in 1:nVisits) {
        dev_out <- dev_out + -2 * 
          log(dbinom(obs_y[i, j], N[i], p) + 1e-6)
      }
    }
    
    returnType(double(0)) 
    return(dev_out)
    
  }
)


