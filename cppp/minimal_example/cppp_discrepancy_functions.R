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