### beta-binomial detection ###

# params requires 2 values 
#   1. lambda - mean abundance across sites
#   2. p - probability of detection
# rho is the correlation parameter on the interval [0, 1); 
# (0 corresponds to the binomial distribution)
simulate_nmixture_betabin <- function(params, nSites, nVisits, rho) {
  
  y <- matrix(NA, nrow = nSites, ncol = nVisits)
  
  N <- rpois(nSites, params$lambda)
  
  for (i in 1:nSites) {
    
    y[i, ] <- VGAM::rbetabinom(nVisits, N[i], params$p, rho)
    
  }
  
  return(y)
}


### negative binomial abundance ###

# params requires 2 values 
#   1. lambda - mean abundance across sites
#   2. p - probability of detection
# r describes the level of overdispersion in the negative binomial
simulate_nmixture_negbin <- function(params, nSites, nVisits, r) {
  
  prob <- r / (1 + r)
  size <- r * params$lambda
  
  y <- matrix(NA, nrow = nSites, ncol = nVisits)
  
  N <- rnbinom(nSites, size = size, prob = prob)
  
  for (i in 1:nSites) {
    
    y[i, ] <- rbinom(nVisits, N[i], params$p)
    
  }
  
  return(y)
}


