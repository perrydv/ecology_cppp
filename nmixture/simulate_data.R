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


### unmodeled variation in N ###

# rho describes the proportion of variation in population sizes occurring 
# among sites
simulate_nmixture_Nvar <- function(params, nSites, nVisits, rho) {
  
  theta1 <- rho * params$lambda
  theta2 <- params$lambda - theta1
  
  v <- rpois(nSites, lambda = theta1)
  
  A <- matrix(NA, nrow = nSites, ncol = nVisits)
  N <- matrix(NA, nrow = nSites, ncol = nVisits)
  y <- matrix(NA, nrow = nSites, ncol = nVisits)
  
  for (i in 1:nSites) {
    
    A[i, ] <- rpois(nVisits, theta2)
    
    N[i, ] <- A[i, ] + v[i]
    
    for (j in 1:nVisits) {
      
      y[i, j] <- rbinom(1, N[i, j], params$p)
      
    }
  }
  
  return(y)
}


### double counting ###

# gamma describes the probability of double counting
simulate_nmixture_dcount <- function(params, nSites, nVisits, gamma) {
  
  N <- rpois(nSites, params$lambda)
  true_count <- matrix(NA, nrow = nSites, ncol = nVisits)
  y <- matrix(NA, nrow = nSites, ncol = nVisits)
  
  for (i in 1:nSites) {
    
    true_count[i, ] <- rbinom(nVisits, N[i], params$p)
    
    for (j in 1:nVisits) {
      
      y[i, j] <- true_count[i, j] + rbinom(1, true_count[i, j], gamma)
      
    }
  }
  
  return(y)
}


