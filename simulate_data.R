
# model 1
simulate_basic <- function(params, nSites, nVisits) {
  
  y <- matrix(NA, nrow = nSites, ncol = nVisits)
  
  z <- rbinom(nSites, 1, params$psi)
  
  for (i in 1:nSites) {
    
    y[i, ] <- rbinom(nVisits, 1, z[i] * params$p)
      
  }
  
  return(y)
}

# model 2
simulate_cov_occ <- function(params, x_site, nSites, nVisits) {
  
  y <- matrix(NA, nrow = nSites, ncol = nVisits)
  
  psi <- rep(NA, nSites)
  
  z <- rep(NA, nSites)
  
  
  for (i in 1:nSites) {
    
    psi[i] <- plogis(params$beta %*% x_site[i, ])
    
    z[i] <- rbinom(1, 1, psi[i])
    
    y[i, ] <- rbinom(nVisits, 1, z[i] * params$p)
    
  }
  
  return(y)
  
} 

# model 3
simulate_cov_occ_det <- function(params, x_site, x_visit, nSites, nVisits) {
  
  y <- matrix(NA, nrow = nSites, ncol = nVisits)
  
  p <- matrix(NA, nrow = nSites, ncol = nVisits)
  
  psi <- rep(NA, nSites)
  
  z <- rep(NA, nSites)
  
  
  for (i in 1:nSites) {
    
    psi[i] <- plogis(params$beta_o %*% x_site[i, ])
    
    z[i] <- rbinom(1, 1, psi[i])
    
    for (j in 1:nVisits) {
      
      p[i, j] <- plogis(params$beta_d %*% x_visit[i, j, ])
      
      y[i, j] <- rbinom(1, 1, z[i] * p[i, j])
      
    }
    
  }
  
  return(y)
  
}

# model 4
simulate_spatial_ranef <- function(params, nRegions, region, nSites, nVisits) {
  
  y <- matrix(NA, nrow = nSites, ncol = nVisits)
  
  psi <- plogis(rnorm(nRegions, params$psi_mean, params$psi_sd))
  
  z <- rep(NA, nSites)
  
  for (i in 1:nSites) {
    
    z[i] <- rbinom(1, 1, psi[region[i]])
    
    y[i, ] <- rbinom(nVisits, 1, z[i] * params$p)
    
  }
  
  return(y)
}
