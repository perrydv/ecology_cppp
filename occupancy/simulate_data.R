### Model 1 - basic model ###

# params requires 2 values 
#   1. psi - Occupancy probability
#   2. p - detection probability

simulate_basic <- function(params, nSites, nVisits) {
  
  y <- matrix(NA, nrow = nSites, ncol = nVisits)
  
  z <- rbinom(nSites, 1, params$psi)
  
  for (i in 1:nSites) {
    
    y[i, ] <- rbinom(nVisits, 1, z[i] * params$p)
      
  }
  
  return(y)
}

### Model 2a - occupancy covariates ###

# params requires 2 values 
#   1. beta - covariate coefficients on occupancy
#   2. p - detection probability
# x_site is the site-level covariate data
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

### Model 2b - occupancy covariates, interaction ###

# params requires 2 values 
#   1. beta - covariate coefficients on occupancy
#   2. p - detection probability
# x_site is the site-level covariate data
simulate_cov_occ_inter <- function(params, x_site, nSites, nVisits) {
  
  y <- matrix(NA, nrow = nSites, ncol = nVisits)
  
  psi <- rep(NA, nSites)
  
  z <- rep(NA, nSites)
  
  
  for (i in 1:nSites) {
    
    psi[i] <- plogis(params$beta[1] + params$beta[2] * x_site[i, 1] +
                       params$beta[3] * x_site[i, 1] * x_site[i, 2])
    
    z[i] <- rbinom(1, 1, psi[i])
    
    y[i, ] <- rbinom(nVisits, 1, z[i] * params$p)
    
  }
  
  return(y)
  
} 

### Model 3 - occupancy & detection covariates ###

# params requires 2 values 
#   1. beta_o - covariate coefficients on occupancy
#   2. beta_p - covariate coefficients on detection
# x_site is the site-level covariate data
# x_visit is the detection-level covariate data
simulate_cov_occ_det <- function(params, x_site, x_visit, nSites, nVisits) {
  
  y <- matrix(NA, nrow = nSites, ncol = nVisits)
  
  p <- matrix(NA, nrow = nSites, ncol = nVisits)
  
  psi <- rep(NA, nSites)
  
  z <- rep(NA, nSites)
  
  
  for (i in 1:nSites) {
    
    psi[i] <- plogis(params$beta_o %*% x_site[i, ])
    
    z[i] <- rbinom(1, 1, psi[i])
    
    for (j in 1:nVisits) {
      
      p[i, j] <- plogis(params$beta_p %*% x_visit[i, j, ])
      
      y[i, j] <- rbinom(1, 1, z[i] * p[i, j])
      
    }
    
  }
  
  return(y)
  
}

### Model 4 - non-independent sites ###

# params requires 3 values 
#   1. psi_mean - mean psi across regions
#   2. psi_sd - sd psi across regions
#   3. p - probability of detection
# region is a vector of length nSites, indicating the region for each site
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


### Model 5 - beta-binomial ###

# params requires 2 values 
#   1. psi - occupancy probability
#   2. p - sd psi across regions
# rho is the correlation paramter on the interval [0, 1); 
# (0 corresponds to the binnomial distribution)
simulate_betabinomial <- function(params, nSites, nVisits, rho) {
  
  y <- rep(NA, nSites)
  
  z <- rbinom(nSites, 1, params$psi)
  
  for (i in 1:nSites) {
    
    y[i] <- VGAM::rbetabinom(1, nVisits, z[i] * params$p, rho)
    
  }
  
  return(y)
}

### Model 6 - mixture on detection###

# params requires 3 values 
#   1. p1 - probability of detection for group 1
#   2. p2 - probability of detection for group 2
#   3. psi - occupancy probability
# pMix is the mixture probability between the two groups

# You can choose to make detection probability change by site or visit. I think
# changing by visit is unrealistic as I don't think the j'th visit to every site
# happens on the same day for them to have the same p but the code is here if 
# we want to explore it.

simulate_det_pMix <- function(params, nSites, nVisits, 
                              pMix, hetSource = "sites") {
  
  y <- matrix(NA, nrow = nSites, ncol = nVisits)
  
  z <- rbinom(nSites, 1, params$psi)
  
  if (hetSource == "sites")
    
    p <- matrix(c(params$p1, params$p2)[rbinom(nSites, 1, pMix) + 1], 
                nrow = nSites, ncol = nVisits)
  
  if (hetSource == "visits")
    
    p <- matrix(c(params$p1, params$p2)[rbinom(nVisits, 1, pMix) + 1], 
                nrow = nSites, ncol = nVisits, byrow = TRUE)
  
  for (i in 1:nSites) {
    
    for (j in 1:nVisits){
      
      y[i, j] <- rbinom(1, 1, z[i] * p[i, j])
      
    }
    
  }
  
  return(y)
}


### Model 7 - Outliers in detection probability
# params requires 2 values 
#   1. p - probability of detection
#   2. psi - occupancy probability
# beta_p models the diff in detection prob between the group and the outliers
# nOutliers defines how many sites are outliers  

simulate_outlier_det <- function(params, nSites, nVisits, beta_p, nOutliers) {
  
  y <- matrix(NA, nrow = nSites, ncol = nVisits)
  
  z <- rbinom(nSites, 1, params$psi)
  
  p <- rep(params$p, nSites)
  
  p[sample(1:nSites, nOutliers, replace = FALSE)] <- 1 / 
    (1 + exp(-(log(params$p / (1 - params$p)) + beta_p)))
  
  for (i in 1:nSites) {
    
    y[i, ] <- rbinom(nVisits, 1, z[i] * p[i])
    
  }
  
  return(y)
  
}

### Model 8 - mixture on occupancy###

# params requires 3 values 
#   1. psi1 - Occupancy of site in group 1
#   2. psi2 - Occupancy of site in group 2
#   3. p - detection probability
# pMix is the mixture probability between the two groups

simulate_occ_pMix <- function(params, nRegions, nSites, nVisits, pMix) {
  
  y <- matrix(NA, nrow = nSites, ncol = nVisits)
  
  psi <- rep(c(params$psi1, params$psi2)[rbinom(nRegions, 1, pMix) + 1],
             each = nSites / nRegions)
  
  z <- rbinom(nSites, 1, psi)
  
  for (i in 1:nSites) {
    
    y[i, ] <- rbinom(nVisits, 1, z[i] * params$p)
    
  }
  
  return(y)
  
}


### Model 9 - Outliers in Occupancy
# params requires 2 values 
#   1. p - probability of detection
#   2. psi - occupancy probability
# beta_o models the difference in occupancy between the group and the outliers
# nOutliers defines how many regions are outliers  

simulate_outlier_occ <- function(params, nRegions, nSites, 
                                 nVisits, beta_o, nOutliers) {
  
  y <- matrix(NA, nrow = nSites, ncol = nVisits)
  
  psi_region <- rep(params$psi, nRegions)
  
  psi_region[sample(1:nRegions, nOutliers, replace = FALSE)] <- 1 / 
    (1 + exp(-(log(params$psi / (1 - params$psi)) + beta_o)))
  
  psi <- rep(psi_region, each = nSites / nRegions)
  
  z <- rbinom(nSites, 1, psi)
  
  for (i in 1:nSites) {
    
    y[i, ] <- rbinom(nVisits, 1, z[i] * params$p)
    
  }
  
  return(y)
  
}
