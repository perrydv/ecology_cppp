# function for calculating p values

calc_pvalues <- function(n_measures, samples) {
  
  # create empty vector 
  vec <- rep(NA, n_measures)
  
  # calculate p values - conditioned on latent state
  vec[1] <- mean(samples[, "D_rep_total"] > samples[, "D_obs_total"])
  vec[2] <- mean(samples[, "chi_rep_total"] > samples[, "chi_obs_total"])
  vec[3] <- mean(samples[, "ratio_rep_total"] > samples[, "ratio_obs_total"])
  vec[4] <- mean(samples[, "tukey_rep_total"] > samples[, "tukey_obs_total"])
  
  # calculate p values - not conditioned on latent state
  vec[5] <- mean(samples[, "D_rep_latent_total"] > samples[, "D_obs_total"])
  vec[6] <- mean(samples[, "chi_rep_latent_total"] > samples[, "chi_obs_total"])
  vec[7] <- mean(
    samples[, "ratio_rep_latent_total"] > samples[, "ratio_obs_total"]
  )
  vec[8] <- mean(
    samples[, "tukey_rep_latent_total"] > samples[, "tukey_obs_total"]
  )
  
  return(vec)
}
