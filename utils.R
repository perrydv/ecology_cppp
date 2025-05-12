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


# function for getting coverage
get_coverage <- function(true_p, true_psi, samples_p_fname, 
                         samples_psi_fname, axis) {
  
  # read in samples
  samples_psi <- readRDS(samples_psi_fname)
  samples_p <- readRDS(samples_p_fname)
  
  # create empty df to store coverage
  coverage <- matrix(NA, nrow = nrow(samples_psi), ncol = 3)
  colnames(coverage) <- c("coverage_p", "coverage_psi", axis)
  
  # get coverage
  for (i in 1:nrow(samples_psi)) {
    
    # calculate credible interval
    ci_p <- hdi(samples_p[i, 3:ncol(samples_p)])[2:3]
    ci_psi <- hdi(samples_psi[i, 3:ncol(samples_psi)])[2:3]
    
    # get p coverage
    if (length(true_p) == 1) {
      coverage[i, "coverage_p"] <- true_p >= ci_p$CI_low && 
        true_p <= ci_p$CI_high 
    } else if (length(true_p) == 2) {
      p <- true_p[1] * (1 - samples_p[i, "pMix"]) + 
        true_p[2] * samples_p[i, "pMix"]
      coverage[i, "coverage_p"] <- p >= ci_p$CI_low && p <= ci_p$CI_high 
    } else {
      p <- true_p[which(unique(samples_p[,"beta_p"]) == samples_p[i, "beta_p"])] 
      coverage[i, "coverage_p"] <- p >= ci_p$CI_low && p <= ci_p$CI_high 
    }
    
    # get psi coverage
    if (length(true_psi) == 1) {
      coverage[i, "coverage_psi"] <- true_psi >= ci_psi$CI_low && 
        true_psi <= ci_psi$CI_high
    } else if (length(true_psi) == 2) {
      psi <- true_psi[1] * (1 - samples_psi[i, "pMix"]) + 
        true_psi[2] * samples_psi[i, "pMix"]
      coverage[i, "coverage_psi"] <- psi >= ci_psi$CI_low && 
        psi <= ci_psi$CI_high 
    } else {
      psi <- true_psi[
        which(unique(samples_psi[,"beta_o"]) == samples_psi[i, "beta_o"])] 
      coverage[i, "coverage_psi"] <- psi >= ci_psi$CI_low && 
        psi <= ci_psi$CI_high 
    }
    
    # store axis
    coverage[i, axis] <- samples_p[i, axis]
    
  }
  
  # get coverage summary
  coverage_summary <- as.data.frame(coverage) %>% 
    group_by(.data[[axis]]) %>% 
    summarize(mean_p = mean(coverage_p),
              mean_psi = mean(coverage_psi))
  
  return(coverage_summary)
}

# function for getting coverage - site occupancy covariates
get_coverage_beta <- function(true_p, true_beta, samples_p_fname, 
                              samples_beta_fname, axis) {
  
  # read in samples
  samples_beta <- readRDS(samples_beta_fname)
  samples_p <- readRDS(samples_p_fname)
  
  # create empty df to store coverage
  coverage <- matrix(NA, nrow = nrow(samples_beta), ncol = 3)
  colnames(coverage) <- c("coverage_p", "coverage_beta", axis)
  
  # get coverage
  for (i in 1:nrow(samples_beta)) {
    
    # calculate credible interval
    ci_p <- hdi(samples_p[i, 3:ncol(samples_p)])[2:3]
    ci_beta <- hdi(samples_beta[i, 3:ncol(samples_beta)])[2:3]
    
    # get p coverage
    coverage[i, "coverage_p"] <- true_p >= ci_p$CI_low && 
      true_p <= ci_p$CI_high 
    
    # get beta coverage
    coverage[i, "coverage_beta"] <- true_beta >= ci_beta$CI_low && 
      true_beta <= ci_beta$CI_high
    
    # store axis
    coverage[i, axis] <- samples_p[i, axis]
    
  }
  
  # get coverage summary
  coverage_summary <- as.data.frame(coverage) %>% 
    group_by(.data[[axis]]) %>% 
    summarize(mean_p = mean(coverage_p),
              mean_beta = mean(coverage_beta))
  
  return(coverage_summary)
}


# function to prepare figure data
prep_fig_data <- function(pvalue_fname, axis_name, axis_values, 
                          null_df_notcond, null_df_cond) {
  
  # alt hypothesis
  pvalues <- as.data.frame(readRDS(pvalue_fname))
  colnames(pvalues) <- c(axis_name, "deviance_yes", "chisq_yes", "likratio_yes",
                         "ftukey_yes", "deviance_no", "chisq_no",
                         "likratio_no", "ftukey_no")
  pvalues_long <- pvalues %>% 
    pivot_longer(cols = -axis_name, 
                 names_to = c("measure", "latent_cond"),
                 values_to = "pvalue",
                 names_sep = "_") %>% 
    mutate(hyp = "alt")
  
  # join null and alternative together
  pvalues_all_notcond <- rbind(
    pvalues_long[pvalues_long$latent_cond == "no",
                 c(axis_name, "measure", "pvalue", "hyp")],
    null_df_notcond %>% mutate(!!axis_name := axis_values[1]),
    null_df_notcond %>% mutate(!!axis_name := axis_values[2]),
    null_df_notcond %>% mutate(!!axis_name := axis_values[3])
  )
  pvalues_all_cond <- rbind(
    pvalues_long[pvalues_long$latent_cond == "yes",
                 c(axis_name, "measure", "pvalue", "hyp")],
    null_df_cond %>% mutate(!!axis_name := axis_values[1]),
    null_df_cond %>% mutate(!!axis_name := axis_values[2]),
    null_df_cond %>% mutate(!!axis_name := axis_values[3])
  )
  
  return(list(pvalues_all_notcond, pvalues_all_cond))
}


# function for creating null/alt p value plots
get_ppp_plots <- function(pvalue_fig_data, title, axis, labels) {
  
  # build the facet formula dynamically
  facet_formula <- as.formula(paste("measure ~", axis))
  
  # create plots
  plot_notcond <- ggplot(pvalue_fig_data[[1]], 
                         aes(x = pvalue, fill = hyp)) +
    geom_histogram(position = "identity", alpha = 0.5) +
    scale_x_continuous(limits = c(-0.05, 1), 
                       breaks = c(0, 0.25, 0.5, 0.75, 1),
                       labels = c("0", "0.25", "0.5", "0.75", "1")) +
    facet_grid(facet_formula, labeller = labeller(
      measure = discrepancy_labels,
      .cols = labels)) +
    labs(x = "ppp", y = "count", fill = "hypothesis") +
    ggtitle(paste0(title, ", ppp not conditioned on latent state")) +
    theme_minimal()
  
  plot_cond <- ggplot(pvalue_fig_data[[2]], 
                      aes(x = pvalue, fill = hyp)) +
    geom_histogram(position = "identity", alpha = 0.5) +
    scale_x_continuous(limits = c(-0.05, 1),
                       breaks = c(0, 0.25, 0.5, 0.75, 1),
                       labels = c("0", "0.25", "0.5", "0.75", "1")) +
    facet_grid(facet_formula, labeller = labeller(
      measure = discrepancy_labels,
      .cols = labels)) +
    labs(x = "ppp", y = "count", fill = "hypothesis") +
    ggtitle(paste0(title, ", ppp conditioned on latent state")) +
    theme_minimal()
  
  return(list(plot_notcond, plot_cond))
}
