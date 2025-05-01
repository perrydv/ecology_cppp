library(tidyverse)
library(bayestestR)

# load auxiliary functions
source("utils.R")

################
# get coverage #
################

p <- 0.3
psi <- 0.6

# non-independent sites
coverage_nonind <- get_coverage(
  true_p = 0.3, true_psi = 0.6, samples_p_fname = "posterior/alt/p_nonind.rds", 
  samples_psi_fname = "posterior/alt/psi_nonind.rds", axis = "sd"
)

# beta-binomial detection
coverage_betabin <- get_coverage(
  true_p = 0.3, true_psi = 0.6, 
  samples_p_fname = "posterior/alt/p_betabinom.rds", 
  samples_psi_fname = "posterior/alt/psi_betabinom.rds", axis = "rho"
)


#############
#############
# ppp plots #
#############
#############

# read in null hypothesis ppp
pvalues_null_notcond <- as.data.frame(readRDS("pvalues/null/basic_latent.rds")) %>% 
  pivot_longer(cols = everything(),
               names_to = "measure",
               values_to = "pvalue") %>% 
  mutate(hyp = "null",
         measure = recode(measure,
                          "chi_sq" = "chisq", "lik_ratio" = "likratio"))
pvalues_null_cond <- as.data.frame(readRDS("pvalues/null/basic.rds")) %>% 
  pivot_longer(cols = everything(),
               names_to = "measure",
               values_to = "pvalue") %>% 
  mutate(hyp = "null",
         measure = recode(measure,
                          "chi_sq" = "chisq", "lik_ratio" = "likratio"))

# create figure labels
discrepancy_labels <- c(
  "chisq" = "chi-squared",
  "deviance" = "deviance",
  "ftukey" = "Freeman-Tukey",
  "likratio" = "likelihood\nratio"
)

#########################
# site non-independence #
#########################

##################
# prep figure data

pvalue_fig_nonind <- prep_fig_data(pvalue_fname = "pvalues/alt/nonind.rds", 
                                   axis_name = "sd", 
                                   axis_values = c(2.5, 5, 10), 
                                   null_df_notcond = pvalues_null_notcond,
                                   null_df_cond = pvalues_null_cond)

# create figure labels
sd_labels <- c(
  "2.5" = "sd = 2.5",
  "5" = "sd = 5",
  "10" = "sd = 10"
)

# create plots
plots_nonind <- get_ppp_plots(pvalue_fig_data = pvalue_fig_nonind, 
                               title = "site non-independence", 
                               axis = "sd", labels = sd_labels)

# save plots
ggsave("figures/ppp/pvalues_nonind_notcond.png", plots_nonind[[1]],
       dpi = 400, width = 7, height = 5)
ggsave("figures/ppp/pvalues_nonind_cond.png", plots_nonind[[2]],
       dpi = 400, width = 7, height = 5)


###########################
# beta-binomial detection #
###########################

##################
# prep figure data

pvalue_fig_betabin <- prep_fig_data(pvalue_fname = "pvalues/alt/betabinom.rds", 
                                    axis_name = "rho", 
                                    axis_values = c(0.1, 0.5, 0.9), 
                                    null_df_notcond = pvalues_null_notcond,
                                    null_df_cond = pvalues_null_cond)
# create figure labels
rho_labels <- c(
  "0.1" = "rho = 0.1",
  "0.5" = "rho = 0.5",
  "0.9" = "rho = 0.9"
)

# create plots
plots_betabin <- get_ppp_plots(pvalue_fig_data = pvalue_fig_betabin, 
                               title = "beta-binomial detection", 
                               axis = "rho", labels = rho_labels)

# save plots
ggsave("figures/ppp/pvalues_betabin_notcond.png", plots_betabin[[1]],
       dpi = 400, width = 7, height = 5)
ggsave("figures/ppp/pvalues_betabin_cond.png", plots_betabin[[2]],
       dpi = 400, width = 7, height = 5)


#####################
# detection mixture #
#####################

##################
# prep figure data

pvalue_fig_detMix <- prep_fig_data(
  pvalue_fname = "pvalues/alt/detectionMix.rds",
  axis_name = "pMix", axis_values = c(0.1, 0.2, 0.5),
  null_df_notcond = pvalues_null_notcond,
  null_df_cond = pvalues_null_cond)

# create figure labels
pMixdet_labels <- c(
  "0.1" = "pMix = 0.1",
  "0.2" = "pMix = 0.2",
  "0.5" = "pMix = 0.5"
)

# create plots
plots_detMix <- get_ppp_plots(pvalue_fig_data = pvalue_fig_detMix, 
                               title = "detection mixture", 
                               axis = "pMix", labels = pMixdet_labels)

# save plots
ggsave("figures/ppp/pvalues_detMix_notcond.png", plots_detMix[[1]],
       dpi = 400, width = 7, height = 5)
ggsave("figures/ppp/pvalues_detMix_cond.png", plots_detMix[[2]],
       dpi = 400, width = 7, height = 5)


#####################
# detection outlier #
#####################

##################
# prep figure data

pvalue_fig_detOut <- prep_fig_data(
  pvalue_fname = "pvalues/alt/detectionOut.rds",
  axis_name = "beta_p", axis_values = c(1, 3, 5),
  null_df_notcond = pvalues_null_notcond,
  null_df_cond = pvalues_null_cond)

# create figure labels
detOut_labels <- c(
  "1" = "beta_p = 1",
  "3" = "beta_p = 3",
  "5" = "beta_p = 5"
)

# create plots
plots_detOut <- get_ppp_plots(pvalue_fig_data = pvalue_fig_detOut, 
                              title = "detection outlier", 
                              axis = "beta_p", labels = detOut_labels)

# save plots
ggsave("figures/ppp/pvalues_detOut_notcond.png", plots_detOut[[1]],
       dpi = 400, width = 7, height = 5)
ggsave("figures/ppp/pvalues_detOut_cond.png", plots_detOut[[2]],
       dpi = 400, width = 7, height = 5)


#####################
# mixture occupancy #
#####################

##################
# prep figure data

pvalue_fig_occMix <- prep_fig_data(
  pvalue_fname = "pvalues/alt/occupancyMix.rds",
  axis_name = "pMix", axis_values = c(0.01, 0.2, 0.5),
  null_df_notcond = pvalues_null_notcond,
  null_df_cond = pvalues_null_cond)

# create figure labels
pMixocc_labels <- c(
  "0.01" = "pMix = 0.01",
  "0.2" = "pMix = 0.2",
  "0.5" = "pMix = 0.5"
)

# create plots
plots_occMix <- get_ppp_plots(pvalue_fig_data = pvalue_fig_occMix, 
                              title = "occupancy mixture", 
                              axis = "pMix", labels = pMixocc_labels)

# save plots
ggsave("figures/ppp/pvalues_occMix_notcond.png", plots_occMix[[1]],
       dpi = 400, width = 7, height = 5)
ggsave("figures/ppp/pvalues_occMix_cond.png", plots_occMix[[2]],
       dpi = 400, width = 7, height = 5)


