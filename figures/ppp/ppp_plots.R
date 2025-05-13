library(tidyverse)
library(bayestestR)

# load auxiliary functions
source("utils.R")

################
# get coverage #
################

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

# detection mixture
coverage_detMix <- get_coverage(
  true_p = c(0.2, 0.8), true_psi = 0.6, 
  samples_p_fname = "posterior/alt/p_detectionMix.rds", 
  samples_psi_fname = "posterior/alt/psi_detectionMix.rds", axis = "pMix"
)

# occupancy interaction
coverage_interaction_beta0 <- get_coverage_beta(
  true_p = 0.1, true_beta = 0, 
  samples_p_fname = "posterior/alt/p_interaction_0.1.rds", 
  samples_beta_fname = "posterior/alt/beta0_interaction_0.1.rds", axis = "beta2")
coverage_interaction_beta1 <- get_coverage_beta(
  true_p = 0.1, true_beta = 0.5, 
  samples_p_fname = "posterior/alt/p_interaction_0.1.rds", 
  samples_beta_fname = "posterior/alt/beta1_interaction_0.1.rds", axis = "beta2")


# detection outlier
p <- 0.2
beta_p <- c(1, 3, 5)
pOutlier <- 0.05
p_alt <- 1 / (1 + exp(-(log(p / (1 - p)) + beta_p)))
p_multi <- p * (1 - pOutlier) + p_alt * pOutlier

coverage_detOut <- get_coverage(
  true_p = p_multi, true_psi = 0.6, 
  samples_p_fname = "posterior/alt/p_detectionOut.rds", 
  samples_psi_fname = "posterior/alt/psi_detectionOut.rds", axis = "beta_p"
)

# occupancy mixture
coverage_occMix <- get_coverage(
  true_p = 0.3, true_psi = c(0.2, 0.8), 
  samples_p_fname = "posterior/alt/p_occupancyMix.rds", 
  samples_psi_fname = "posterior/alt/psi_occupancyMix.rds", axis = "pMix"
)

# occupancy outlier
psi <- 0.3
beta_o <- c(0.1, 0.5, 1)
pOutlier <- 0.05
psi_alt <- 1 / (1 + exp(-(log(psi / (1 - psi)) + beta_o)))
psi_multi <- psi_alt * (1 - pOutlier) + psi_alt * pOutlier

coverage_occOut <- get_coverage(
  true_p = 0.2, true_psi = psi_multi, 
  samples_p_fname = "posterior/alt/p_occupancyOut.rds", 
  samples_psi_fname = "posterior/alt/psi_occupancyOut.rds", axis = "beta_o"
)


#############
#############
# ppp plots #
#############
#############

###############
# no covariates

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

###############
# occupancy covariates

# read in null hypothesis ppp
pvalues_null_notcond <- as.data.frame(readRDS("pvalues/null/interaction_0.1.rds")) %>% 
  pivot_longer(cols = everything(),
               names_to = "measure",
               values_to = "pvalue") %>% 
  filter(measure %in% c("deviance_latent", "chi_sq_latent", "lik_ratio_latent",
                        "ftukey_latent")) %>% 
  mutate(hyp = "null",
         measure = recode(measure,
                          "chi_sq_latent" = "chisq", "ftukey_latent" = "ftukey",
                          "deviance_latent" = "deviance",
                          "lik_ratio_latent" = "likratio"))
pvalues_null_cond <- as.data.frame(readRDS("pvalues/null/interaction_0.1.rds")) %>% 
  pivot_longer(cols = everything(),
               names_to = "measure",
               values_to = "pvalue") %>% 
  filter(measure %in% c("deviance", "chi_sq", "lik_ratio", "ftukey")) %>% 
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
# occupancy interaction #
#########################

##################
# prep figure data

pvalue_fig_inter <- prep_fig_data(pvalue_fname = "pvalues/alt/interaction_0.1.rds", 
                                   axis_name = "beta2", 
                                   axis_values = c(1, 3, 5), 
                                   null_df_notcond = pvalues_null_notcond,
                                   null_df_cond = pvalues_null_cond)

# create figure labels
beta2_labels <- c(
  "1" = "beta2 = 1",
  "3" = "beta2 = 3",
  "5" = "beta2 = 5"
)

# create plots
plots_inter <- get_ppp_plots(pvalue_fig_data = pvalue_fig_inter, 
                              title = "covariate interaction", 
                              axis = "beta2", labels = beta2_labels)

# save plots
ggsave("figures/ppp/pvalues_inter_notcond_p0.1.png", plots_inter[[1]],
       dpi = 400, width = 7, height = 5)
ggsave("figures/ppp/pvalues_inter_cond_p0.1.png", plots_inter[[2]],
       dpi = 400, width = 7, height = 5)

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


#####################
# outlier occupancy #
#####################

##################
# prep figure data

pvalue_fig_occOut <- prep_fig_data(
  pvalue_fname = "pvalues/alt/occupancyOut.rds",
  axis_name = "beta_o", axis_values = c(0.1, 0.5, 1),
  null_df_notcond = pvalues_null_notcond,
  null_df_cond = pvalues_null_cond)

# create figure labels
occOut_labels <- c(
  "0.1" = "beta_o = 0.1",
  "0.5" = "beta_o = 0.5",
  "1" = "beta_o = 1"
)

# create plots
plots_occOut <- get_ppp_plots(pvalue_fig_data = pvalue_fig_occOut, 
                              title = "occupancy outlier", 
                              axis = "beta_o", labels = occOut_labels)

# save plots
ggsave("figures/ppp/pvalues_occOut_notcond.png", plots_occOut[[1]],
       dpi = 400, width = 7, height = 5)
ggsave("figures/ppp/pvalues_occOut_cond.png", plots_occOut[[2]],
       dpi = 400, width = 7, height = 5)


