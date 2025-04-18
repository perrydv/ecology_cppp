library(tidyverse)
library(bayestestR)

################
# get coverage #
################

p <- 0.3
psi <- 0.6

# read in samples
samples_psi_nonind <- readRDS("posterior/alt/psi_nonind.rds")
samples_p_nonind <- readRDS("posterior/alt/p_nonind.rds")
samples_psi_betabinom <- readRDS("posterior/alt/psi_betabinom.rds")
samples_p_betabinom <- readRDS("posterior/alt/p_betabinom.rds")

# create empty df to store coverage
coverage_nonind <- matrix(NA, nrow = nrow(samples_psi_nonind), ncol = 3)
coverage_betabin <- matrix(NA, nrow = nrow(samples_psi_betabinom), ncol = 3)
colnames(coverage_nonind) <- c("coverage_p", "coverage_psi", "sd")
colnames(coverage_betabin) <- c("coverage_p", "coverage_psi", "rho")

# get coverage
for (i in 1:nrow(samples_psi_nonind)) {
  
  # non ind
  ci_nonind_p <- hdi(samples_p_nonind[i, 3:802])[2:3]
  ci_nonind_psi <- hdi(samples_psi_nonind[i, 3:802])[2:3]
  coverage_nonind[i, "coverage_p"] <- p >= ci_nonind_p$CI_low && 
    p <= ci_nonind_p$CI_high
  coverage_nonind[i, "coverage_psi"] <- psi >= ci_nonind_psi$CI_low && 
    psi <= ci_nonind_psi$CI_high
  coverage_nonind[i, "sd"] <- samples_p_nonind[i, "sd"]
  
  # betabinom
  ci_betabinom_p <- hdi(samples_p_betabinom[i, 3:802])[2:3]
  ci_betabinom_psi <- hdi(samples_psi_betabinom[i, 3:802])[2:3]
  coverage_betabin[i, "coverage_p"] <- p >= ci_betabinom_p$CI_low && 
    p <= ci_betabinom_p$CI_high
  coverage_betabin[i, "coverage_psi"] <- psi >= ci_betabinom_psi$CI_low && 
    psi <= ci_betabinom_psi$CI_high
  coverage_betabin[i, "rho"] <- samples_p_betabinom[i, "rho"]
  
}

# get coverage summaries
coverage_nonind_summary <- as.data.frame(coverage_nonind) %>% 
  group_by(sd) %>% 
  summarize(mean_p = mean(coverage_p),
            mean_psi = mean(coverage_psi))

coverage_betabin_summary <- as.data.frame(coverage_betabin) %>% 
  group_by(rho) %>% 
  summarize(mean_p = mean(coverage_p),
            mean_psi = mean(coverage_psi))

#########################
# site non-independence #
#########################

# read in p values
pvalues_nonind <- as.data.frame(readRDS("pvalues/alt/nonind.rds"))
colnames(pvalues_nonind) <- c("sd", "deviance_yes", "chisq_yes", "likratio_yes",
                              "ftukey_yes", "deviance_no", "chisq_no",
                              "likratio_no", "ftukey_no")
pvalues_nonind_long <- pvalues_nonind %>% 
  pivot_longer(cols = -sd, 
               names_to = c("measure", "latent_cond"),
               values_to = "pvalue",
               names_sep = "_")

discrepancy_labels <- c(
  "chisq" = "chi-squared",
  "deviance" = "deviance",
  "ftukey" = "Freeman-Tukey",
  "likratio" = "likelihood ratio"
)

sd_labels <- c(
  "2.5" = "sd = 2.5",
  "5" = "sd = 5",
  "10" = "sd = 10"
)

fill_labels <- c(
  "no" = expression(y^{rep} ~ "|" ~ theta^L * "," ~ theta^T),
  "yes" = expression(y^{rep} * theta^{L~","~rep} ~ "|" ~ theta^T)
)

plot_nonind <- ggplot(pvalues_nonind_long, 
                      aes(x = pvalue, fill = latent_cond)) +
  geom_histogram(position = "identity", alpha = 0.5) +
  scale_x_continuous(limits = c(0, 1),
                     breaks = c(0, 0.25, 0.5, 0.75, 1),
                     labels = c("0", "0.25", "0.5", "0.75", "1")) +
  scale_fill_manual(
    values = c("no" = "#1f77b4",
               "yes" = "#ff7f0e"),
    labels = fill_labels) +
  facet_grid(measure ~ sd, labeller = labeller(
    measure = discrepancy_labels, sd = sd_labels)) +
  labs(x = "p value", y = "count", fill = "latent") +
  theme_minimal()
ggsave("figures/alt_pvalues_nonind.png", dpi = 400, plot_nonind,
       width = 7, height = 5)


###########################
# beta-binomial detection #
###########################

# read in p values
pvalues_betabin <- as.data.frame(readRDS("pvalues/alt/betabinom.rds"))
colnames(pvalues_betabin) <- c("rho", "deviance_yes", "chisq_yes", "likratio_yes",
                              "ftukey_yes", "deviance_no", "chisq_no",
                              "likratio_no", "ftukey_no")
pvalues_betabin_long <- pvalues_betabin %>% 
  pivot_longer(cols = -rho, 
               names_to = c("measure", "latent_cond"),
               values_to = "pvalue",
               names_sep = "_")

rho_labels <- c(
  "0.1" = "rho = 0.1",
  "0.5" = "rho = 0.5",
  "0.9" = "rho = 0.9"
)

plot_betabinom <- ggplot(pvalues_betabin_long, 
                      aes(x = pvalue, fill = latent_cond)) +
  geom_histogram(position = "identity", alpha = 0.5) +
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1),
                     labels = c("0", "0.25", "0.5", "0.75", "1")) +
  scale_fill_manual(
    values = c("no" = "#1f77b4",
               "yes" = "#ff7f0e"),
    labels = fill_labels) +
  facet_grid(measure ~ rho, labeller = labeller(
    measure = discrepancy_labels, rho = rho_labels)) +
  labs(x = "p value", y = "count", fill = "latent") +
  theme_minimal()
ggsave("figures/alt_pvalues_betabinom.png", dpi = 400, plot_betabinom,
       width = 7, height = 5)


