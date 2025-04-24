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

#########################
# site non-independence #
#########################

##################
# read in p values

# alt hypothesis
pvalues_nonind <- as.data.frame(readRDS("pvalues/alt/nonind.rds"))
colnames(pvalues_nonind) <- c("sd", "deviance_yes", "chisq_yes", "likratio_yes",
                              "ftukey_yes", "deviance_no", "chisq_no",
                              "likratio_no", "ftukey_no")
pvalues_nonind_long <- pvalues_nonind %>% 
  pivot_longer(cols = -sd, 
               names_to = c("measure", "latent_cond"),
               values_to = "pvalue",
               names_sep = "_") %>% 
  mutate(hyp = "alt")

# join null and alternative together
pvalues_nonind_all_notcond <- rbind(
  pvalues_nonind_long[pvalues_nonind_long$latent_cond == "no",
                      c("sd", "measure", "pvalue", "hyp")],
  pvalues_null_notcond %>% mutate(sd = 2.5),
  pvalues_null_notcond %>% mutate(sd = 5),
  pvalues_null_notcond %>% mutate(sd = 10))
pvalues_nonind_all_cond <- rbind(
  pvalues_nonind_long[pvalues_nonind_long$latent_cond == "yes",
                      c("sd", "measure", "pvalue", "hyp")],
  pvalues_null_cond %>% mutate(sd = 2.5),
  pvalues_null_cond %>% mutate(sd = 5),
  pvalues_null_cond %>% mutate(sd = 10))

# create figure labels
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

# create plots
plot_nonind_notcond <- ggplot(pvalues_nonind_all_notcond, 
                      aes(x = pvalue, fill = hyp)) +
  geom_histogram(position = "identity", alpha = 0.5) +
  scale_x_continuous(limits = c(0, 1),
                     breaks = c(0, 0.25, 0.5, 0.75, 1),
                     labels = c("0", "0.25", "0.5", "0.75", "1")) +
  facet_grid(measure ~ sd, labeller = labeller(
    measure = discrepancy_labels, sd = sd_labels)) +
  labs(x = "ppp", y = "count", fill = "hypothesis") +
  ggtitle("site non-independence, ppp not conditioned on latent state") +
  theme_minimal()

plot_nonind_cond <- ggplot(pvalues_nonind_all_cond, 
                              aes(x = pvalue, fill = hyp)) +
  geom_histogram(position = "identity", alpha = 0.5) +
  scale_x_continuous(limits = c(0, 1),
                     breaks = c(0, 0.25, 0.5, 0.75, 1),
                     labels = c("0", "0.25", "0.5", "0.75", "1")) +
  facet_grid(measure ~ sd, labeller = labeller(
    measure = discrepancy_labels, sd = sd_labels)) +
  labs(x = "ppp", y = "count", fill = "hypothesis") +
  ggtitle("site non-independence, ppp conditioned on latent state") +
  theme_minimal()

# save plots
ggsave("figures/pvalues_nonind_notcond.png", dpi = 400, plot_nonind_notcond,
       width = 7, height = 5)
ggsave("figures/pvalues_nonind_cond.png", dpi = 400, plot_nonind_cond,
       width = 7, height = 5)


###########################
# beta-binomial detection #
###########################

##################
# read in p values

# alt hypothesis
pvalues_betabin <- as.data.frame(readRDS("pvalues/alt/betabinom.rds"))
colnames(pvalues_betabin) <- c("rho", "deviance_yes", "chisq_yes", "likratio_yes",
                              "ftukey_yes", "deviance_no", "chisq_no",
                              "likratio_no", "ftukey_no")
pvalues_betabin_long <- pvalues_betabin %>% 
  pivot_longer(cols = -rho, 
               names_to = c("measure", "latent_cond"),
               values_to = "pvalue",
               names_sep = "_") %>% 
  mutate(hyp = "alt")

# join null and alternative together
pvalues_betabin_all_notcond <- rbind(
  pvalues_betabin_long[pvalues_betabin_long$latent_cond == "no",
                      c("rho", "measure", "pvalue", "hyp")],
  pvalues_null_notcond %>% mutate(rho = 0.1),
  pvalues_null_notcond %>% mutate(rho = 0.5),
  pvalues_null_notcond %>% mutate(rho = 0.9))
pvalues_betabin_all_cond <- rbind(
  pvalues_betabin_long[pvalues_betabin_long$latent_cond == "yes",
                      c("rho", "measure", "pvalue", "hyp")],
  pvalues_null_cond %>% mutate(rho = 0.1),
  pvalues_null_cond %>% mutate(rho = 0.5),
  pvalues_null_cond %>% mutate(rho = 0.9))

# create figure labels
rho_labels <- c(
  "0.1" = "rho = 0.1",
  "0.5" = "rho = 0.5",
  "0.9" = "rho = 0.9"
)

# create plots
plot_betabin_notcond <- ggplot(pvalues_betabin_all_notcond, 
                              aes(x = pvalue, fill = hyp)) +
  geom_histogram(position = "identity", alpha = 0.5) +
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1),
                     labels = c("0", "0.25", "0.5", "0.75", "1")) +
  facet_grid(measure ~ rho, labeller = labeller(
    measure = discrepancy_labels, rho = rho_labels)) +
  labs(x = "ppp", y = "count", fill = "hypothesis") +
  ggtitle("beta-binomial detection, ppp not conditioned on latent state") +
  theme_minimal()

plot_betabin_cond <- ggplot(pvalues_betabin_all_cond, 
                           aes(x = pvalue, fill = hyp)) +
  geom_histogram(position = "identity", alpha = 0.5) +
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1),
                     labels = c("0", "0.25", "0.5", "0.75", "1")) +
  facet_grid(measure ~ rho, labeller = labeller(
    measure = discrepancy_labels, rho = rho_labels)) +
  labs(x = "ppp", y = "count", fill = "hypothesis") +
  ggtitle("beta-binomial detection, ppp conditioned on latent state") +
  theme_minimal()

# save plots
ggsave("figures/pvalues_betabin_notcond.png", dpi = 400, plot_betabin_notcond,
       width = 7, height = 5)
ggsave("figures/pvalues_betabin_cond.png", dpi = 400, plot_betabin_cond,
       width = 7, height = 5)


