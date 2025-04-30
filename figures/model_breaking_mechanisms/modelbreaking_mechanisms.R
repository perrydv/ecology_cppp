library(tidyverse)
library(patchwork)

# non independent sites
psi <- 0.6

sd_2.5 <- plogis(rnorm(10000, psi, 2.5))
sd_5 <- plogis(rnorm(10000, psi, 5))
sd_10 <- plogis(rnorm(10000, psi, 10))

plot_sd2.5 <- ggplot() +
  geom_histogram(aes(x = sd_2.5)) +
  geom_vline(aes(xintercept = psi), color = "red", linetype = "dashed") +
  labs(x = "region psi", y = "count") +
  ggtitle("sd = 2.5") +
  theme_minimal()

plot_sd5 <- ggplot() +
  geom_histogram(aes(x = sd_5)) +
  geom_vline(aes(xintercept = psi), color = "red", linetype = "dashed") +
  labs(x = "region psi", y = "count") +
  ggtitle("sd = 5") +
  theme_minimal()

plot_sd10 <- ggplot() +
  geom_histogram(aes(x = sd_10)) +
  geom_vline(aes(xintercept = psi), color = "red", linetype = "dashed") +
  labs(x = "region psi", y = "count") +
  ggtitle("sd = 10") +
  theme_minimal()

final_sd_plot <- plot_sd2.5 + plot_sd5 + plot_sd10 + plot_layout(nrow = 1)
ggsave("figures/model_breaking_mechanisms/modelbreak_nonind.png", 
       dpi = 400, width = 7, height = 3)

# beta binomial detection
p <- 0.3
nVisits <- 50

rho_0.1 <- VGAM::rbetabinom(10000, nVisits, p, 0.1)
rho_0.5 <- VGAM::rbetabinom(10000, nVisits, p, 0.5)
rho_0.9 <- VGAM::rbetabinom(10000, nVisits, p, 0.9)

plot_rho0.1 <- ggplot() +
  geom_histogram(aes(x = rho_0.1)) +
  geom_vline(aes(xintercept = p * nVisits), 
             color = "red", linetype = "dashed") +
  labs(x = "detections / 50 visits", y = "count") +
  ggtitle("rho = 0.1") +
  theme_minimal()

plot_rho0.5 <- ggplot() +
  geom_histogram(aes(x = rho_0.5)) +
  geom_vline(aes(xintercept = p * nVisits), 
             color = "red", linetype = "dashed") +
  labs(x = "detections / 50 visits", y = "count") +
  ggtitle("rho = 0.5") +
  theme_minimal()

plot_rho0.9 <- ggplot() +
  geom_histogram(aes(x = rho_0.9)) +
  geom_vline(aes(xintercept = p * nVisits), 
             color = "red", linetype = "dashed") +
  labs(x = "detections / 50 visits", y = "count") +
  ggtitle("rho = 0.9") +
  theme_minimal()

final_rho_plot <- plot_rho0.1 + plot_rho0.5 + plot_rho0.9 + 
  plot_layout(nrow = 1)
ggsave("figures/model_breaking_mechanisms/modelbreak_betabin.png", 
       dpi = 400, width = 7, height = 3)


# mixture on detection (heterogeneity across sites)
p1 <- 0.2
p2 <- 0.8
pMix <- c(0.1, 0.2, 0.5)
nSites <- 1000
nVisits <- 50

pvec0.1 <- c(p1, p2)[rbinom(nSites, 1, pMix[1]) + 1]
pvec0.2 <- c(p1, p2)[rbinom(nSites, 1, pMix[2]) + 1]
pvec0.5 <- c(p1, p2)[rbinom(nSites, 1, pMix[3]) + 1]

plot_pMix0.1 <- ggplot() +
  geom_histogram(aes(x = sapply(pvec0.1, function(p) rbinom(1, nVisits, p)))) +
  labs(x = "detections / 50 visits", y = "count") +
  ggtitle("pMix = 0.1") +
  theme_minimal()

plot_pMix0.2 <- ggplot() +
  geom_histogram(aes(x = sapply(pvec0.2, function(p) rbinom(1, nVisits, p)))) +
  labs(x = "detections / 50 visits", y = "count") +
  ggtitle("pMix = 0.2") +
  theme_minimal()

plot_pMix0.5 <- ggplot() +
  geom_histogram(aes(x = sapply(pvec0.5, function(p) rbinom(1, nVisits, p)))) +
  labs(x = "detections / 50 visits", y = "count") +
  ggtitle("pMix = 0.5") +
  theme_minimal()

final_mixturedet_plot <- plot_pMix0.1 + plot_pMix0.2 + plot_pMix0.5 + 
  plot_layout(nrow = 1)
ggsave("figures/model_breaking_mechanisms/modelbreak_mixture_detection.png", 
       final_mixturedet_plot, dpi = 400, 
       width = 7, height = 3)

# detection outliers
p <- 0.2
beta_p <- c(1, 3, 5)
nSites <- 1000
nOutliers <- 0.05 * nSites
nVisits <- 50

p_vec1 <- p_vec3 <- p_vec5 <-rep(p, nSites)

p_vec1[sample(1:nSites, nOutliers, replace = F)] <- 1 / 
  (1 + exp(-(log(p / (1 - p)) + beta_p[1])))
p_vec3[sample(1:nSites, nOutliers, replace = F)] <- 1 / 
  (1 + exp(-(log(p / (1 - p)) + beta_p[2])))
p_vec5[sample(1:nSites, nOutliers, replace = F)] <- 1 / 
  (1 + exp(-(log(p / (1 - p)) + beta_p[3])))

plot_pOut1 <- ggplot() +
  geom_histogram(aes(x = sapply(p_vec1, function(p) rbinom(1, nVisits, p)))) +
  labs(x = "detections / 50 visits", y = "count") +
  ggtitle("beta = 1") +
  theme_minimal()

plot_pOut3 <- ggplot() +
  geom_histogram(aes(x = sapply(p_vec3, function(p) rbinom(1, nVisits, p)))) +
  labs(x = "detections / 50 visits", y = "count") +
  ggtitle("beta = 3") +
  theme_minimal()

plot_pOut5 <- ggplot() +
  geom_histogram(aes(x = sapply(p_vec5, function(p) rbinom(1, nVisits, p)))) +
  labs(x = "detections / 50 visits", y = "count") +
  ggtitle("beta = 5") +
  theme_minimal()

final_outdet_plot <- plot_pOut1 + plot_pOut3 + plot_pOut5 + 
  plot_layout(nrow = 1)
ggsave("figures/model_breaking_mechanisms/modelbreak_outlier_detection.png", 
       final_outdet_plot, dpi = 400, 
       width = 7, height = 3)

# mixture on occupancy
psi1 <- 0.2
psi2 <- 0.8
pMix <- c(0.01, 0.2, 0.5)
nSites <- 1000
nRegions <- nSites / 5
nVisits <- 50

psivec0.01 <- rbinom(nSites, 1, 
                    p = rep(c(psi1, psi2)[rbinom(nRegions, 1, pMix[1]) + 1],
                            each = nSites / nRegions))
psivec0.2 <- rbinom(nSites, 1, 
                    p = rep(c(psi1, psi2)[rbinom(nRegions, 1, pMix[2]) + 1],
                            each = nSites / nRegions))
psivec0.5 <- rbinom(nSites, 1, 
                    p = rep(c(psi1, psi2)[rbinom(nRegions, 1, pMix[3]) + 1],
                            each = nSites / nRegions))

plot_oMix0.01 <- ggplot() +
  geom_histogram(aes(x = psivec0.01)) +
  labs(x = "z", y = "count") +
  ggtitle("pMix = 0.01") +
  theme_minimal()

plot_oMix0.2 <- ggplot() +
  geom_histogram(aes(x = psivec0.2)) +
  labs(x = "z", y = "count") +
  ggtitle("pMix = 0.2") +
  theme_minimal()

plot_oMix0.5 <- ggplot() +
  geom_histogram(aes(x = psivec0.5)) +
  labs(x = "z", y = "count") +
  ggtitle("pMix = 0.5") +
  theme_minimal()

final_mixtureocc_plot <- plot_oMix0.01 + plot_oMix0.2 + plot_oMix0.5 + 
  plot_layout(nrow = 1)
ggsave("figures/model_breaking_mechanisms/modelbreak_mixture_occ.png", 
       final_mixtureocc_plot, dpi = 400, 
       width = 7, height = 3)


# occupancy outliers
p <- 0.2
psi <- 0.3
beta_o <- c(0.1, 0.5, 1)
nSites <- 1000
nRegions <- nSites / 5
nOutliers <- 0.05 * nRegions
nVisits <- 50

psi_region <- rep(psi, nRegions)

psi_vec0.1 <- rep(psi_region[sample(1:nRegions, nOutliers, replace = F)] <- 1 / 
                  (1 + exp(-(log(psi / (1 - psi)) + beta_o[1]))), 
                each = nSites/nRegions)

psi_vec0.5 <- rep(psi_region[sample(1:nRegions, nOutliers, replace = F)] <- 1 / 
                  (1 + exp(-(log(psi / (1 - psi)) + beta_o[2]))), 
                each = nSites/nRegions)

psi_vec1 <- rep(psi_region[sample(1:nRegions, nOutliers, replace = F)] <- 1 / 
                  (1 + exp(-(log(psi / (1 - psi)) + beta_o[3]))), 
                each = nSites/nRegions)

plot_oOut0.1 <- ggplot() +
  geom_histogram(aes(x = rbinom(nSites, 1, psi_vec0.1))) +
  labs(x = "z", y = "count") +
  ggtitle("beta = 0.1") +
  theme_minimal()

plot_oOut0.5 <- ggplot() +
  geom_histogram(aes(x = rbinom(nSites, 1, psi_vec0.5))) +
  labs(x = "z", y = "count") +
  ggtitle("beta = 0.5") +
  theme_minimal()

plot_oOut1 <- ggplot() +
  geom_histogram(aes(x = rbinom(nSites, 1, psi_vec1))) +
  labs(x = "z", y = "count") +
  ggtitle("beta = 1") +
  theme_minimal()

final_outlierocc_plot <- plot_oOut0.1 + plot_oOut0.5 + plot_oOut1 + 
  plot_layout(nrow = 1)
ggsave("figures/model_breaking_mechanisms/modelbreak_outlier_occ.png", 
       final_outlierocc_plot, dpi = 400, 
       width = 7, height = 3)

