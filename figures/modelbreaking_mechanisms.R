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
ggsave("figures/modelbreak_mechanism_nonind.png", dpi = 400, 
       width = 7, height = 3)

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
ggsave("figures/modelbreak_mechanism_betabin.png", dpi = 400, 
       width = 7, height = 3)


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
ggsave("figures/modelbreak_mechanism_mixture_detection.png", 
       final_mixturedet_plot, dpi = 400, 
       width = 7, height = 3)



