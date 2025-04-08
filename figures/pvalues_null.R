library(tidyverse)

# read in p values
pvalues <- rbind(
  as.data.frame(readRDS("pvalues/null/basic.rds")) %>% 
    pivot_longer(cols = everything(), names_to = "discrepancy",
                 values_to = "p_value") %>% 
    mutate(model = "basic", latent = "no"),
  as.data.frame(readRDS("pvalues/null/cov_occ.rds")) %>% 
    pivot_longer(cols = everything(), names_to = "discrepancy",
                 values_to = "p_value") %>% 
    mutate(model = "cov_occ", latent = "no"),
  as.data.frame(readRDS("pvalues/null/cov_occ_det.rds")) %>% 
    pivot_longer(cols = everything(), names_to = "discrepancy",
                 values_to = "p_value") %>% 
    mutate(model = "cov_occ_det", latent = "no"),
  as.data.frame(readRDS("pvalues/null/sp_ranef.rds")) %>% 
    pivot_longer(cols = everything(), names_to = "discrepancy",
                 values_to = "p_value") %>% 
    mutate(model = "sp_ranef", latent = "no"),
  as.data.frame(readRDS("pvalues/null/basic_latent.rds")) %>% 
    pivot_longer(cols = everything(), names_to = "discrepancy",
                 values_to = "p_value") %>% 
    mutate(model = "basic", latent = "yes"),
  as.data.frame(readRDS("pvalues/null/cov_occ_latent.rds")) %>% 
    pivot_longer(cols = everything(), names_to = "discrepancy",
                 values_to = "p_value") %>% 
    mutate(model = "cov_occ", latent = "yes"),
  as.data.frame(readRDS("pvalues/null/cov_occ_det_latent.rds")) %>% 
    pivot_longer(cols = everything(), names_to = "discrepancy",
                 values_to = "p_value") %>% 
    mutate(model = "cov_occ_det", latent = "yes"),
  as.data.frame(readRDS("pvalues/null/sp_ranef_latent.rds")) %>% 
    pivot_longer(cols = everything(), names_to = "discrepancy",
                 values_to = "p_value") %>% 
    mutate(model = "sp_ranef", latent = "yes")
)

model_labels <- c(
  "basic" = "basic",
  "cov_occ" = "occupancy\ncovariates",
  "cov_occ_det" = "occupancy +\ndetection\ncovariates",
  "sp_ranef" = "site random effect"
)

discrepancy_labels <- c(
  "chi_sq" = "chi-squared",
  "deviance" = "deviance",
  "ftukey" = "Freeman-Tukey",
  "lik_ratio" = "likelihood ratio"
)

fill_labels <- c(
  "no" = expression(y^{rep} ~ "|" ~ theta^L * "," ~ theta^T),
  "yes" = expression(y^{rep} * theta^{L~","~rep} ~ "|" ~ theta^T)
)


plot_wlatent <- ggplot(pvalues, aes(x = p_value, fill = latent)) +
  geom_histogram() +
  scale_y_continuous(breaks = c(0, 50, 100)) +
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1),
                     labels = c("0", "0.25", "0.5", "0.75", "1")) +
  scale_fill_manual(
    values = c("no" = "#1f77b4",
               "yes" = "#ff7f0e"),
    labels = fill_labels
  ) +
  facet_grid(model ~ discrepancy, labeller = labeller(
    model = model_labels,
    discrepancy = discrepancy_labels
  )) +
  labs(x = "p value", y = "count") +
  theme_minimal()
ggsave("figures/null_pvalues_wlatent.png", dpi = 400, plot_wlatent,
       width = 7, height = 5)

plot <- ggplot(pvalues[pvalues$latent == "no",], 
               aes(x = p_value)) +
  geom_histogram() +
  scale_y_continuous(breaks = c(0, 50, 100)) +
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1),
                     labels = c("0", "0.25", "0.5", "0.75", "1")) +
  facet_grid(model ~ discrepancy, labeller = labeller(
    model = model_labels,
    discrepancy = discrepancy_labels
  )) +
  labs(x = "p value", y = "count") +
  theme_minimal()
ggsave("figures/null_pvalues.png", dpi = 400, plot,
       width = 7, height = 5)
