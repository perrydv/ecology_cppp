library(tidyverse)

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

ggplot(data = pvalues_nonind_long) +
  geom_histogram(aes(x = pvalue, color))

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


plot_wlatent <- ggplot(pvalues, aes(x = p_value, fill = latent)) +
  geom_histogram(position = "identity", alpha = 0.5) +
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
