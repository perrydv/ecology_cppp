library(nimble)
library(tidyverse)
library(patchwork)

# load auxiliary functions
source("utils.R")


#################
# beta-binomial #
#################


# read in data
all_data <- readRDS("occupancy/saved_outputs/output_betabin.rds")

# print all plots
get_cppp_plot(plot_type = c("density", "dot_coverage", "dot_bias", "power"),
              all_data,
              param = c("p", "psi", "all_param"),
              breakage_axis_name = "rho", cdtn = c(TRUE, FALSE),
              print = TRUE)

# save all plots
get_cppp_plot(plot_type = c("density", "dot_coverage", "dot_bias", "power"), 
              all_data,
              param = c("p", "psi", "all_param"),
              breakage_axis_name = "rho", cdtn = c(TRUE, FALSE),
              print = FALSE, save = TRUE,
              filepath = "figures/occupancy/betabin")

# print one plot
get_cppp_plot(plot_type = "dot_coverage", all_data, param = "p",
              breakage_axis_name = "rho", cdtn = FALSE, print = TRUE)


##################
# site-level cov #
##################


# read in data
all_data <- readRDS("occupancy/saved_outputs/output_sitecov.rds")

# save all plots
get_cppp_plot(plot_type = c("density", "dot_coverage", "dot_bias", "power"), 
              all_data,
              param = c("p", "beta[1]", "beta[2]", "all_param"),
              breakage_axis_name = "beta2", cdtn = FALSE,
              print = FALSE, save = TRUE,
              filepath = "figures/occupancy/sitecov")
