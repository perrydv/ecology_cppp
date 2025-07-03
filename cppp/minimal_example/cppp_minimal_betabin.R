library(nimble)
library(tidyverse)
library(patchwork)

# load CPPP functions
source("cppp/sally_code/calculateCPPP_original.R")

# load functions to simulate data
source("simulate_data.R")

# load minimal model
source("cppp/minimal_example/cppp_minimal_model.R")

# load discrepancy functions
source("cppp/minimal_example/cppp_discrepancy_functions.R")

# simulation constants
# data size
nSites <- 50
nVisits <- 6
nDatasets <- 50

# parameter values
p <- 0.3
rho <- c(0.01, 0.1, 0.25, 0.5)
psi <- 0.6

# MCMC 
niter <- 5000
nburnin <- 1000
thin <- 5
nCalibrationReplicates <- 100

# discrepancy function list
discrepancyFunctions <- list(ratioDiscFunction, tukeyDiscFunction, 
                             devianceDiscFunction)

# discrepancy function arguments
args <- list(nVisits = nVisits, dataNames = "y",
             latent_occ = "z", prob_detection = "p")
discrepancyFunctionsArgs <- list(args, args, args)

############
# get cppp #
############

# create empty list of dataframes to store ppp and cppp
ppp_out <- list()
cppp_out <- list()
for (i in 1:2) {
  ppp_out[[i]] <- array(NA, dim = c(length(rho), nDatasets, 
                                    length(discrepancyFunctions)))
  cppp_out[[i]] <- array(NA, dim = c(length(rho), nDatasets, 
                                     length(discrepancyFunctions)))
}
names(ppp_out) <- c("true", "false")
names(cppp_out) <- c("true", "false")

# create model instance and compile
constants <- list(nSites = nSites, nVisits = nVisits)
model_uncompiled <- nimbleModel(
  model_minimal, constants = constants,
  data = list(y = simulate_betabinomial(params = list(psi = psi, p = p),
                                        nSites, nVisits, rho[1])))

model <- compileNimble(model_uncompiled)

# configure MCMC
mcmc_conf <- configureMCMC(model_uncompiled, monitors = c("psi", "p", "z"))

# build MCMC
mcmc <- buildMCMC(mcmc_conf)

# compile mcmc
compiled_mcmc <- compileNimble(mcmc, project = model)

# simulate n datasets
for (n in 1:nDatasets) {
  
  # loop through rho
  for (i in 1:length(rho)) {
    
    # simulate data
    y <- simulate_betabinomial(params = list(psi = psi, p = p),
                               nSites, nVisits, rho[i])
    
    # add data
    model_uncompiled$y <- y
    model$y <- y
    
    # add inits
    inits <- function(y) list(psi = runif(1, 0, 1), p = runif(1, 0, 1), 
                              z = pmin(y, 1)) 
    
    model_uncompiled$setInits(inits(model_uncompiled$y))
    model$setInits(inits(model$y))
    
    # generate samples
    MCMCOutput <- runMCMC(compiled_mcmc, niter = niter, 
                          nburnin = nburnin, thin = thin)
    
    ###############
    # cppp inputs #
    ###############
    
    dataNames <- "y"
    
    condition_on_latent_states <- c(TRUE, FALSE)
    
    # calculate for conditioning vs. not conditioning on latent state
    for (j in 1:length(condition_on_latent_states)) {
      
      if (condition_on_latent_states[j]) {
        # if conditioning on latent state
        paramNames <- c("p", "psi", "z")
        samples <- MCMCOutput
      } else {
        # if *not* conditioning on latent state
        paramNames <- c("p", "psi")
        samples <- MCMCOutput[, paramNames]
      }
      
      # run calibration
      out_cal <- runCalibration(
        model = model_uncompiled, dataNames = dataNames, paramNames = paramNames, 
        origMCMCSamples = samples, discrepancyFunctions = discrepancyFunctions, 
        discrepancyFunctionsArgs = discrepancyFunctionsArgs,
        nCalibrationReplicates = nCalibrationReplicates,
        returnSamples = FALSE, returnDiscrepancies = FALSE) 
      
      # add ppp
      ppp_out[[j]][i, n, ] <- out_cal$obsPPP
      
      # calculate cppp
      for (k in 1:length(out_cal$obsPPP)) {
        cppp_out[[j]][i, n, k] <- mean(out_cal$repPPP[k, ] <= out_cal$obsPPP[k])
      }
    }
  }
}


data_ratio <- data.frame(
  rho = c(rep(rho[1], 48), rep(rho[2], 48), rep(rho[3], 48), rep(rho[4], 48),
          rep(rho[1], 48), rep(rho[2], 48), rep(rho[3], 48), rep(rho[4], 48)), 
  ratio = c(cppp_out$true[1, 1:24, 1], ppp_out$true[1, 1:24, 1],
            cppp_out$true[2, 1:24, 1], ppp_out$true[2, 1:24, 1],
            cppp_out$true[3, 1:24, 1], ppp_out$true[3, 1:24, 1],
            cppp_out$true[4, 1:24, 1], ppp_out$true[4, 1:24, 1],
            cppp_out$false[1, 1:24, 1], ppp_out$false[1, 1:24, 1],
            cppp_out$false[2, 1:24, 1], ppp_out$false[2, 1:24, 1],
            cppp_out$false[3, 1:24, 1], ppp_out$false[3, 1:24, 1],
            cppp_out$false[4, 1:24, 1], ppp_out$false[4, 1:24, 1]),
  condition = c(rep(TRUE, 48 * 4), rep(FALSE, 48 * 4)),
  pvalue = c(rep("cppp", 24), rep("ppp", 24), rep("cppp", 24), rep("ppp", 24), 
             rep("cppp", 24), rep("ppp", 24), rep("cppp", 24), rep("ppp", 24),
             rep("cppp", 24), rep("ppp", 24), rep("cppp", 24), rep("ppp", 24), 
             rep("cppp", 24), rep("ppp", 24), rep("cppp", 24), rep("ppp", 24))
)

data_tukey <- data.frame(
  rho = c(rep(rho[1], 48), rep(rho[2], 48), rep(rho[3], 48), rep(rho[4], 48),
          rep(rho[1], 48), rep(rho[2], 48), rep(rho[3], 48), rep(rho[4], 48)), 
  tukey = c(cppp_out$true[1, 1:24, 2], ppp_out$true[1, 1:24, 2],
            cppp_out$true[2, 1:24, 2], ppp_out$true[2, 1:24, 2],
            cppp_out$true[3, 1:24, 2], ppp_out$true[3, 1:24, 2],
            cppp_out$true[4, 1:24, 2], ppp_out$true[4, 1:24, 2],
            cppp_out$false[1, 1:24, 2], ppp_out$false[1, 1:24, 2],
            cppp_out$false[2, 1:24, 2], ppp_out$false[2, 1:24, 2],
            cppp_out$false[3, 1:24, 2], ppp_out$false[3, 1:24, 2],
            cppp_out$false[4, 1:24, 2], ppp_out$false[4, 1:24, 2]),
  condition = c(rep(TRUE, 48 * 4), rep(FALSE, 48 * 4)),
  pvalue = c(rep("cppp", 24), rep("ppp", 24), rep("cppp", 24), rep("ppp", 24), 
             rep("cppp", 24), rep("ppp", 24), rep("cppp", 24), rep("ppp", 24),
             rep("cppp", 24), rep("ppp", 24), rep("cppp", 24), rep("ppp", 24), 
             rep("cppp", 24), rep("ppp", 24), rep("cppp", 24), rep("ppp", 24))
)

data_deviance <- data.frame(
  rho = c(rep(rho[1], 48), rep(rho[2], 48), rep(rho[3], 48), rep(rho[4], 48),
          rep(rho[1], 48), rep(rho[2], 48), rep(rho[3], 48), rep(rho[4], 48)), 
  deviance = c(cppp_out$true[1, 1:24, 3], ppp_out$true[1, 1:24, 3],
               cppp_out$true[2, 1:24, 3], ppp_out$true[2, 1:24, 3],
               cppp_out$true[3, 1:24, 3], ppp_out$true[3, 1:24, 3],
               cppp_out$true[4, 1:24, 3], ppp_out$true[4, 1:24, 3],
               cppp_out$false[1, 1:24, 3], ppp_out$false[1, 1:24, 3],
               cppp_out$false[2, 1:24, 3], ppp_out$false[2, 1:24, 3],
               cppp_out$false[3, 1:24, 3], ppp_out$false[3, 1:24, 3],
               cppp_out$false[4, 1:24, 3], ppp_out$false[4, 1:24, 3]),
  condition = c(rep(TRUE, 48 * 4), rep(FALSE, 48 * 4)),
  pvalue = c(rep("cppp", 24), rep("ppp", 24), rep("cppp", 24), rep("ppp", 24), 
             rep("cppp", 24), rep("ppp", 24), rep("cppp", 24), rep("ppp", 24),
             rep("cppp", 24), rep("ppp", 24), rep("cppp", 24), rep("ppp", 24), 
             rep("cppp", 24), rep("ppp", 24), rep("cppp", 24), rep("ppp", 24))
)

plot_true_ratio <- ggplot(data = data_ratio[data_ratio$condition == TRUE, ], 
                          aes(x = as.factor(rho), y = ratio, 
                              color = as.factor(pvalue))) +
  geom_jitter(position = position_jitterdodge(dodge.width = 0.8)) +
  labs(x = "rho", y = "p-value", color = "") +
  ggtitle("ratio") +
  theme_minimal()

plot_true_tukey <- ggplot(data = data_tukey[data_tukey$condition == TRUE, ], 
                          aes(x = as.factor(rho), y = tukey, 
                              color = as.factor(pvalue))) +
  geom_jitter(position = position_jitterdodge(dodge.width = 0.8)) +
  labs(x = "rho", y = "p-value", color = "") +
  ggtitle("freeman-tukey") +
  theme_minimal()

plot_true_deviance <- ggplot(data = data_deviance[data_deviance$condition == TRUE, ], 
                             aes(x = as.factor(rho), y = deviance, 
                                 color = as.factor(pvalue))) +
  geom_jitter(position = position_jitterdodge(dodge.width = 0.8)) +
  labs(x = "rho", y = "p-value", color = "") +
  ggtitle("deviance") +
  theme_minimal()

plot_true <- plot_true_ratio + plot_true_tukey + 
  plot_true_deviance + plot_layout(ncol = 1)
ggsave("cppp/minimal_example/figures/betabin_true.png", plot_true, 
       dpi = 400, height = 7, width = 4)


plot_false_ratio <- ggplot(data = data_ratio[data_ratio$condition == FALSE, ], 
                          aes(x = as.factor(rho), y = ratio, 
                              color = as.factor(pvalue))) +
  geom_jitter(position = position_jitterdodge(dodge.width = 0.8)) +
  labs(x = "rho", y = "p-value", color = "") +
  ggtitle("ratio") +
  theme_minimal()

plot_false_tukey <- ggplot(data = data_tukey[data_tukey$condition == FALSE, ], 
                          aes(x = as.factor(rho), y = tukey, 
                              color = as.factor(pvalue))) +
  geom_jitter(position = position_jitterdodge(dodge.width = 0.8)) +
  labs(x = "rho", y = "p-value", color = "") +
  ggtitle("freeman-tukey") +
  theme_minimal()

plot_false_deviance <- ggplot(data = data_deviance[data_deviance$condition == FALSE, ], 
                             aes(x = as.factor(rho), y = deviance, 
                                 color = as.factor(pvalue))) +
  geom_jitter(position = position_jitterdodge(dodge.width = 0.8)) +
  labs(x = "rho", y = "p-value", color = "") +
  ggtitle("deviance") +
  theme_minimal()

plot_false <- plot_false_ratio + plot_false_tukey + 
  plot_false_deviance + plot_layout(ncol = 1)
ggsave("cppp/minimal_example/figures/betabin_false.png", plot_false, 
       dpi = 400, height = 7, width = 4)

