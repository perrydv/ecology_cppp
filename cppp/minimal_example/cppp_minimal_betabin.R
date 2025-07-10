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
nDatasets <-100

# parameter values
p <- 0.3
rho <- c(0.01, 0.1, 0.25, 0.5)
psi <- 0.6

# MCMC 
niter <- 5000
nburnin <- 1000
thin <- 5
nCalibrationReplicates <- 100


############
# get cppp #
############

# create empty list of dataframes to store ppp and cppp
ppp_out <- list()
cppp_out <- list()
coverage <- list()


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
compiled_mcmc <- compileNimble(mcmc, project = model, resetFunctions = TRUE)

condition_on_latent_states <- c(TRUE, FALSE)

# calculate for conditioning vs. not conditioning on latent state
for (j in 1:length(condition_on_latent_states)) {
  
  if (condition_on_latent_states[j]) {
    # if conditioning on latent state
    dataNames <- "y"
    paramNames <- c("p", "psi", "z")
    simNodes <- unique(c(model$expandNodeNames(dataNames), 
                         model$getDependencies(paramNames, includeData = FALSE, 
                                               self = FALSE)))
    paramIndices <- 1:(nSites + 2)
    
    # discrepancy function list
    discrepancyFunctions <- list(ratioDiscFunction, devianceDiscFunction, 
                                 chisqDiscFunction, tukeyDiscFunction)
    
    discrepancyNames <- c("Ratio", "Deviance", "Chi-Square", "Freeman-Tukey")
    
    # discrepancy function arguments
    args <- list(nVisits = nVisits, dataNames = "y",
                 latent_occ = "z", prob_detection = "p", prob_occupancy = "psi")
    discrepancyFunctionsArgs <- list(args, args, args, args)
    
    ppp_out[[j]] <- array(NA, dim = c(length(rho), nDatasets,
                                      length(discrepancyFunctions)),
                          dimnames = list(rho, 1:nDatasets, discrepancyNames))
    cppp_out[[j]] <- array(NA, dim = c(length(rho), nDatasets,
                                       length(discrepancyFunctions)),
                           dimnames = list(rho, 1:nDatasets, discrepancyNames))
    coverage[[j]] <- array(NA, dim = c(length(rho), nDatasets, 2),
                           dimnames = list(rho, 1:nDatasets, c("psi", "p")))
    
    
  } else {
    # if *not* conditioning on latent state
    dataNames <- c("y", "z")
    paramNames <- c("p", "psi")
    simNodes <- unique(c(model$expandNodeNames(dataNames), 
                         model$getDependencies(paramNames, includeData = FALSE, 
                                               self = FALSE)))
    paramIndices <- 1:2
    
    # discrepancy function list
    discrepancyFunctions <- list(ratioDiscFunction_NoLatent, 
                                 devianceDiscFunction_NoLatent,
                                 chisqDiscFunction_NoLatent,
                                 tukeyDiscFunction_NoLatent)
    
    discrepancyNames <- c("Ratio", "Deviance", "Chi-Square", "Freeman-Tukey")
    # discrepancy function arguments
    args <- list(nVisits = nVisits, dataNames = "y",
                 latent_occ = "z", prob_detection = "p", prob_occupancy = "psi")
    discrepancyFunctionsArgs <- list(args, args, args,args)
    
    ppp_out[[j]] <- array(NA, dim = c(length(rho), nDatasets,
                                      length(discrepancyFunctions)),
                          dimnames = list(rho, 1:nDatasets, discrepancyNames))
    cppp_out[[j]] <- array(NA, dim = c(length(rho), nDatasets,
                                       length(discrepancyFunctions)),
                           dimnames = list(rho, 1:nDatasets, discrepancyNames))
    coverage[[j]] <- array(NA, dim = c(length(rho), nDatasets, 2),
                           dimnames = list(rho, 1:nDatasets, c("psi", "p")))
    
  }
  
  
  # create compiled objects for calculating ppp and cppp
  modelCalcDisc <- calcDiscrepancies(model = model_uncompiled,
                                     dataNames = dataNames,
                                     paramNames = paramNames,
                                     paramIndices = paramIndices,
                                     simNodes = simNodes,
                                     discrepancyFunctions = discrepancyFunctions,
                                     discrepancyFunctionsArgs = discrepancyFunctionsArgs)
  
  cModelCalcDisc <- compileNimble(modelCalcDisc, project = model_uncompiled,
                                  resetFunctions = TRUE)
  
  setAndSimPP <- setAndSimNodes(model = model_uncompiled, 
                                nodes = paramNames, 
                                simNodes = simNodes)
  
  cSetAndSimPP <- compileNimble(setAndSimPP, project = model_uncompiled,
                                resetFunctions = TRUE)
  
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
      
      
      
      mcmc_psi_q <- quantile(MCMCOutput[, "psi"], c(0.025, 0.975))
      mcmc_p_q <- quantile(MCMCOutput[, "p"], c(0.025, 0.975))
      coverage[[j]][i, n, 1] <- psi >= mcmc_psi_q[1] && psi <= mcmc_psi_q[2]
      coverage[[j]][i, n, 2] <- p >= mcmc_p_q[1] && psi <= mcmc_p_q[2]      

      
      ###############
      # cppp inputs #
      ###############
      
      if (condition_on_latent_states[j]) {
        # if conditioning on latent state
        samples <- MCMCOutput
      } else {
        # if *not* conditioning on latent state
        samples <- MCMCOutput[, paramNames]
      }
      
      # run calibration
      out_cal <- runCalibration_sim(
        model = model_uncompiled, dataNames = dataNames, 
        paramNames = paramNames, 
        origMCMCSamples = samples, cModelCalcDisc = cModelCalcDisc, 
        cMcmc = compiled_mcmc, cSetAndSimPP = cSetAndSimPP,
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

data_ppp <- rbind(as.data.frame.table(ppp_out[[1]]) %>% 
                    mutate(method = "ppp", condition = TRUE),
                          as.data.frame.table(ppp_out[[2]]) %>% 
                    mutate(method = "ppp", condition = FALSE)) %>% 
  setNames(c("rho", "sim", "discrepancy", "pvalue", "method", "condition"))

data_cppp <- rbind(as.data.frame.table(cppp_out[[1]]) %>% 
                     mutate(method = "cppp", condition = TRUE),
                  as.data.frame.table(cppp_out[[2]]) %>% 
                    mutate(method = "cppp", condition = FALSE)) %>% 
  setNames(c("rho", "sim", "discrepancy", "pvalue", "method", "condition"))

data_coverage = rbind(as.data.frame.table(coverage[[1]]) %>% 
                        mutate(condition = TRUE),
      as.data.frame.table(coverage[[2]]) %>% mutate(condition = FALSE)) %>% 
  setNames(c("rho", "sim", "par", "coverage", "condition")) %>% 
  pivot_wider(names_from = par, values_from = coverage)


alldata <- bind_rows(data_ppp, data_cppp) %>% left_join(data_coverage)

density_condition <- alldata %>% 
  filter(condition == T) %>% 
  ggplot(aes(x = pvalue, 
             fill = as.factor(method), group = method)) +
  geom_density(alpha = 0.7) +
  xlim(c(0, 1)) +
  facet_grid(rho ~ discrepancy, scales = "free") +
  scale_color_manual(values = c("black", NA)) +
  scale_x_continuous(breaks = c(0, 0.5, 1),
                     labels = c("0", "0.5", "1")) +
  labs(x = "p-value", y = "count", fill = "") +
  ggtitle("conditioned on latent state") +
  geom_hline(yintercept = 0.05, linetype = 2) +
  theme_minimal(base_family = "Arial")
ggsave("cppp/minimal_example/figures/density_true.png",
       density_condition, dpi = 400, height = 6, width = 6)

density_nocondition <- alldata %>% 
  filter(condition == F) %>% 
  ggplot(aes(x = pvalue, 
             fill = as.factor(method), group = method)) +
  geom_density(alpha = 0.7) +
  xlim(c(0, 1)) +
  facet_grid(rho ~ discrepancy, scales = "free")+
  scale_color_manual(values = c("black", NA))+
  scale_x_continuous(breaks = c(0, 0.5, 1),
                     labels = c("0", "0.5", "1")) +
  labs(x = "p-value", y = "count", fill = "") +
  ggtitle("not conditioned on latent state") +
  labs(x = "rho", y = "p-value", color = "") +
  geom_hline(yintercept = 0.05, linetype = 2) +
  theme_minimal(base_family = "Arial")
ggsave("cppp/minimal_example/figures/density_false.png",
       density_nocondition, dpi = 400, height = 6, width = 6)

coverpsi_condition <- alldata %>% 
  filter(condition == T) %>% 
  ggplot(aes(x = as.factor(rho), y = pvalue, 
             shape = as.factor(method), color = psi, group = method)) +
  geom_jitter(position = position_jitterdodge(dodge.width = 0.8)) +
  facet_grid(discrepancy ~ .) +
  scale_color_manual(values = c("black", NA)) +
  labs(x = "rho", y = "p-value", color = "coverage", shape = "") +
  geom_hline(yintercept = 0.05, linetype = 2) +
  ggtitle("psi coverage, conditioned on latent state") +
  theme_minimal(base_family = "Arial")
ggsave("cppp/minimal_example/figures/cover_psi_true.png",
       coverpsi_condition, dpi = 400, height = 6, width = 6)


coverpsi_nocondition <- alldata %>% 
  filter(condition == F) %>% 
  ggplot(aes(x = as.factor(rho), y = pvalue, 
             shape = as.factor(method), color = psi, group = method)) +
  geom_jitter(position = position_jitterdodge(dodge.width = 0.8)) +
  facet_grid(discrepancy ~ .) +
  scale_color_manual(values = c("black", NA)) +
  labs(x = "rho", y = "p-value", color = "coverage", shape = "") +
  geom_hline(yintercept = 0.05, linetype = 2) +
  ggtitle("psi coverage, not conditioned on latent state") +
  theme_minimal(base_family = "Arial")
ggsave("cppp/minimal_example/figures/cover_psi_false.png",
       coverpsi_nocondition, dpi = 400, height = 6, width = 6)

alldata_long <- alldata %>% 
  mutate(pvalue_disc = cut(pvalue, c(0, 0.05, 1), include.lowest = T)) %>% 
  pivot_longer(c(p, psi)) %>% 
  group_by(rho, discrepancy, method, condition, pvalue_disc, name, value) %>% 
  tally() %>% 
  ungroup() %>% 
  complete(rho, discrepancy, method, condition, pvalue_disc, name, value, 
           fill = list(n = 0))

power_plot <- alldata_long %>% 
  filter(name == "psi") %>% 
  ggplot()+
  geom_bar(aes(x = rho, y = n, fill = interaction(pvalue_disc, value)), 
           stat = "identity")+
  facet_grid(discrepancy ~ condition + method)+
  scale_fill_brewer(palette = "Set1",
                    labels = c("True Positive", "False Negative", 
                               "False Positive", "True Negative"), 
                    name = "Power") +
  theme_minimal(base_family = "Arial")
ggsave("cppp/minimal_example/figures/power.png",
       power_plot, dpi = 400, height = 6, width = 6)         
