#########################################
# function for running cppp simulations #
#########################################

run_cppp_simulations <- function(
    constants, # model constants (list)
    simulated_data, # simulated data (array dimensions: nDatasets, length(breakage_axis), nSites)
    model_uncompiled,
    MCMC_monitors, # params to monitor in MCMC
    breakage_axis, # vector
    data_name_list, # conditioned vs. not conditioned
    param_name_list, # conditioned vs. not conditioned
    param_indices_list, # conditioned vs. not conditioned
    discrepancyFunctions,
    discrepancyNames,
    discrepancyFunctionsArgs,
    coverage_params, # named list with true param values
    init_function, # function to gerenate initial values
    nDatasets,
    niter,
    nburnin,
    thin,
    nCalibrationReplicates
  ) {
  
  # create empty list of dataframes to store ppp and cppp
  ppp_out <- list()
  cppp_out <- list()
  coverage <- list()
  
  # compile model
  model <- compileNimble(model_uncompiled)
  
  # configure MCMC
  mcmc_conf <- configureMCMC(model_uncompiled, monitors = MCMC_monitors)
  
  # build MCMC
  mcmc <- buildMCMC(mcmc_conf)
  
  # compile mcmc
  compiled_mcmc <- compileNimble(mcmc, project = model, resetFunctions = TRUE)
  
  condition_on_latent_states <- c(TRUE, FALSE)
  
  # calculate for conditioning vs. not conditioning on latent state
  for (j in 1:length(condition_on_latent_states)) {
    
    if (condition_on_latent_states[j]) {
      # if conditioning on latent state
      dataNames <- data_name_list[[j]]
      paramNames <- param_name_list[[j]]
      paramIndices <- param_indices_list[[j]]
      simNodes <- unique(c(model$expandNodeNames(dataNames), 
                           model$getDependencies(paramNames, includeData = FALSE, 
                                                 self = FALSE)))
      
      # create empty arrays
      ppp_out[[j]] <- array(NA, dim = c(length(breakage_axis), nDatasets,
                                        length(discrepancyFunctions)),
                            dimnames = list(breakage_axis, 1:nDatasets, 
                                            discrepancyNames))
      cppp_out[[j]] <- array(NA, dim = c(length(breakage_axis), nDatasets,
                                         length(discrepancyFunctions)),
                             dimnames = list(breakage_axis, 1:nDatasets, 
                                             discrepancyNames))
      coverage[[j]] <- array(NA, dim = c(length(breakage_axis), nDatasets, 2),
                             dimnames = list(breakage_axis, 1:nDatasets, 
                                             names(coverage_params)))
      
      
    } else {
      # if *not* conditioning on latent state
      dataNames <- data_name_list[[j]]
      paramNames <- param_name_list[[j]]
      paramIndices <- param_indices_list[[j]]
      simNodes <- unique(c(model$expandNodeNames(dataNames), 
                           model$getDependencies(paramNames, includeData = FALSE, 
                                                 self = FALSE)))
      
      # create empty arrays
      ppp_out[[j]] <- array(NA, dim = c(length(breakage_axis), nDatasets,
                                        length(discrepancyFunctions)),
                            dimnames = list(breakage_axis, 1:nDatasets, 
                                            discrepancyNames))
      cppp_out[[j]] <- array(NA, dim = c(length(breakage_axis), nDatasets,
                                         length(discrepancyFunctions)),
                             dimnames = list(breakage_axis, 1:nDatasets, 
                                             discrepancyNames))
      coverage[[j]] <- array(NA, dim = c(length(breakage_axis), nDatasets, 
                                         length(coverage_params)),
                             dimnames = list(breakage_axis, 1:nDatasets, 
                                             names(coverage_params)))
      
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
      for (i in 1:length(breakage_axis)) {
        
        # add simulated data **note: assumes data is named y
        model_uncompiled$y <- simulated_data[n, i, ]
        model$y <- simulated_data[n, i, ]
        
        # add inits
        model_uncompiled$setInits(init_function(model_uncompiled$y))
        model$setInits(init_function(model$y))
        
        # generate samples
        MCMCOutput <- runMCMC(compiled_mcmc, niter = niter, 
                              nburnin = nburnin, thin = thin)
        
        # get coverage
        for (p in 1:length(coverage_params)) {
          bounds <- quantile(MCMCOutput[, names(coverage_params)[p]], 
                             c(0.025, 0.975))
          coverage[[j]][i, n, p] <- coverage_params[p] >= bounds[1] && 
            coverage_params[p] <= bounds[2]
        }     
        
        
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
  
  # convert from lists to arrays
  data_ppp <- rbind(as.data.frame.table(ppp_out[[1]]) %>% 
                      mutate(method = "ppp", condition = TRUE),
                    as.data.frame.table(ppp_out[[2]]) %>% 
                      mutate(method = "ppp", condition = FALSE)) %>% 
    setNames(c("breakage_axis", "sim", "discrepancy", 
               "pvalue", "method", "condition"))
  
  data_cppp <- rbind(as.data.frame.table(cppp_out[[1]]) %>% 
                       mutate(method = "cppp", condition = TRUE),
                     as.data.frame.table(cppp_out[[2]]) %>% 
                       mutate(method = "cppp", condition = FALSE)) %>% 
    setNames(c("breakage_axis", "sim", "discrepancy", 
               "pvalue", "method", "condition"))
  
  data_coverage <- rbind(as.data.frame.table(coverage[[1]]) %>% 
                          mutate(condition = TRUE),
                        as.data.frame.table(coverage[[2]]) %>% 
                          mutate(condition = FALSE)) %>% 
    setNames(c("breakage_axis", "sim", "par", "coverage", "condition")) %>% 
    pivot_wider(names_from = par, values_from = coverage)
  
  all_data <- bind_rows(data_ppp, data_cppp) %>% left_join(data_coverage)
  
  return(all_data)
}


##############
# cppp plots #
##############

# density plot
get_cppp_density_plot <- function(all_data, cdtn) {
  
  title <- ifelse(cdtn, "conditioned on latent state",
                  "not conditioned on latent state")
  
  plot <- all_data %>% 
    filter(condition == cdtn) %>% 
    ggplot(aes(x = pvalue, 
               fill = as.factor(method), group = method)) +
    geom_density(alpha = 0.7) +
    xlim(c(0, 1)) +
    facet_grid(breakage_axis ~ discrepancy, scales = "free") +
    scale_color_manual(values = c("black", NA)) +
    scale_x_continuous(breaks = c(0, 0.5, 1),
                       labels = c("0", "0.5", "1")) +
    labs(x = "p-value", y = "count", fill = "") +
    ggtitle(title) +
    geom_hline(yintercept = 0.05, linetype = 2) +
    theme_minimal(base_family = "Arial")
  
  return(plot)
}

# dot plot
get_cppp_dot_plot <- function(all_data, param, breakage_axis_name, cdtn) {
  
  title <- ifelse(cdtn, paste0(param, " coverage, conditioned on latent state"),
                  paste0(param, " coverage, conditioned on latent state"))
  
  plot <- all_data %>% 
    filter(condition == cdtn) %>% 
    ggplot(aes(x = as.factor(breakage_axis), y = pvalue, 
               shape = as.factor(method), color = !!sym(param), 
               group = method)) +
    geom_jitter(position = position_jitterdodge(dodge.width = 0.8)) +
    facet_grid(discrepancy ~ .) +
    scale_color_manual(values = c("black", NA)) +
    labs(x = breakage_axis_name, y = "p-value", 
         color = "coverage", shape = "") +
    geom_hline(yintercept = 0.05, linetype = 2) +
    ggtitle(title) +
    theme_minimal(base_family = "Arial")
  
  return(plot)
}

# power plot
get_cppp_power_plot <- function(all_data, param, breakage_axis_name) {
  
  # convert data to long
  all_data_long <- all_data %>% 
    mutate(pvalue_disc = cut(pvalue, c(0, 0.05, 1), 
                             include.lowest = T)) %>% 
    pivot_longer(c(p, psi, all_param)) %>% 
    group_by(breakage_axis, discrepancy, method, 
             condition, pvalue_disc, name, value) %>% 
    tally() %>% 
    ungroup() %>% 
    complete(breakage_axis, discrepancy, method, 
             condition, pvalue_disc, name, value, 
             fill = list(n = 0))
  
  plot <- all_data_long %>% 
    filter(name == param) %>% 
    ggplot()+
    geom_bar(aes(x = breakage_axis, y = n, 
                 fill = interaction(pvalue_disc, value)), 
             stat = "identity")+
    facet_grid(discrepancy ~ condition + method) +
    labs(x = breakage_axis_name) +
    ggtitle(paste(param, " coverage")) +
    scale_fill_brewer(palette = "Set1",
                      labels = c("True Positive", "False Negative", 
                                 "False Positive", "True Negative"), 
                      name = "Power") +
    theme_minimal(base_family = "Arial")
  
  return(plot)
}
