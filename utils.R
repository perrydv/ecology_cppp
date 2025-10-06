#########################################
# function for running cppp simulations #
#########################################

run_cppp_simulations <- function(
  constants, # model constants (list)
  simulated_data, # simulated data (list)
  sim_data_names,
  model_uncompiled,
  mcmc_monitors, # params to monitor in MCMC
  breakage_axis, # vector
  data_name_list, # conditioned vs. not conditioned
  param_name_list, # conditioned vs. not conditioned
  # param_indices_list, # conditioned vs. not conditioned
  discrepancyFunctions,
  discrepancyNames,
  discrepancyFunctionsArgs,
  coverage_params, # named list with true param values
  init_function, # function to gerenate initial values
  init_args,
  nDatasets,
  niter,
  nburnin,
  thin,
  nchain,
  nCalibrationReplicates,
  condition_on_latent_states,
  MCMCcontrol = list(niter = 500,
                     thin = 1,
                     nburnin = 0)
) {

  # create empty list of dataframes to store ppp and cppp
  ppp_out <- list()
  cppp_out <- list()
  coverage <- list()
  bias <- list()

  # compile model
  model <- compileNimble(model_uncompiled)

  # configure MCMC
  mcmc_conf <- configureMCMC(model_uncompiled, monitors = mcmc_monitors)

  # build MCMC
  mcmc <- buildMCMC(mcmc_conf)

  # compile mcmc
  compiled_mcmc <- compileNimble(mcmc, project = model, resetFunctions = TRUE)

  # calculate for conditioning vs. not conditioning on latent state
  for (j in seq_along(condition_on_latent_states)) {
    
    print(paste0("condition: ", j))
    
    if (condition_on_latent_states[j]) {
      # if conditioning on latent state
      dataNames <- data_name_list[[j]]
      paramNames <- model$expandNodeNames(param_name_list[[j]])
      # paramIndices <- param_indices_list[[j]]
      simNodes <- unique(c(model$expandNodeNames(dataNames),
                           model$getDependencies(paramNames,
                                                 includeData = FALSE,
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
      bias[[j]] <- array(NA, dim = c(length(breakage_axis), nDatasets, 
                                     length(coverage_params)),
                         dimnames = list(breakage_axis, 1:nDatasets, 
                                         names(coverage_params)))
      
      
    } else {
      # if *not* conditioning on latent state
      dataNames <- data_name_list[[j]]
      paramNames <- model$expandNodeNames(param_name_list[[j]])
      # paramIndices <- param_indices_list[[j]]
      simNodes <- unique(c(model$expandNodeNames(dataNames), 
                           model$getDependencies(paramNames, 
                                                 includeData = FALSE, 
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
      bias[[j]] <- array(NA, dim = c(length(breakage_axis), nDatasets, 
                                     length(coverage_params)),
                         dimnames = list(breakage_axis, 1:nDatasets, 
                                         names(coverage_params)))
      
    }
    
    
    # create compiled objects for calculating ppp and cppp
    modelCalcDisc <- calcDiscrepancies(model = model_uncompiled, 
                                       dataNames = dataNames,
                                       paramNames = paramNames,
                                       #paramIndices = paramIndices,
                                       simNodes = simNodes,
                                       discrepancyFunctions = 
                                         discrepancyFunctions,
                                       discrepancyFunctionsArgs = 
                                         discrepancyFunctionsArgs)
    
    cModelCalcDisc <- compileNimble(modelCalcDisc, project = model_uncompiled,
                                    resetFunctions = TRUE)
    
    setAndSimPP <- setAndSimNodes(model = model_uncompiled, 
                                  nodes = paramNames, 
                                  simNodes = simNodes)
    
    cSetAndSimPP <- compileNimble(setAndSimPP, project = model_uncompiled,
                                  resetFunctions = TRUE)
    
    # simulate n datasets
    for (n in 1:nDatasets) {
      
      print(paste0("dataset: ", n))
      
      # loop through breakage axis
      for (i in seq_along(breakage_axis)) {
        
        print(paste0("breakage axis: ", i))
        
        simulated_data_sub <- list()
        for (a in seq_along(sim_data_names)) {
          dims <- length(dim(simulated_data[[a]]))
          if (dims == 3) {
            simulated_data_sub[[a]] <- simulated_data[[a]][n, i, ] 
          } else if (dims == 4) {
            simulated_data_sub[[a]] <- simulated_data[[a]][n, i, , ] 
          }
        }
        
        # add data and inits
        add_data_inits(model_uncompiled, model, simulated_data_sub, 
                       sim_data_names,
                       init_function, init_args)
        
        # generate samples
        mcmc_out <- runMCMC(compiled_mcmc, niter = niter, 
                            nburnin = nburnin, thin = thin,
                            nchain = nchain)
        MCMCOutput <- do.call(rbind, mcmc_out)
        
        
        # get coverage and bias
        for (p in seq_along(coverage_params)) {
          bounds <- quantile(MCMCOutput[, names(coverage_params)[p]], 
                             c(0.025, 0.975),
                             na.rm = TRUE)
          coverage[[j]][i, n, p] <- coverage_params[p] >= bounds[1] && 
            coverage_params[p] <= bounds[2]
          bias[[j]][i, n, p] <- abs(
            mean(MCMCOutput[, names(coverage_params)[p]],
                 na.rm = TRUE) - coverage_params[p]
          )
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
        out_cal <- runCalibration_sim(model = model_uncompiled,
                                      dataNames = dataNames, 
                                      paramNames = paramNames,
                                      origMCMCSamples = samples,
                                      cModelCalcDisc = cModelCalcDisc,
                                      cMcmc = compiled_mcmc,
                                      cSetAndSimPP = cSetAndSimPP,
                                      nCalibrationReplicates = 
                                        nCalibrationReplicates,
                                      MCMCcontrol,
                                      returnSamples = TRUE,
                                      returnDiscrepancies = TRUE) 
        
        print("cal done")
        
        # add ppp
        ppp_out[[j]][i, n, ] <- out_cal$obsPPP
        
        # calculate cppp
        if (length(out_cal$obsPPP) > 1) {
          for (k in seq_along(out_cal$obsPPP)) {
            cppp_out[[j]][i, n, k] <- mean(out_cal$repPPP[k, ] <= 
                                             out_cal$obsPPP[k])
          }
        } else {
          cppp_out[[j]][i, n, 1] <- mean(out_cal$repPPP <= 
                                           out_cal$obsPPP)
        }
      }
    }
  }
  
  if (all(c("TRUE", "FALSE") %in% condition_on_latent_states)) {
    
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
    
    data_bias <- rbind(as.data.frame.table(bias[[1]]) %>% 
                         mutate(condition = TRUE),
                       as.data.frame.table(bias[[2]]) %>% 
                         mutate(condition = FALSE)) %>% 
      setNames(c("breakage_axis", "sim", "par", "bias", "condition")) %>% 
      pivot_wider(names_from = par, values_from = bias,
                  names_glue = "{par}_bias")
    
    all_data <- bind_rows(data_ppp, data_cppp) %>% left_join(data_coverage) %>% 
      left_join(data_bias)
  } else {
    
    # convert from lists to arrays
    data_ppp <- as.data.frame.table(ppp_out[[1]]) %>% 
      mutate(method = "ppp", condition = condition_on_latent_states) %>% 
      setNames(c("breakage_axis", "sim", "discrepancy", 
                 "pvalue", "method", "condition"))
    
    data_cppp <- as.data.frame.table(cppp_out[[1]]) %>% 
      mutate(method = "cppp", condition = condition_on_latent_states) %>% 
      setNames(c("breakage_axis", "sim", "discrepancy", 
                 "pvalue", "method", "condition"))
    
    data_coverage <- as.data.frame.table(coverage[[1]]) %>% 
      mutate(condition = condition_on_latent_states) %>% 
      setNames(c("breakage_axis", "sim", "par", "coverage", "condition")) %>% 
      pivot_wider(names_from = par, values_from = coverage)
    
    data_bias <- as.data.frame.table(bias[[1]]) %>% 
      mutate(condition = condition_on_latent_states) %>% 
      setNames(c("breakage_axis", "sim", "par", "bias", "condition")) %>% 
      pivot_wider(names_from = par, values_from = bias,
                  names_glue = "{par}_bias")
    
    all_data <- bind_rows(data_ppp, data_cppp) %>% left_join(data_coverage) %>% 
      left_join(data_bias)
  }
  
  return(all_data)
}

add_data_inits <- function(model_uncompiled, model, simulated_data,
                           sim_data_names, init_function, init_args) {
  
  # add data
  for (i in seq_along(sim_data_names)) {
    
    values(model, sim_data_names[i]) <- as.vector(simulated_data[[i]])
    values(model_uncompiled, 
           sim_data_names[i]) <- as.vector(simulated_data[[i]])
    
  }
  
  # add initial values
  model$setInits(init_function(simulated_data, init_args))
  model_uncompiled$setInits(init_function(simulated_data, init_args))
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
    geom_vline(xintercept = 0.05, linetype = 2) +
    theme_minimal(base_family = "Arial")
  
  return(plot)
}

# dot plot - bias
get_cppp_dot_bias_plot <- function(all_data, param, breakage_axis_name, cdtn) {
  
  param_full <- paste0(param, "_bias")
  
  title <- ifelse(cdtn, paste0(param, " bias, conditioned on latent state"),
                  paste0(param, " bias, not conditioned on latent state"))
  
  plot <- all_data %>% 
    filter(condition == cdtn) %>% 
    ggplot(aes(x = as.factor(breakage_axis), y = pvalue, 
               shape = as.factor(method), color = !!sym(param_full), 
               group = method)) +
    geom_jitter(position = position_jitterdodge(dodge.width = 0.8)) +
    facet_grid(discrepancy ~ .) +
    #scale_color_manual(values = c("black", NA)) +
    labs(x = breakage_axis_name, y = "p-value", 
         color = "bias", shape = "") +
    geom_hline(yintercept = 0.05, linetype = 2) +
    ggtitle(title) +
    theme_minimal(base_family = "Arial")
  
  return(plot)
}

# dot plot - coverage
get_cppp_dot_coverage_plot <- function(all_data, param, 
                                       breakage_axis_name, cdtn) {
  
  title <- ifelse(cdtn, paste0(param, " coverage, conditioned on latent state"),
                  paste0(param, " coverage, not conditioned on latent state"))
  
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
                             include.lowest = TRUE)) %>% 
    pivot_longer(param) %>% 
    group_by(breakage_axis, discrepancy, method, 
             condition, pvalue_disc, name, value) %>% 
    tally() %>% 
    ungroup() %>% 
    complete(breakage_axis, discrepancy, method, 
             condition, pvalue_disc, name, value, 
             fill = list(n = 0))
  
  plot <- all_data_long %>% 
    filter(name == param) %>% 
    ggplot() +
    geom_bar(aes(x = breakage_axis, y = n, 
                 fill = interaction(pvalue_disc, value)), 
             stat = "identity") +
    facet_grid(discrepancy ~ condition + method) +
    labs(x = breakage_axis_name) +
    ggtitle(paste(param)) +
    scale_fill_brewer(palette = "Set1",
                      labels = c("True Positive", "False Negative", 
                                 "False Positive", "True Negative"), 
                      name = "Power") +
    theme_minimal(base_family = "Arial")
  
  return(plot)
}

get_cppp_plot <- function(plot_type, all_data, param = NULL, 
                          breakage_axis_name = NULL, cdtn = NULL, 
                          print = TRUE, save = FALSE,
                          filepath = NULL) {
  
  if ("density" %in% plot_type) {
    
    for (i in seq_along(cdtn)) {
      
      plot <- get_cppp_density_plot(all_data, cdtn[i])
      if (print) {
        print(plot)
      }
      if (save) {
        ggsave(paste0(filepath, "/density_", cdtn[i], ".png"),
               plot, dpi = 400, height = 6, width = 6)
      }
    }
    
  } 
  
  if ("dot_coverage" %in% plot_type) {
    
    for (i in seq_along(cdtn)) {
      
      for (j in seq_along(param)) {
        
        plot <- get_cppp_dot_coverage_plot(all_data, param[j], 
                                           breakage_axis_name, cdtn[i])
        if (print) {
          print(plot)
        }
        if (save) {
          ggsave(paste0(filepath, "/cover_", param[j], "_", cdtn[i], ".png"),
                 plot, dpi = 400, height = 6, width = 6)
        }
        
      }
    }
  } 
  
  if ("dot_bias" %in% plot_type) {
    
    param_bias <- param[param != "all_param"]
    
    for (i in seq_along(cdtn)) {
      
      for (j in seq_along(param_bias)) {
        
        plot <- get_cppp_dot_bias_plot(all_data, param_bias[j], 
                                       breakage_axis_name, cdtn[i])
        if (print) {
          print(plot)
        }
        if (save) {
          ggsave(paste0(filepath, "/bias_", param_bias[j], "_", 
                        cdtn[i], ".png"),
                 plot, dpi = 400, height = 6, width = 6)
        }
        
      }
    }
  } 
  
  if ("power" %in% plot_type) {
    
    for (i in seq_along(param)) {
      
      plot <- get_cppp_power_plot(all_data, param[i], breakage_axis_name)
      if (print) {
        print(plot)
      }
      if (save) {
        ggsave(paste0(filepath, "/power_", param[i], ".png"),
               plot, dpi = 400, height = 6, width = 7)
      }
      
    }
  }
}
