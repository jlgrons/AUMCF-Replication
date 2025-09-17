# -----------------------------------------------------------------------------
# Main script to run the simulation settings.
# -----------------------------------------------------------------------------

# Convenience function for data generation.

# Generates data under the scenarios in the paper, including:
  ## Case 1: Constant treatment effect (Under the Null).
    ## Case 1a: No covariate effect.
    ## Case 1b: Covariate effect.
  ## Case 2: Constant treatment effect (Under the Alternative).
  ## Case 3: Time-varying treatment effect.
    ## Case 3a: Under the Null.
    ## Case 3b: Under the Alternative.

SimulateData <- function(params, calc_truth = FALSE){

  if (!calc_truth){
    
    # Keep censoring and n rate if not calculating the true parameter value.
    censoring_rate <- params$censor
    
    n <- params$n
    
  }else{
    
    # Set censoring rate to 0 to calculate the true parameter value. 
    censoring_rate <- 0
    
    # Increase arm size.
    n <- 2000 
     }
  
  if(params$experiment != 3){   
    
    # Generate for Case 1 or 2: Constant treatment effect.
    data <- SimData(n = n,
                    censoring_rate = censoring_rate, 
                    base_death_rate = params$BaseDeath0,
                    base_event_rate_0 = params$BaseEvent0, 
                    base_event_rate_1 = params$BaseEvent1,
                    frailty_variance = params$frailtyVar,
                    tau = params$time,
                    adjust = params$adjusted)

    
    }else if(params$experiment == 3){
    
    # Generate for Case 3: Time-varying treatment effect.
    data <- SimData_TV(n = n,
                    lambda_cens = censoring_rate,
                    lambda_death = params$BaseDeath0,
                    tau = params$time, 
                    beta = params$TV_effect,
                    symmetric = params$symmetric)
    }

}
  
  
# Convenience function to loop across simulation replicates.
SimulationLoop <- function(i) {
  
  set.seed(i)
  
  # Generate data.
  data <- SimulateData(params)
  
  # Run all methods.
  results <- RunAllMethods(data, params$time, params$adjusted)
}


# -----------------------------------------------------------------------------
# Run the simulation and format the output.
# -----------------------------------------------------------------------------
t0 <- proc.time()
output <- lapply(1:params$reps, SimulationLoop) 
sim <- do.call(rbind, output)


# -----------------------------------------------------------------------------
# Compute the true parameter values.  
# -----------------------------------------------------------------------------
# Method names. 
method_names <-unique(sim$type)

if(params$experiment == 1 | (params$experiment == 3 & params$symmetric == TRUE) ){
  
  # Null Case.
  if(params$adjusted == 0){
    
    truth_values <- c(0, rep(1, length(method_names) - 1))
    truth <- setNames(truth_values, method_names)
    
  } else{
    
    truth_values <- c(0, rep(1, length(method_names) - 2), 0 )
    truth <- setNames(truth_values, method_names)
    
  }
  
  
}else{
  
  all_truth_values <- c()
  
  for(i in 1:1){ # Adjust if you need more - plot mean(truth) vs # replicates for all parameters.
    
    # Non-Null Case.
    big_data <- SimulateData(params, calc_truth = TRUE)
    
    # Return results.
    truth_results <- RunAllMethods(big_data, params$time)
    truth_values <- truth_results[, "value"]
    truth_i <- setNames(truth_values, method_names)
    
    all_truth_values <- rbind(all_truth_values, truth_i)
  }
  
  truth <- colMeans(all_truth_values)
  
}


# -----------------------------------------------------------------------------
# Summarize the results. 
# -----------------------------------------------------------------------------
# Add true value to the simulation.
sim_augmented <- sim %>%
  mutate(true_value = truth[type])
## NOTE: I also suggest you save all the simulation data here.

# Calculate summary stats.
summary_table <- sim_augmented %>%
  group_by(type) %>%
  summarise(
    pt_est = mean(value),
    bias = mean(value) - first(true_value),
    prob_reject_H0 = mean(p_value < 0.05, na.rm = TRUE),
    ase = mean(se),
    ese = sd(value),
    cov_p = mean(lower <= first(true_value) & first(true_value) <= upper),
    .groups = "drop"
  )

# Save summary stats + true value.
summary_table_tv <- summary_table %>%
  mutate(true_value = truth[type])
out <- data.frame(summary_table_tv)
print(out)

# Store simulation settings.
out$n <- params$n
out$time <- params$time
out$rep <- params$reps

# Output directory. 
out_stem <- paste0(params$out)


# Save results. 
if(params$experiment != 3 ){
  
sim_file <- paste0(out_stem, "sim",
                   "N", params$n, 
                   "_T", params$time,
                   "_l0", params$BaseEvent0,
                   "_l1", params$BaseEvent1,
                   "_F", params$frailtyVar,
                   "_c", params$censor,
                   ".rds")
}else{
  
  sim_file <- paste0(out_stem, "sim",
                     "N", params$n, 
                     "_T", params$time,
                     "_c", params$censor,
                     "_TV", round(params$TV_effect, 2) * -1,
                     ".rds")
  
}

saveRDS(object = sim_augmented, file = sim_file)

# -----------------------------------------------------------------------------
# Runtime.
# -----------------------------------------------------------------------------
t1 <- proc.time()
cat(t1-t0, "\n")



