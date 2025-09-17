# Purpose: Simulate data with a time-varying treatment effect.
# Updated: 2025-06-12

#' Simulate Observation Time and Status
#'
#' @param lambda_cens Rate for the time-to censoring.
#' @param lambda_death Rate for the time-to death.
#' @param tau Truncation time.
#' @return Data.frame.
.ObsTimeStatus <- function(
    lambda_cens,
    lambda_death,
    tau
) {
  
  # Time to censoring.
  if (lambda_cens > 0) {
    time_cens <- stats::rexp(n = 1, rate = lambda_cens)
    time_cens <- min(time_cens, tau)
  } else {
    time_cens <- tau
  }
  
  # Time to death.
  if (lambda_death > 0) {
    time_death <- stats::rexp(n = 1, rate = lambda_death)
  } else {
    time_death <- Inf
  }
  
  # Final status.
  obs_time <- ifelse(time_death <= time_cens, time_death, time_cens)
  status <- ifelse(obs_time == time_death, 2, 0)
  out <- data.frame(
    time = obs_time,
    status = status
  )
  return(out)
}


#' Find Max Lambda
#'
#' @param lambda_fn A *function* of time that returns the rate for the 
#'   recurrent event process.
#' @param tau Truncation time.
#' @param n_grid Number of grid points.
#' @return Scalar maximum event rate.
.LambdaMax <- function(lambda_fn, tau, n_grid = 100) {
  grid <- seq(from = 0, to = tau, length.out = n_grid)
  rates <- lambda_fn(grid)
  max_lambda <- max(rates)
  return(max_lambda)
}


#' Simulate for a Single Subject
#' 
#' Simulates data from a non-homogeneous Poisson process by thinning
#' a homogeneous Poisson process with rate = lambda_max.
#'
#' @param lambda_fn A *function* of time that returns the rate for the 
#'   recurrent event process.
#' @param lambda_max Maximum value of event rate.
#' @param lambda_cens Rate for the time-to censoring.
#' @param lambda_death Rate for the time-to death.
#' @param tau The truncation time.
#' @return Data.frame.
SimSubj <- function(
    lambda_fn,
    lambda_max,
    lambda_cens = 0.25,
    lambda_death = 0.25,
    tau = 4
) {
  
  # Time and status.
  df <- .ObsTimeStatus(
    lambda_cens = lambda_cens, 
    lambda_death = lambda_death,
    tau = tau
  )
  
  # Poisson process.
  event_times <- numeric(0)
  t <- 0
  while (t < df$time) {
    
    # Simulate gap time.
    gap <- stats::rexp(n = 1, rate = lambda_max)
    
    # Increment homogeneous process.
    t <- t + gap
    if (t > df$time) {break}
    
    # Decide whether to accept proposal.
    pi <- min(lambda_fn(t) / lambda_max, 1)
    accept <- stats::rbinom(n = 1, size = 1, prob = pi)
    if (accept == 1) {
      event_times <- c(event_times, t)
    }
    
  }
  
  # Output.
  if (length(event_times) > 0) {
    out <- data.frame(
      time = event_times,
      status = 1
    )
    out <- rbind(out, df)
  } else {
    out <- df
  }
  return(out)
}


#' Simulate arm
#' 
#' @param lambda_fn A *function* of time that returns the rate for the 
#'   recurrent event process.
#' @param n Sample size.
#' @param lambda_cens Rate for the time-to censoring.
#' @param lambda_death Rate for the time-to death.
#' @param lambda_max Maximum value of event rate.
#' @param tau The truncation time.
#' @return Data.frame.
SimArm <- function(
    lambda_fn,
    n,
    lambda_max = NULL,
    lambda_cens = 0.25,
    lambda_death = 0.25,
    tau = 4
) {
  
  lambda_fn <- Vectorize(lambda_fn)
  
  # Find lambda max.
  if (is.null(lambda_max)) {
    lambda_max <- .LambdaMax(lambda_fn, tau)
  }
  
  # Generate data.
  data <- lapply(seq_len(n), function(idx) {
    
    df <- SimSubj(
      lambda_fn = lambda_fn,
      lambda_max = lambda_max,
      lambda_cens = lambda_cens,
      lambda_death = lambda_death,
      tau = tau
    )
    df$idx <- idx
    return(df)
    
  })
  data <- do.call(rbind, data)
  return(data)
}


#' Simulate data
#' 
#' Set the rate functions for the treatment and control arms within this function.
#' 
#' @param n Sample size.
#' @param lambda_cens Rate for the time-to censoring.
#' @param lambda_death Rate for the time-to death.
#' @param tau The truncation time.
#' @return Data.frame.
SimData_TV <- function(
    n,
    lambda_cens = 0.25,
    lambda_death = 0.25,
    tau = 4.0, 
    beta = log(1), 
    symmetric = TRUE
) {
  
  # Base event rate.
  lambda_base <- 1.0
  
  # Break point, after which the rate for the treatment arm switches
  # from lambda_base to lambda_base * exp(beta).
  if (tau > 1){
    bp <- 1.0
  } else {
    bp <- 0.5
  }
  
  
  # Rate for treatment arm.
  rate_1 <- function(time) {
    if (time <= bp) {
      out <- lambda_base
    } else {
      out <- lambda_base * exp(beta)
    }
    return(out)
  }
  
  # Rate for control arm.
  if(symmetric){ # If under H0, then generate the same way for the control arm. 
    rate_0 <- function(time) {
      if (time <= bp) {
        out <- lambda_base
      } else {
        out <- lambda_base * exp(beta)
      }
      return(out)
    }
  }else{
    rate_0 <- function(time) {
      out <- lambda_base
      return(out)
    }
  }

  # Simulate data for treatment arm.
  data_1 <- SimArm(
    lambda_fn = rate_1,
    n = n,
    lambda_cens = lambda_cens,
    lambda_death = lambda_death,
    tau = tau
  )
  data_1$arm <- 1
  
  # Simulate data for control arm.
  data_0 <- SimArm(
    lambda_fn = rate_0,
    n = n,
    lambda_cens = lambda_cens,
    lambda_death = lambda_death,
    tau = tau
  )
  data_0$arm <- 0
  data_0$idx <- data_0$idx + n
  
  # Final data set.
  data <- rbind(data_1, data_0)
  data <- data[, c("idx", "arm", "time", "status")]
  return(data)
  
}


SimDataArm <- function(n, censoring_rate, base_death_rate, base_event_rate, 
                    frailty_variance, tau,
                    beta_death = 0, beta_event = 0, adjust = FALSE,
                    min_death_rate = 0.05, min_event_rate = 0.05) {
  
  # n: number of subjects
  # base_event_rate: baseline recurrent event rate
  # frailty_variance: frailty variance
  # tau: end of follow-up
  # base_death_rate: baseline death hazard
  # censoring_rate: baseline censoring hazard
  # beta_death: covariate effect on death
  # beta_event: covariates effect on event
  
  # Add adjustment covariates if required. 
  if(adjust){
    
    covar <- rnorm(n)
    
  } else {
    
    covar <- rep(0, n)
  }
  
  # Calculate subject-specific event rate.  
  base_event_rate_i <- base_event_rate * exp(beta_event * covar)
  death_i <- base_death_rate * exp(beta_death * covar) # Note: Typo fixed 6/18.
  
  # Apply floor to death and event rates. 
  base_event_rate_i <- pmax(base_event_rate_i, min_death_rate)
  death_i <- pmax(death_i, min_event_rate)
  
  # Simulate Gamma frailty (mean = 1, variance = frailty_variance).
  if(frailty_variance > 0){
    
    frailty <- rgamma(n, shape = 1/frailty_variance, rate = 1/frailty_variance)
    
  } else {
    
    frailty <- rep(1, n)
    
  }
  
  # Apply frailty. 
  base_event_rate_i <-  base_event_rate_i * frailty 
  death_i <- death_i * frailty 
  cens_i <- censoring_rate  
  
  id_vec <- integer()
  time_vec <- numeric()
  status_vec <- integer()
  
  id <- 1
  for (i in seq_len(n)) {
    # Simulate death and censoring times
    t_death <- rexp(1, rate = death_i[i])
    
    if(cens_i != 0){
      t_cens <- rexp(1, rate = cens_i)
    } else {
      t_cens <- tau + 1
    }
    
    t_obs <- min(tau, t_death, t_cens)
    
    # Determine event type at end of follow-up
    if (t_obs == t_death) {
      status_terminal <- 2  # death
    } else {
      status_terminal <- 0  # censored or administrative
    }
    
    # Simulate recurrent events in [0, t_obs]
    num_events <- rpois(1, base_event_rate_i[i] * t_obs)
    if (num_events > 0) {
      event_times <- sort(runif(num_events, min = 0, max = t_obs))
    } else {
      event_times <- numeric(0)
    }
    
    # Store recurrent events
    id_vec <- c(id_vec, rep(id, length(event_times)))
    time_vec <- c(time_vec, event_times)
    status_vec <- c(status_vec, rep(1, length(event_times)))  # 1 = recurrent
    
    
    # Store terminal event (death or censoring)
    id_vec <- c(id_vec, id)
    time_vec <- c(time_vec, t_obs)
    status_vec <- c(status_vec, status_terminal)
    
    id <- id + 1
  }
  
  if(adjust){
    
    covars <- setNames(covar, unique(id_vec))
    
    sim_data_0 <- data.frame(idx = id_vec,
                           time = time_vec,
                           status = status_vec)
    
    sim_data <- sim_data_0 %>%
      mutate(covar = covars[idx])
    
  }else{
    
    sim_data <- data.frame(idx = id_vec,
                           time = time_vec,
                           status = status_vec)
  }
  
 return(sim_data) 
 
}


SimData <- function(n, censoring_rate, base_death_rate,
                    base_event_rate_0, base_event_rate_1,
                    frailty_variance, tau, adjust,
                    base_death_rate_1 = NULL) {
  
  if(adjust == 1){
    
    # Parameters set for simulation in the paper. 
    beta_death_0 <- log(0.5)
    beta_death_1 <- log(0.5)
    beta_event_0 <- log(2)
    beta_event_1 <- log(2)
    
  } else if (adjust == 2){
    
    # Parameters set for simulation in the paper. 
    beta_death_0 <- 0
    beta_death_1 <- 0
    beta_event_0 <- 0
    beta_event_1 <- 0
    
  } else {
    
    beta_death_0 <- 0
    beta_death_1 <- 0
    beta_event_0 <- 0
    beta_event_1 <- 0
    
  }
  
  if(is.null(base_death_rate_1)){
    
    base_death_rate_1 = base_death_rate
    
  }
  
  data_arm_0 <- SimDataArm(n, censoring_rate, base_death_rate,
             base_event_rate_0, 
             frailty_variance, tau, 
             beta_death_0, beta_event_0, adjust)
  data_arm_0$arm <- 0
  
  data_arm_1 <- SimDataArm(n, censoring_rate, base_death_rate_1,
                           base_event_rate_1, 
                           frailty_variance, tau, 
                           beta_death_1, beta_event_1, adjust)
  data_arm_1$arm <- 1
  data_arm_1$idx <- data_arm_1$idx + n
 
  data <- rbind(data_arm_0, data_arm_1)
}







