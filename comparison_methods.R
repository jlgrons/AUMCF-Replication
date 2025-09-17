##########################################################################
############## All functions will return a data frame with ##############
#   value    se   lower    upper   p_value   type
#   (pt est)       (95% ci      )            (name of the model)
##########################################################################

# These functions assume that the data has the following columns:
# idx: Subject index.
# time: Time of event, censoring, or death.
# status: Coded 0 for censoring, 1 for an event, 2 for death.
# arm: Coded 0 for control and 1 for treatment arm. 

library(dplyr)
library(survival)
library(reReg)
library(MASS)
library(frailtypack)
library(WR)

#####################################################
# Cox Proportional Hazards Model (First Event Only) #
#####################################################
coxPHmodel <- function(data, composite = TRUE){
  
  # Extract first event for all subjects.
  first_event_data <- data %>%
    group_by(idx) %>%
    arrange(time) %>%
    slice(1) %>%
    ungroup()
  
 
  if(composite){
    
    # Treat death or event as the first event + fit model.
    fit_coxph <- coxph(Surv(time, status != 0) ~ arm, 
                       data = first_event_data)
    model_type <- "coxph"
    
  } else {
    
    # Treat death as censoring + fit model.
    fit_coxph <- coxph(Surv(time, status == 1) ~ arm, 
                       data = first_event_data)
    model_type <- "coxph_nc"
  }
  
  
  # Summarize results. 
  s_coxp <- summary(fit_coxph)
  coxp_hr <- s_coxp$coefficients[1,2]
  coxp_se <- s_coxp$coefficients[1,3]
  coxp_ci_l <- s_coxp$conf.int[1,3]
  coxp_ci_u <- s_coxp$conf.int[1,4]
  coxp_p <- s_coxp$coefficients[1,5]
  
  result <- data.frame(value = coxp_hr,
                       se = coxp_se,
                       lower = coxp_ci_l,
                       upper = coxp_ci_u,
                       p_value = coxp_p,
                       type  = model_type)
  return(result)
}


##############
# LWYY Model #
##############
LWYYmodel <- function(data){
  
  # Get data in the proper format.
  data_reReg <- data %>%
    mutate(recur_status = ifelse(status == 1, 1, 0), # Indicator for an event.
           terminal_status = ifelse(status == 2, 1, 0)) # Indicator for death.

  # Fit model.
  fit_lwyy <- tryCatch(reReg(
    Recur(time, idx, recur_status, terminal_status) ~ arm,
    data = data_reReg,
    model = "cox.LWYY"
    ), error=function(e) NULL)
  
  # # Get data in start-stop format. 
  # data_sorted <-data %>%
  #   arrange(idx, time) %>%
  #   group_by(idx) %>%
  #   mutate(start = lag(time, default = 0),  # Start time is last event or 0.
  #          stop = time) %>%                 # Stop time is current time.
  #   ungroup()
  # 
  # # Add indicator for the recurrent event; treat death as censoring.  
  # data_sorted <- data_sorted %>%
  #   mutate(event_indicator = ifelse(status == 1, 1, 0))
  # 
  # # Fit model.
  # fit_lwyy <- coxph(
  #   Surv(start, stop, event_indicator) ~ arm + cluster(idx),
  #   data = data_sorted)
  
  # Summarize results. 
  if(!is.null(fit_lwyy)){
    
    s_lwyy <- summary(fit_lwyy)
    lwyy_coef <- s_lwyy$coefficients.rec[1,1]
    lwyy_se <- s_lwyy$coefficients.rec[1,2]
    lwyy_hr <- exp(lwyy_coef)
    lwyy_ci_l <- exp(lwyy_coef - 1.96 * lwyy_se)
    lwyy_ci_u <- exp(lwyy_coef + 1.96 * lwyy_se)
    lwyy_p <- s_lwyy$coefficients[1,4]
    
    result <- data.frame(value = lwyy_hr,
                         se = lwyy_se,
                         lower = lwyy_ci_l,
                         upper = lwyy_ci_u,
                         p_value = lwyy_p,
                         type  = "lwyy")
  }else{
    
    result <- data.frame(value = NA,
                         se = NA,
                         lower = NA,
                         upper = NA,
                         p_value = NA,
                         type  = "lwyy")
  }     
  
  return(result)
}


###########################
# Negative Binomial Model #
###########################
NegBmodel <- function(data){
  
  # Get data in the proper format. 
  nb_data <- data %>%
    group_by(idx, arm) %>%
    summarise(
      count = sum(status == 1), # Number of recurrent events.
      followup = max(time), # Total follow-up time.
      .groups = "drop"
    )
  
  # Fit model.
  fit_nb <- glm.nb(count ~ arm + offset(log(followup)), 
                   data = nb_data)
  
  # Summarize results. 
  s_nb <- summary(fit_nb)
  nb_coef <- s_nb$coefficients[2,1]
  nb_se <- s_nb$coefficients[2,2]
  nb_rr <- exp(nb_coef)
  nb_ci_l <- exp(nb_coef - 1.96 * nb_se)
  nb_ci_u <- exp(nb_coef + 1.96 * nb_se)
  nb_p <- s_nb$coefficients[2,4]
  
  result <- data.frame(value = nb_rr,
                       se = nb_se,
                       lower = nb_ci_l,
                       upper = nb_ci_u,
                       p_value = nb_p,
                       type  = "nb")
  return(result)
}


########################
# Joint Frailty Model #
#######################
FRYmodel <- function(data){
  
  # Get data in the proper format. 
  df_frailty <- data %>%
    filter(status %in% c(0, 1)) %>%
    arrange(idx, time) %>%
    group_by(idx) %>%
    mutate(start = lag(time, default = 0)) %>%
    ungroup()
  
  # Fit model.
  fit_frailty <- frailtyPenal(
    formula = Surv(start, time, status) ~ cluster(idx) + arm,
    data = df_frailty,
    recurrentAG = TRUE, 
    n.knots = 7,
    kappa = 1e-2
  )
  
  # Summarize results. 
  frailty_coef <- fit_frailty$coef["arm"]
  frailty_se <- sqrt(fit_frailty$varH)
  frailty_hr <- exp(frailty_coef)
  frailty_ci_l <- exp(frailty_coef - 1.96 * frailty_se)
  frailty_ci_u <- exp(frailty_coef + 1.96 * frailty_se)
  frailty_p <- fit_frailty$beta_p.value["arm"]
  
  result <- data.frame(value = frailty_hr,
                       se = frailty_se,
                       lower = frailty_ci_l,
                       upper = frailty_ci_u,
                       p_value = frailty_p,
                       type  = "frailty")
  return(result)
}


##############
# Win Ratio #
#############
WRstat <- function(data, reduced = TRUE){
  
  # Recode death + event for the purposes of the package.
  data <- data %>%
    mutate(status = case_when(
      status == 1 ~ 2, # 2 is the indicator of event.
      status == 2 ~ 1, # 1 is indicator of death.
      TRUE ~ status
    ))
  
  # Fit recurrent event win ratio statistics.
  # Details here: https://pmc.ncbi.nlm.nih.gov/articles/PMC9246892/.
  wr_rec_all <- WRrec(data[, "idx"],
                      data[, "time"],
                      data[, "status"],
                      data[, "arm"],
                      strata = NULL, 
                      naive = TRUE)
  
  # Summarize last-event-assisted win ratio (LWR) results. 
  result_LWR <- data.frame(value = exp(wr_rec_all$log.WR),
                           se = wr_rec_all$se * exp(wr_rec_all$log.WR),
                           lower = exp(wr_rec_all$log.WR - 1.96 * wr_rec_all$se),
                           upper = exp(wr_rec_all$log.WR + 1.96 * wr_rec_all$se),
                           p_value = wr_rec_all$pval,
                           type  = "wr_LWR")
  
  # Format data to calculate the standard win ratio (SWR).
  # LWR reduces to SWR with a non-recurrent, non-fatal event.  
  o <- order(data$idx, data$time)
  dat_o <- data[o,]
  dat_oc <-dat_o[!duplicated(dat_o[c("idx","status")]),] # Keep first + terminal event.
 
  # Calculate the standard win ratio.
  wr_rec_std <- WRrec(dat_oc[, "idx"],
                      dat_oc[, "time"],
                      dat_oc[, "status"],
                      dat_oc[, "arm"],
                      strata = NULL, 
                      naive = FALSE)
  
  # Summarize standard win ratio results. 
  result_STD <- data.frame(value = exp(wr_rec_std$log.WR),
                           se = wr_rec_std$se * exp(wr_rec_std$log.WR),
                           lower = exp(wr_rec_std$log.WR - 1.96 * wr_rec_std$se),
                           upper = exp(wr_rec_std$log.WR + 1.96 * wr_rec_std$se),
                           p_value = wr_rec_std$pval,
                           type  = "wr_STD")
  
  result <- rbind(result_LWR, result_STD)
  
  # If desired, add in other recurrent event win ratio statistics.
  if(!reduced){
    
    wr_z_FWR <- wr_rec_all$log.WR.FI/wr_rec_all$se.FI
    p_val_FWR <- 2 * (1 - pnorm(abs(wr_z_FWR)))
    result_FWR <- data.frame(value = exp(wr_rec_all$log.WR.FI),
                             se = wr_rec_all$se.FI * exp(wr_rec_all$log.WR.FI),
                             lower = exp(wr_rec_all$log.WR.FI - 1.96 * wr_rec_all$se.FI),
                             upper = exp(wr_rec_all$log.WR.FI + 1.96 * wr_rec_all$se.FI),
                             p_value =  p_val_FWR,
                             type  = "wr_FWR")
    
    wr_z_NWR <- wr_rec_all$log.WR.naive/wr_rec_all$se.naive
    p_val_NWR <- 2 * (1 - pnorm(abs(wr_z_NWR)))
    result_NWR <- data.frame(value = exp(wr_rec_all$log.WR.naive),
                             se = wr_rec_all$se.naive * exp(wr_rec_all$log.WR.naive),
                             lower = exp(wr_rec_all$log.WR.naive - 1.96 * wr_rec_all$se.naive),
                             upper = exp(wr_rec_all$log.WR.naive + 1.96 * wr_rec_all$se.naive),
                             p_value =  p_val_NWR,
                             type  = "wr_NWR")
    
    result <- rbind(result_LWR, result_FWR, result_NWR, result_STD)
    
    
  }
  
  return(result)
  
}


###############################
# Runs all comparison methods #
###############################
RunAllMethods <- function(data, tau, adjusted = FALSE){
  
  # Try to calculate AUMCF.
  boot <- try(
    MCC::CompareAUCs(data, tau = tau)
  )
  
  if (class(boot) != "try-error") {
    
    # Comparison Methods. 
    coxp <- coxPHmodel(data)
    coxp_nc <- coxPHmodel(data, composite = FALSE)
    lwyy <- LWYYmodel(data)
    nb <- NegBmodel(data)
    frailty <- FRYmodel(data)
    wr <- WRstat(data)
    
    # AUMCF (Difference).
    aumcf_diff <- data.frame(
      value = boot@CIs$observed[1],
      se = boot@CIs$se[1],
      lower = boot@CIs$lower[1],
      upper = boot@CIs$upper[1],
      p_value = boot@Pvals$p[1],
      type = "AUMCF_diff"
    )
    
    # AUMCF (Ratio).
    aumcf_ratio <- data.frame(
      value = boot@CIs$observed[2],
      se = boot@CIs$se[2],
      lower = boot@CIs$lower[2],
      upper = boot@CIs$upper[2],
      p_value = boot@Pvals$p[2],
      type = "AUMCF_ratio"
    )
    
    # Return results.
    results <- rbind(aumcf_diff, aumcf_ratio, coxp, coxp_nc, 
                     lwyy, nb, frailty, wr)
  }
  
  # Run adjusted analysis if necessary.
  if(adjusted == 1 | adjusted == 2){
    
    boot_adj <- try(
      MCC::CompareAUCs(data, tau = params$time, 
                       covars = data %>% dplyr::select(covar))
    )
    
    if (class(boot_adj) != "try-error") {
      aucmf_diff_adj <- data.frame(
        value = boot_adj@CIs$observed[1],
        se = boot_adj@CIs$se[1],
        lower = boot_adj@CIs$lower[1],
        upper = boot_adj@CIs$upper[1],
        p_value = boot_adj@Pvals$p[1],
        type = "adjAUMCF_diff"
      )
      results <- rbind(results, aucmf_diff_adj)
      
    }
  }
  
  
  rownames(results) <- 1:nrow(results)
  
  return(results)
  
}







