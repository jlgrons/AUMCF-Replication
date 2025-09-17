# AUMCF-Simulation

This repository contains the code and results for the paper *"Nonparametric estimation of the total treatment effect with multiple outcomes in the presence of terminal events"* by Jessica Gronsbell, Zachary R. McCaw, Isabelle-Emmanuella Nogues, Xiangshan Kong, Tianxi Cai, Lu Tian, and LJ Wei. You can find the preprint [here](https://arxiv.org/abs/2412.09304).

## Repository Structure

comparison_methods.R: Functions for comparison methods.

data_generation.R: Functions for data generation.

run_simulation.R: Main script to run a simulation.

example_file.R: Helper script that shows how to pass parameter values for various simulation settings from the paper.


## Required packages

```r
install.packages(c("dplyr", "frailtypack", "MCC", "optparse", "parallel", "survival", "reReg", "WR"))
```

## Simple Example

The code below shows to generate data and obtain results. 

```r
# Packages.
library(optparse)
library(MCC)
library(parallel)
library(dplyr)
library(survival)
library(reReg)
library(MASS)
library(frailtypack)
library(WR)

# Functions for comparison methods.
source("comparison_methods.R")

# Functions for data generation.
source("data_generation.R")

# Simulation parameters.
params <- list(
  n = 200, # Sample size.
  censor = 0.2, # Censoring rate.
  time = 4, # Observation window.
  frailtyVar = 0, # Frailty variance.
  # Baseline rate of terminal events in both treatment arms.
  BaseDeath0 = 0.2,
  # Baseline rate of recurrent events in each treatment arm.
  BaseEvent0 = 1.0,
  BaseEvent1 = 1.0,
  # Indicator for whether to do an adjusted analysis.
  adjusted = 0
)


# Generate data.
set.seed(92047)
data <- SimData(n = params$n,
                censoring_rate = params$censor, 
                base_death_rate = params$BaseDeath0,
                base_event_rate_0 = params$BaseEvent0, 
                base_event_rate_1 = params$BaseEvent1,
                frailty_variance = params$frailtyVar,
                tau = params$time,
                adjust = params$adjusted)

# Run AUMCF analysis + comparison methods.
 results <- RunAllMethods(data, params$time, params$adjusted)
 results
```

The code ouputs the AUMCF analysis (**AUMCF_diff**), the Cox proportional hazards model (**coxph**), the LWYY method (**LWYY**), negative binomial regression (**nb**), a frailty model (**frailty**), the last-event-assisted win ratio (**wr_LWR**), and the standard win ratio (**wr_STD**).

```r
>  results
      value         se      lower    upper   p_value       type
1 0.3569739 0.50752825 -0.6377632 1.351711 0.4818328 AUMCF_diff
2 1.0204374 0.10791614  0.8259010 1.260796 0.8512899      coxph
3 1.0464455 0.06828976  0.9153500 1.196316 0.5061771       lwyy
4 1.0473404 0.07076926  0.9116913 1.203172 0.5133772         nb
5 1.0023491 0.07075777  0.8725469 1.151461 0.9735471    frailty
6 1.0350103 0.11884749  0.8264216 1.296247 0.7644218     wr_LWR
7 1.1194387 0.13183433  0.8886975 1.410089 0.3380397     wr_STD
```
