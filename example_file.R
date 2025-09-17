# ---------------------------------------------------------------------- #
# Example of how to run some of the simulation cases from the paper.
# ---------------------------------------------------------------------- #
 
 ## Case 1: Constant treatment effect (Under the Null).
    ## Case 1a: No covariate effect (Adjusted = 0).
    ## Case 1b: Covariate effect (Adjusted = 1/2: Informative/Uninformative).
  ## Case 2: Constant treatment effect (Under the Alternative).
  ## Case 3: Time-varying treatment effect.
    ## Case 3a: Under the Null (i.e., symmetric across arms).
    ## Case 3b: Under the Alternative.

# Working directory.
setwd('~/Documents/GitHub/AUMCF-simulation')
rm(list = ls())

# Packages.
library(optparse)
library(MCC)
library(parallel)
library(dplyr)

# Functions for comparison methods.
source("comparison_methods.R")

# Functions for data generation.
source("data_generation.R")

#Number of simulation replicates. 
my_reps <- 100

# ------------------------- #
# Case 1a: Null, Unadjusted #
# ------------------------- #
params <- list(
  n = 200,  # Vary to be 100, 200, or 400.
  censor = 0.2,
  time = 4, # Vary to be 1 or 4.
  frailtyVar = 0, # Vary to be 0 or 3.
  BaseDeath0 = 0.2,
  BaseDeath1 = 0.2,
  BaseEvent0 = 1.0,
  BaseEvent1 = 1.0,
  adjusted = 0,
  symmetric = FALSE,
  TV_effect = 0,
  experiment = 1,
  reps = my_reps,
  out = "Test/"
)

source('run_simulation.R')


# -------------------------------------- #
# Case 1b: Null, Adjusted - Informative  #
# -------------------------------------- #
params <- list(
  n = 200, # Vary to be 100, 200, or 400.
  censor = 0.2,
  time = 4, # Vary to be 1 or 4.
  frailtyVar = 0, # Vary to be 0 or 3.
  BaseDeath0 = 0.2,
  BaseDeath1 = 0.2,
  BaseEvent0 = 1.0,
  BaseEvent1 = 1.0,
  adjusted = 1,
  symmetric = FALSE,
  TV_effect = 0,
  experiment = 1,
  reps = my_reps,
  out = "Test/"
)

source('run_simulation.R')


# ---------------------------------------- #
# Case 1b: Null, Adjusted - Uninformative  #
# ---------------------------------------- #
params <- list(
  n = 200, # Vary to be 100, 200, or 400.
  censor = 0.2,
  time = 4, # Vary to be 1 or 4.
  frailtyVar = 0, # Vary to be 0 or 3.
  BaseDeath0 = 0.2,
  BaseDeath1 = 0.2,
  BaseEvent0 = 1.0,
  BaseEvent1 = 1.0,
  adjusted = 2,
  symmetric = FALSE,
  TV_effect = 0,
  experiment = 1,
  reps = my_reps,
  out = "Test/"
)

source('run_simulation.R')


# ---------------------------- #
# Case 2: Non-Null, Unadjusted #
# ---------------------------- #
params <- list(
  n = 200,  # Vary to be 100, 200, or 400.
  censor = 0.2,
  time = 4, # Vary to be 1 or 4.
  frailtyVar = 0, # Vary to be 0 or 3.
  BaseDeath0 = 0.2,
  BaseDeath1 = 0.2,
  BaseEvent0 = 1.4,
  BaseEvent1 = 1.0,
  adjusted = 0,
  symmetric = FALSE,
  TV_effect = 0,
  experiment = 2,
  reps = my_reps,
  out = "Test/"
)

source('run_simulation.R')


# ----------------------------------- #
# Case 3a: Symmetric (Null) TV effect #
# ----------------------------------- #
params <- list(
  n = 200, # Vary to be 100, 200, or 400.
  censor = 0.2,
  time = 4, # Vary to be 1, 4, or 6.
  BaseDeath0 = 0.2,
  BaseDeath1 = 0.2,
  TV_effect = log(0.5),
  symmetric = TRUE, 
  adjusted = 0,
  experiment = 3,
  reps = my_reps,
  out = "Test/"
)

source('run_simulation.R')


# ------------------------------------------- #
# Case 3b: Non-Symmetric (Non-Null) TV effect #
# ------------------------------------------- #
params <- list(
  n = 200, # Vary to be 100, 200, or 400.
  censor = 0.2,
  time = 4, # Vary to be 1, 4, or 6.
  BaseDeath0 = 0.2,
  BaseDeath1 = 0.2,
  TV_effect = log(0.5),
  symmetric = FALSE,
  adjusted = 0,
  experiment = 3,
  reps = my_reps,
  out = "Test/"
)

source('run_simulation.R')

