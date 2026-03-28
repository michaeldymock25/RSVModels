
rm(list = ls())

library(deSolve)
library(parallel)

source("R/functions.R")
source("R/parameters.R")
source("R/models.R")

set.seed(6824767)

max_time <- 200
N_sim <- 2
y0 <- initial_values(mod = "base", size_months = parms$size_months, N_sim = N_sim)
out <- mod_base(y0 = y0, max_time = max_time, parms = parms, N_sim = N_sim)
inc <- extract_incidence(mod = "base", out = out, times = (max_time - 4):max_time)
