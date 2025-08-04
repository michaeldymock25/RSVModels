
rm(list = ls())

set.seed(6824767)

library(deSolve)

source("R/functions.R")
source("R/parameters.R")
source("R/models.R")

times <- seq(0, 400, 0.25)
N_sim <- 2
y0 <- initial_values(mod = "base", size_months = parms$size_months)
out <- mod_base(y0 = y0, times = times, parms = parms, N_sim = N_sim)
agg_out <- aggregate_output(mod = "base", out = out, times = times, age_years = parms$age_years)
