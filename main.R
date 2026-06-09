
rm(list = ls())

library(deSolve)
library(parallel)

source("R/functions.R")
source("R/parameters.R")
source("R/models.R")

set.seed(6824767)

max_time <- 20
N_sim <- 2
y0_base <- initial_values(mod = "base", size_months = parms$size_months, N_sim = N_sim)
out_base <- mod_base(y0 = y0_base, max_time = max_time, parms = parms, N_sim = N_sim)
inc_base <- extract_incidence(mod = "base", out = out_base, times = (max_time - 1):max_time)

y0_vax <- initial_values(mod = "vax", size_months = parms$size_months, N_sim = N_sim)
out_vax <- mod_vax(y0 = y0_vax, max_time = max_time, parms = parms, N_sim = N_sim)
inc_vax <- extract_incidence(mod = "vax", out = out_vax, times = (max_time - 1):max_time)

y0_mab <- initial_values(mod = "mab", size_months = parms$size_months, N_sim = N_sim)
out_mab <- mod_mab(y0 = y0_mab, max_time = max_time, parms = parms, N_sim = N_sim)
inc_mab <- extract_incidence(mod = "mab", out = out_mab, times = (max_time - 1):max_time)

