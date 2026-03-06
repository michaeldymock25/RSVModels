
set.seed(6824767)

N_sim <- 10000
pop <- readRDS("Data/auspop.rds")
pop <- data.frame(age = pop$lower.age.limit, population = pop$population)
mixing <- as.matrix(read.csv("Data/mixing.csv", header = TRUE), header = TRUE)

parms <- list(total_pop   = sum(pop$population),
              age_years   = c(seq(0, 5, 1/12), seq(10, 75, 5)),
              age_months  = 12*c(seq(0, 5, 1/12), seq(10, 75, 5)),
              size_months = c(rep(pop$population[1]/60, 60), pop$population[-c(1,nrow(pop))]),
              mixing      = mixing,
              b0          = 0.087,
              b1          = -0.193,
              phi         = 1.536,
              omega       = c(rep(1.00, 59), rep(0.35, 16)),
              delta       = 1/(rgamma(N_sim, 16, 4)/(365/12)),
              gamma       = 1/(rgamma(N_sim, 81, 9)/(365/12)),
              nu          = 1/(rgamma(N_sim, 529, 2.3)/(365/12)),
              r_sigma     = runif(N_sim, 0.5, 1),
              rho_V       = rep(0.8, N_sim),
              p_vax       = sample(3:9, size = N_sim, replace = TRUE),
              kappa_V     = rbeta(N_sim, 10, 4.286),
              rho_M       = rep(0.8, N_sim),
              p_mab       = sample(3:9, size = N_sim, replace = TRUE),
              kappa_M     = rbeta(N_sim, 10, 1.765))

saveRDS(parms, file = "Data/parameters.rds")
