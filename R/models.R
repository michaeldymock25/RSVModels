
#' @title mod_base
#' @description Wrapper function for base model
#' @import deSolve
#' @import parallel
#' @param y0 Matrix of initial values of the state variables (rows are age groups and columns are variables)
#' @param max_time Maximum time in months to run the model (includes burn in)
#' @param parms Named list of model parameters
#' @param N_sim Number of simulations (draws from the parameter prior distributions)
#' @param batch_size Size of each batch for parallel simulations. Default is 100.
#' @param ncores Number of cores to run in parallel. Default is one.
#' @param thin Integer used to thin stored outputs. Must be a factor of max_time. Note that this will not influence the transmission modelling. Default is one.
#' @param minimal Logical indicator. If set to TRUE, then only minimal output is stored (incidence only). Default is FALSE.
#' @return Array of named outputs containing simulations, times, ages and states
#' @name mod_base
#' @export
mod_base <- function(y0, max_time, parms, N_sim, batch_size = 100, ncores = 1, thin = 1, minimal = FALSE){
  deSolve_func <- function(t, y, parms){
    y <- matrix(y, nrow = 75, ncol = 5)

    with(as.list(c(y, parms)),{
      S <- y[, 1]
      E <- y[, 2]
      I <- y[, 3]
      R <- y[, 4]
      N <- S + E + I + R

      temp <- omega*I/N
      s <- sweep(mixing, MARGIN = 2, temp, FUN = "*")
      lambda <- b0*(1 + b1*cos(2*pi*t/12 + phi))*rowSums(s)

      infect <- lambda*sigma*S
      dS <- nu*R - infect
      dE <- infect - delta*E
      dI <- delta*E - gamma*I
      dR <- gamma*I - nu*R
      dinc <- infect

      return(list(c(dS, dE, dI, dR, dinc)))
    })
  }

  N_batch <- ceiling(N_sim/batch_size)
  batch_lens <- sapply(1:N_batch, function(i) ifelse(i < N_batch, batch_size, N_sim - (N_batch - 1)*batch_size))

  out_l <- mclapply(1:N_batch, function(i){
    if(minimal){
      tmp <- array(NA, dim = c(batch_lens[i], max_time/thin, 75, 1),
                   dimnames = list("simulation" = 1:batch_lens[i],
                                   "time" = 1:(max_time/thin),
                                   "age" = parms$age_years,
                                   "variable" = "Incidence"))
    } else {
      tmp <- array(NA, dim = c(batch_lens[i], max_time/thin, 75, 5),
                   dimnames = list("simulation" = 1:batch_lens[i],
                                   "time" = 1:(max_time/thin),
                                   "age" = parms$age_years,
                                   "variable" = c("S", "E", "I", "R", "Incidence")))
    }
    sim_ref <- (sum(c(0,batch_lens)[1:i]) + 1):sum(batch_lens[1:i])
    for(sim in 1:batch_lens[i]){
      y0_tmp <- y0[sim_ref[sim],,]
      parms_tmp <- list(size_months = parms$size_months,
                        mixing      = parms$mixing,
                        b0          = parms$b0,
                        b1          = parms$b1,
                        phi         = parms$phi,
                        omega       = parms$omega,
                        delta       = parms$delta[sim_ref[sim]],
                        gamma       = parms$gamma[sim_ref[sim]],
                        nu          = parms$nu[sim_ref[sim]],
                        sigma       = c(1 - exp(-parms$r_sigma[sim_ref[sim]]*(1:12)), rep(1, 75 - 12)))
      for(tt in 1:max_time){
        mod_out <- ode(y = y0_tmp,
                       times = seq(from = 0, to = 1, by = 0.25),
                       func = deSolve_func,
                       parms = parms_tmp)
        mod_out <- matrix(as.vector(mod_out[nrow(mod_out),-1]), nrow = 75, ncol = 5)
        if(minimal){
          if(tt %% thin == 0) tmp[sim,tt/thin,,] <- mod_out[,5]
        } else {
          if(tt %% thin == 0) tmp[sim,tt/thin,,] <- mod_out
        }
        y0_tmp <- matrix(0, nrow = 75, ncol = 4)
        y0_tmp[1, 1] <- parms_tmp$size_months[1]/12
        y0_tmp[2:60, 1:4] <- mod_out[1:59, 1:4]
        y0_tmp[61, 1:4] <- mod_out[60, 1:4] + 59/60*mod_out[61, 1:4]
        y0_tmp[62:75, 1:4] <- 1/60*mod_out[61:(75 - 1), 1:4] + 59/60*mod_out[62:75, 1:4]
        y0_tmp <- cbind(y0_tmp, matrix(0, nrow = 75, ncol = 1))
      }
    }
    return(tmp)
  }, mc.cores = ncores)
  if(minimal){
    out <- array(NA, dim = c(N_sim, max_time/thin, 75, 1),
                 dimnames = list("simulation" = 1:N_sim,
                                 "time" = 1:(max_time/thin),
                                 "age" = parms$age_years,
                                 "variable" = "Incidence"))
  } else {
    out <- array(NA, dim = c(N_sim, max_time/thin, 75, 5),
                 dimnames = list("simulation" = 1:N_sim,
                                 "time" = 1:(max_time/thin),
                                 "age" = parms$age_years,
                                 "variable" = c("S", "E", "I", "R", "Incidence")))
  }
  for(i in 1:N_batch) out[(sum(c(0,batch_lens)[1:i]) + 1):sum(batch_lens[1:i]),,,] <- out_l[[i]]
  gc()

  return(out)
}

#' @title mod_vax
#' @description Wrapper function for vaccine model
#' @import deSolve
#' @import parallel
#' @param y0 Matrix of initial values of the state variables (rows are age groups and columns are variables)
#' @param max_time Maximum time in months to run the model (includes burn in)
#' @param parms Named list of model parameters
#' @param N_sim Number of simulations (draws from the parameter prior distributions)
#' @param batch_size Size of each batch for parallel simulations. Default is 100.
#' @param ncores Number of cores to run in parallel. Default is one.
#' @param thin Integer used to thin stored outputs. Must be a factor of max_time. Note that this will not influence the transmission modelling. Default is one.
#' @param minimal Logical indicator. If set to TRUE, then only minimal output is stored (vaccine recipients and incidence only). Default is FALSE.
#' @return Array of named outputs containing simulations, times, ages and states
#' @name mod_vax
#' @export
mod_vax <- function(y0, max_time, parms, N_sim, batch_size = 100, ncores = 1, thin = 1, minimal = FALSE){
  deSolve_func <- function(t, y, parms){
    y <- matrix(y, nrow = 75, ncol = 6)

    with(as.list(c(y, parms)),{
      S <- y[, 1]
      E <- y[, 2]
      I <- y[, 3]
      R <- y[, 4]
      V <- y[, 5]
      N <- S + E + I + R + V

      temp <- omega*I/N
      s <- sweep(mixing, MARGIN = 2, temp, FUN = "*")
      lambda <- b0*(1 + b1*cos(2*pi*t/12 + phi))*rowSums(s)

      infectS <- lambda*sigma*S
      infectV <- lambda*(1-rho_V)*V
      dS <- nu*R - infectS
      dE <- infectS + infectV - delta*E
      dI <- delta*E - gamma*I
      dR <- gamma*I - nu*R
      dV <- -infectV
      dinc <- infectS + infectV

      return(list(c(dS, dE, dI, dR, dV, dinc)))
    })
  }

  N_batch <- ceiling(N_sim/batch_size)
  batch_lens <- sapply(1:N_batch, function(i) ifelse(i < N_batch, batch_size, N_sim - (N_batch - 1)*batch_size))

  out_l <- mclapply(1:N_batch, function(i){
    if(minimal){
      tmp <- array(NA, dim = c(batch_lens[i], max_time/thin, 75, 2),
                   dimnames = list("simulation" = 1:batch_lens[i],
                                   "time" = 1:(max_time/thin),
                                   "age" = parms$age_years,
                                   "variable" = c("V", "Incidence")))
    } else {
      tmp <- array(NA, dim = c(batch_lens[i], max_time/thin, 75, 6),
                   dimnames = list("simulation" = 1:batch_lens[i],
                                   "time" = 1:(max_time/thin),
                                   "age" = parms$age_years,
                                   "variable" = c("S", "E", "I", "R", "V", "Incidence")))
    }
    sim_ref <- (sum(c(0,batch_lens)[1:i]) + 1):sum(batch_lens[1:i])
    for(sim in 1:batch_lens[i]){
      y0_tmp <- y0[sim_ref[sim],,]
      parms_tmp <- list(size_months = parms$size_months,
                        mixing      = parms$mixing,
                        b0          = parms$b0,
                        b1          = parms$b1,
                        phi         = parms$phi,
                        omega       = parms$omega,
                        delta       = parms$delta[sim_ref[sim]],
                        gamma       = parms$gamma[sim_ref[sim]],
                        nu          = parms$nu[sim_ref[sim]],
                        sigma       = c(1 - exp(-parms$r_sigma[sim_ref[sim]]*(1:12)), rep(1, 75 - 12)),
                        rho_V       = parms$rho_V[sim_ref[sim]],
                        dur_V       = parms$dur_V[sim_ref[sim]],
                        kappa_V     = parms$kappa_V[sim_ref[sim]])
      for(tt in 1:max_time){
        mod_out <- ode(y = y0_tmp,
                       times = seq(from = 0, to = 1, by = 0.25),
                       func = deSolve_func,
                       parms = parms_tmp)
        mod_out <- matrix(as.vector(mod_out[nrow(mod_out),-1]), nrow = 75, ncol = 6)
        if(minimal){
          if(tt %% thin == 0) tmp[sim,tt/thin,,] <- mod_out[,c(5,6)]
        } else {
          if(tt %% thin == 0) tmp[sim,tt/thin,,] <- mod_out
        }
        y0_tmp <- matrix(0, nrow = 75, ncol = 5)
        y0_tmp[1, 1] <- (1 - parms_tmp$kappa_V)*parms_tmp$size_months[1]/12
        y0_tmp[1, 5] <- parms_tmp$kappa_V*parms_tmp$size_months[1]/12
        y0_tmp[2:parms_tmp$dur_V, 1:5] <- mod_out[1:(parms_tmp$dur_V - 1), 1:5]
        y0_tmp[parms_tmp$dur_V + 1, 1] <- mod_out[parms_tmp$dur_V, 1] + mod_out[parms_tmp$dur_V, 5]
        y0_tmp[parms_tmp$dur_V + 1, 2:4] <- mod_out[parms_tmp$dur_V, 2:4]
        y0_tmp[(parms_tmp$dur_V + 2):60, 1:4] <- mod_out[(parms_tmp$dur_V + 1):59, 1:4]
        y0_tmp[61, 1:4] <- mod_out[60, 1:4] + 59/60*mod_out[61, 1:4]
        y0_tmp[62:75, 1:4] <- 1/60*mod_out[61:(75 - 1), 1:4] + 59/60*mod_out[62:75, 1:4]
        y0_tmp <- cbind(y0_tmp, matrix(0, nrow = 75, ncol = 1))
      }
    }
    return(tmp)
  }, mc.cores = ncores)
  if(minimal){
    out <- array(NA, dim = c(N_sim, max_time/thin, 75, 2),
                 dimnames = list("simulation" = 1:N_sim,
                                 "time" = 1:(max_time/thin),
                                 "age" = parms$age_years,
                                 "variable" = c("V", "Incidence")))
  } else {
    out <- array(NA, dim = c(N_sim, max_time/thin, 75, 6),
                 dimnames = list("simulation" = 1:N_sim,
                                 "time" = 1:(max_time/thin),
                                 "age" = parms$age_years,
                                 "variable" = c("S", "E", "I", "R", "V", "Incidence")))
  }
  for(i in 1:N_batch) out[(sum(c(0,batch_lens)[1:i]) + 1):sum(batch_lens[1:i]),,,] <- out_l[[i]]
  gc()

  return(out)
}

#' @title mod_mab
#' @description Wrapper function for monoclonal antibody model
#' @import deSolve
#' @import parallel
#' @param y0 Matrix of initial values of the state variables (rows are age groups and columns are variables)
#' @param max_time Maximum time in months to run the model (includes burn in)
#' @param parms Named list of model parameters
#' @param N_sim Number of simulations (draws from the parameter prior distributions)
#' @param batch_size Size of each batch for parallel simulations. Default is 100.
#' @param ncores Number of cores to run in parallel. Default is one.
#' @param thin Integer used to thin stored outputs. Must be a factor of max_time. Note that this will not influence the transmission modelling. Default is one.
#' @param minimal Logical indicator. If set to TRUE, then only minimal output is stored (monoclonal antibody recipients and incidence only). Default is FALSE.
#' @return Array of named outputs containing simulations, times, ages and states.
#' @name mod_mab
#' @export
mod_mab <- function(y0, max_time, parms, N_sim, batch_size = 100, ncores = 1, thin = 1, minimal = FALSE){
  deSolve_func <- function(t, y, parms){
    y <- matrix(y, nrow = 75, ncol = 6)

    with(as.list(c(y, parms)),{
      S <- y[, 1]
      E <- y[, 2]
      I <- y[, 3]
      R <- y[, 4]
      M <- y[, 5]
      N <- S + E + I + R + M

      temp <- omega*I/N
      s <- sweep(mixing, MARGIN = 2, temp, FUN = "*")
      lambda <- b0*(1 + b1*cos(2*pi*t/12 + phi))*rowSums(s)

      infectS <- lambda*sigma*S
      infectM <- lambda*(1-rho_M)*M
      dS <- nu*R - infectS
      dE <- infectS + infectM - delta*E
      dI <- delta*E - gamma*I
      dR <- gamma*I - nu*R
      dM <- -infectM
      dinc <- infectS + infectM

      return(list(c(dS, dE, dI, dR, dM, dinc)))
    })
  }

  N_batch <- ceiling(N_sim/batch_size)
  batch_lens <- sapply(1:N_batch, function(i) ifelse(i < N_batch, batch_size, N_sim - (N_batch - 1)*batch_size))

  out_l <- mclapply(1:N_batch, function(i){
    if(minimal){
      tmp <- array(NA, dim = c(batch_lens[i], max_time/thin, 75, 2),
                   dimnames = list("simulation" = 1:batch_lens[i],
                                   "time" = 1:(max_time/thin),
                                   "age" = parms$age_years,
                                   "variable" = c("M", "Incidence")))
    } else {
      tmp <- array(NA, dim = c(batch_lens[i], max_time/thin, 75, 6),
                   dimnames = list("simulation" = 1:batch_lens[i],
                                   "time" = 1:(max_time/thin),
                                   "age" = parms$age_years,
                                   "variable" = c("S", "E", "I", "R", "M", "Incidence")))
    }
    sim_ref <- (sum(c(0,batch_lens)[1:i]) + 1):sum(batch_lens[1:i])
    for(sim in 1:batch_lens[i]){
      y0_tmp <- y0[sim_ref[sim],,]
      parms_tmp <- list(size_months = parms$size_months,
                        mixing      = parms$mixing,
                        b0          = parms$b0,
                        b1          = parms$b1,
                        phi         = parms$phi,
                        omega       = parms$omega,
                        delta       = parms$delta[sim_ref[sim]],
                        gamma       = parms$gamma[sim_ref[sim]],
                        nu          = parms$nu[sim_ref[sim]],
                        sigma       = c(1 - exp(-parms$r_sigma[sim_ref[sim]]*(1:12)), rep(1, 75 - 12)),
                        rho_M       = parms$rho_M[sim_ref[sim]],
                        dur_M       = parms$dur_M[sim_ref[sim]],
                        kappa_M     = parms$kappa_M[sim_ref[sim]])
      for(tt in 1:max_time){
        mod_out <- ode(y = y0_tmp,
                       times = seq(from = 0, to = 1, by = 0.25),
                       func = deSolve_func,
                       parms = parms_tmp)
        mod_out <- matrix(as.vector(mod_out[nrow(mod_out),-1]), nrow = 75, ncol = 6)
        if(minimal){
          if(tt %% thin == 0) tmp[sim,tt/thin,,] <- mod_out[,c(5,6)]
        } else {
          if(tt %% thin == 0) tmp[sim,tt/thin,,] <- mod_out
        }
        y0_tmp <- matrix(0, nrow = 75, ncol = 5)
        y0_tmp[1, 1] <- (1 - parms_tmp$kappa_M)*parms_tmp$size_months[1]/12
        y0_tmp[1, 5] <- parms_tmp$kappa_M*parms_tmp$size_months[1]/12
        y0_tmp[2:parms_tmp$dur_M, 1:5] <- mod_out[1:(parms_tmp$dur_M - 1), 1:5]
        y0_tmp[parms_tmp$dur_M + 1, 1] <- mod_out[parms_tmp$dur_M, 1] + mod_out[parms_tmp$dur_M, 5]
        y0_tmp[parms_tmp$dur_M + 1, 2:4] <- mod_out[parms_tmp$dur_M, 2:4]
        y0_tmp[(parms_tmp$dur_M + 2):60, 1:4] <- mod_out[(parms_tmp$dur_M + 1):59, 1:4]
        y0_tmp[61, 1:4] <- mod_out[60, 1:4] + 59/60*mod_out[61, 1:4]
        y0_tmp[62:75, 1:4] <- 1/60*mod_out[61:(75 - 1), 1:4] + 59/60*mod_out[62:75, 1:4]
        y0_tmp <- cbind(y0_tmp, matrix(0, nrow = 75, ncol = 1))
      }
    }
    return(tmp)
  }, mc.cores = ncores)
  if(minimal){
    out <- array(NA, dim = c(N_sim, max_time/thin, 75, 2),
                 dimnames = list("simulation" = 1:N_sim,
                                 "time" = 1:(max_time/thin),
                                 "age" = parms$age_years,
                                 "variable" = c("M", "Incidence")))
  } else {
    out <- array(NA, dim = c(N_sim, max_time/thin, 75, 6),
                 dimnames = list("simulation" = 1:N_sim,
                                 "time" = 1:(max_time/thin),
                                 "age" = parms$age_years,
                                 "variable" = c("S", "E", "I", "R", "M", "Incidence")))
  }
  for(i in 1:N_batch) out[(sum(c(0,batch_lens)[1:i]) + 1):sum(batch_lens[1:i]),,,] <- out_l[[i]]
  gc()

  return(out)
}
