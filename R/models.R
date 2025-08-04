
#' @title mod_base
#' @description Wrapper function for base model
#' @import deSolve
#' @param y0 Matrix of initial values of the state variables (rows are age groups and columns are variables)
#' @param times Vector of times to evaluate model
#' @param parms Named list of model parameters
#' @param N_sim Number of simulations (draws from the parameter prior distributions)
#' @return Matrix of outputs (rows are times and columns are variables)
#' @name mod_base
#' @export
mod_base <- function(y0, times, parms, N_sim){
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

  out <- array(NA, dim = c(N_sim, max(times), 75, 5),
               dimnames = list("simulation" = 1:N_sim,
                               "time" = 1:max(times),
                               "age" = parms$age_vect_years,
                               "variable" = c("S", "E", "I", "R", "Incidence")))
  for(sim in 1:N_sim){
    y0_tmp <- y0
    parms_tmp <- list(size_months = parms$size_months,
                      mixing      = parms$mixing,
                      b0          = parms$b0,
                      b1          = parms$b1,
                      phi         = parms$phi,
                      omega       = parms$omega,
                      delta       = parms$delta[sim],
                      gamma       = parms$gamma[sim],
                      nu          = parms$nu[sim],
                      sigma       = c(1 - exp(-parms$r_sigma[sim]*(1:12)), rep(1, 75 - 12)),
                      rho_V       = parms$rho_V[sim],
                      p_vax       = parms$p_vax[sim],
                      kappa_V     = parms$kappa_V[sim])
    for(tt in 1:max(times)){
      mod_out <- ode(y = y0_tmp,
                     times = seq(from = tt - 1, to = tt, by = 0.25),
                     func = deSolve_func,
                     parms = parms_tmp)
      out[sim,tt,,] <- matrix(as.vector(mod_out[nrow(mod_out),-1]), nrow = 75, ncol = 5)

      y0_tmp <- matrix(0, nrow = 75, ncol = 4)
      y0_tmp[1, 1] <- parms_tmp$size_months[1]/12
      y0_tmp[2:60, 1:4] <- out[sim, tt, 1:59, 1:4]
      y0_tmp[61, 1:4] <- out[sim, tt, 60, 1:4] + 59/60*out[sim, tt, 61, 1:4]
      y0_tmp[62:75, 1:4] <- 1/60*out[sim, tt, 61:(75 - 1), 1:4] + 59/60*out[sim, tt, 62:75, 1:4]
      y0_tmp <- cbind(y0_tmp, matrix(0, nrow = 75, ncol = 1))
    }
  }

  out <- lapply(1:N_sim, function(sim) cbind(1:max(times), t(apply(out[sim,,,], 1, as.vector))))
  return(out)
}

#' @title mod_vax
#' @description Wrapper function for vaccine model
#' @import deSolve
#' @param y0 Matrix of initial values of the state variables (rows are age groups and columns are variables)
#' @param times Vector of times to evaluate model
#' @param parms Named list of model parameters
#' @param N_sim Number of simulations (draws from the parameter prior distributions)
#' @return Matrix of outputs (rows are times and columns are variables)
#' @name mod_vax
#' @export
mod_vax <- function(y0, times, parms, N_sim){
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

  out <- array(NA, dim = c(N_sim, max(times), 75, 6),
               dimnames = list("simulation" = 1:N_sim,
                               "time" = 1:max(times),
                               "age" = parms$age_vect_years,
                               "variable" = c("S", "E", "I", "R", "V", "Incidence")))
  for(sim in 1:N_sim){
    y0_tmp <- y0
    parms_tmp <- list(size_months = parms$size_months,
                      mixing      = parms$mixing,
                      b0          = parms$b0,
                      b1          = parms$b1,
                      phi         = parms$phi,
                      omega       = parms$omega,
                      delta       = parms$delta[sim],
                      gamma       = parms$gamma[sim],
                      nu          = parms$nu[sim],
                      sigma       = c(1 - exp(-parms$r_sigma[sim]*(1:12)), rep(1, 75 - 12)),
                      rho_V       = parms$rho_V[sim],
                      p_vax       = parms$p_vax[sim],
                      kappa_V     = parms$kappa_V[sim])
    for(tt in 1:max(times)){
      mod_out <- ode(y = y0_tmp,
                     times = seq(from = tt - 1, to = tt, by = 0.25),
                     func = deSolve_func,
                     parms = parms_tmp)
      out[sim,tt,,] <- matrix(as.vector(mod_out[nrow(mod_out),-1]), nrow = 75, ncol = 6)

      y0_tmp <- matrix(0, nrow = 75, ncol = 5)
      y0_tmp[1, 1] <- (1 - parms_tmp$kappa_V)*parms_tmp$size_months[1]/12
      y0_tmp[1, 5] <- parms_tmp$kappa_V*parms_tmp$size_months[1]/12
      y0_tmp[2:parms_tmp$p_vax, 1:5] <- out[sim, tt, 1:(parms_tmp$p_vax - 1), 1:5]
      y0_tmp[parms_tmp$p_vax + 1, 1] <- out[sim, tt, parms_tmp$p_vax, 1] + out[sim, tt, parms_tmp$p_vax, 5]
      y0_tmp[parms_tmp$p_vax + 1, 2:4] <- out[sim, tt, parms_tmp$p_vax, 2:4]
      y0_tmp[(parms_tmp$p_vax + 2):60, 1:4] <- out[sim, tt, (parms_tmp$p_vax + 1):59, 1:4]
      y0_tmp[61, 1:4] <- out[sim, tt, 60, 1:4] + 59/60*out[sim, tt, 61, 1:4]
      y0_tmp[62:75, 1:4] <- 1/60*out[sim, tt, 61:(75 - 1), 1:4] + 59/60*out[sim, tt, 62:75, 1:4]
      y0_tmp <- cbind(y0_tmp, matrix(0, nrow = 75, ncol = 1))
    }
  }

  out <- lapply(1:N_sim, function(sim) cbind(1:max(times), t(apply(out[sim,,,], 1, as.vector))))
  return(out)
}

#' @title mod_mab
#' @description Wrapper function for monoclonal antibody model
#' @import deSolve
#' @param y0 Matrix of initial values of the state variables (rows are age groups and columns are variables)
#' @param times Vector of times to evaluate model
#' @param parms Named list of model parameters
#' @param N_sim Number of simulations (draws from the parameter prior distributions)
#' @return Matrix of outputs (rows are times and columns are variables)
#' @name mod_mab
#' @export
mod_mab <- function(y0, times, parms, N_sim){
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
      dV <- -infectM
      dinc <- infectS + infectM

      return(list(c(dS, dE, dI, dR, dM, dinc)))
    })
  }

  out <- array(NA, dim = c(N_sim, max(times), 75, 6),
               dimnames = list("simulation" = 1:N_sim,
                               "time" = 1:max(times),
                               "age" = parms$age_vect_years,
                               "variable" = c("S", "E", "I", "R", "M", "Incidence")))
  for(sim in 1:N_sim){
    y0_tmp <- y0
    parms_tmp <- list(size_months = parms$size_months,
                      mixing      = parms$mixing,
                      b0          = parms$b0,
                      b1          = parms$b1,
                      phi         = parms$phi,
                      omega       = parms$omega,
                      delta       = parms$delta[sim],
                      gamma       = parms$gamma[sim],
                      nu          = parms$nu[sim],
                      sigma       = c(1 - exp(-parms$r_sigma[sim]*(1:12)), rep(1, 75 - 12)),
                      rho_M       = parms$rho_M[sim],
                      p_mab       = parms$p_mab[sim],
                      kappa_M     = parms$kappa_M[sim])
    for(tt in 1:max(times)){
      mod_out <- ode(y = y0_tmp,
                     times = seq(from = tt - 1, to = tt, by = 0.25),
                     func = deSolve_func,
                     parms = parms_tmp)
      out[sim,tt,,] <- matrix(as.vector(mod_out[nrow(mod_out),-1]), nrow = 75, ncol = 6)

      y0_tmp <- matrix(0, nrow = 75, ncol = 5)
      y0_tmp[1, 1] <- (1 - parms_tmp$kappa_M)*parms_tmp$size_months[1]/12
      y0_tmp[1, 5] <- parms_tmp$kappa_M*parms_tmp$size_months[1]/12
      y0_tmp[2:parms_tmp$p_mab, 1:5] <- out[sim, tt, 1:(parms_tmp$p_mab - 1), 1:5]
      y0_tmp[parms_tmp$p_mab + 1, 1] <- out[sim, tt, parms_tmp$p_mab, 1] + out[sim, tt, parms_tmp$p_mab, 5]
      y0_tmp[parms_tmp$p_mab + 1, 2:4] <- out[sim, tt, parms_tmp$p_mab, 2:4]
      y0_tmp[(parms_tmp$p_mab + 2):60, 1:4] <- out[sim, tt, (parms_tmp$p_mab + 1):59, 1:4]
      y0_tmp[61, 1:4] <- out[sim, tt, 60, 1:4] + 59/60*out[sim, tt, 61, 1:4]
      y0_tmp[62:75, 1:4] <- 1/60*out[sim, tt, 61:(75 - 1), 1:4] + 59/60*out[sim, tt, 62:75, 1:4]
      y0_tmp <- cbind(y0_tmp, matrix(0, nrow = 75, ncol = 1))
    }
  }

  out <- lapply(1:N_sim, function(sim) cbind(1:max(times), t(apply(out[sim,,,], 1, as.vector))))
  return(out)
}
