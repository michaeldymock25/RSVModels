
#' @title get_contact_matrix
#' @description Generates social contact mixing matrix based on the Australian population
#' @import conmat
#' @param pop Table containing population data. Must contain columns called "age" and "population" containing the lower age limit for and population size, respectively.
#' @param seed Optional parameter used for reproducibility. Defaults to 6824767.
#' @return Returns mixing matrix.
#' @name get_contact_matrix
#' @export
get_contact_matrix <- function(pop, seed = 6824767){
  set.seed(seed)
  suppressPackageStartupMessages(library("conmat", character.only = TRUE))
  pop_conmat <- conmat::as_conmat_population(pop, age = age, population = population)
  setting_models <- conmat::fit_setting_contacts(contact_data_list = conmat::get_polymod_setting_data(countries = "United Kingdom"),
                                                 population = conmat::get_polymod_population(countries = "United Kingdom"))
  mixing_agg <- conmat::predict_setting_contacts(population = pop_conmat, contact_model = setting_models, age_breaks = pop_conmat$age)
  mixing_agg <- mixing_agg$all
  mixing_agg_sym <- matrix(0, nrow = 16, ncol = 16)
  for(i in 1:16) for(j in 1:16) mixing_agg_sym[i,j] <- (mixing_agg[i,j] + mixing_agg[j, i])/2

  mixing <- matrix(0, ncol = 75, nrow = 75)
  colnames(mixing) <- rownames(mixing) <- c(seq(0, 5, 1/12), seq(10, 75, 5))
  mixing[1:60, 1:60] <- mixing_agg_sym[1,1]/60
  for(j in 1:60) mixing[61:75, j] <- mixing_agg_sym[2:16, 1]
  for(j in 61:75) mixing[1:60, j] <- mixing_agg_sym[1, j - 60 + 1]/60
  mixing[61:75, 61:75] <- mixing_agg_sym[2:16, 2:16]
  mixing <- t(mixing)
  mixing <- mixing*365.25/12

  return(mixing)
}

#' @title initial_values
#' @description Generates matrix of initial values to be used as input to a transmission model
#' @param mod The specific transmission model that requires initial values (base, vax or mab).
#' @param size_months Population size for each age group. Is a parameter in parms.
#' @param N_sim Number of simulations (draws from the parameter prior distributions)
#' @return Returns matrix of initial values. Rows are age groups and columns are states. Final column contains incidence (zero initially).
#' @name initial_values
#' @export
initial_values <- function(mod, size_months, N_sim){
  num_var <- ifelse(mod == "base", 5, ifelse(mod %in% c("vax", "mab"), 6, NA))
  y0 <- array(data = 0, dim = c(N_sim, 75, num_var))
  y0[,,1] <- rep(size_months*0.99, each = N_sim)
  y0[,,3] <- rep(size_months*0.01, each = N_sim)
  return(y0)
}

#' @title extract_incidence
#' @description Extracts incidence estimates from the model output for a given time period.
#' @param mod The specific transmission model that requires initial values (base, vax or mab).
#' @param out Transmission model output.
#' @param times Time period for incidence extraction
#' @param total Logical indictor. If TRUE then will sum over time period, otherwise will return disaggregated. Default is TRUE.
#' @return Returns matrix of incidence (rows are simulations and columns are age groups).
#' @name extract_incidence
#' @export
extract_incidence <- function(mod, out, times, total = TRUE){
  col <- ifelse(mod == "base", 5, 6)
  inc <- out[,times,,col]
  if(total) inc <- apply(inc, c(1,3), sum)
  return(inc)
}

