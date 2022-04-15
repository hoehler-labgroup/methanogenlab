#'Calculates Gibbs free energy
#'
#' `gibbs.step()` calculates the Gibbs free energy at each step in the
#' model for a given Q and standard Gibbs free energy.
#' @param standard.gibbs The standard Gibbs free energy of formation of the reaction in J/mol
#' @param Q Current reaction quotient
#' @inheritParams methanogenesis.time
gibbs.step <- function(standard.gibbs,Q, temperature){
  gibbs.step <- (standard.gibbs + 8.314462 * temperature * log(Q)) / 1000
  return(gibbs.step)
}

#' Calculates standard Gibbs free energy for a given temperature and pressure
#'
#'`standard.gibbs()` determines the standard Gibbs free energy to be used in `gibbs.step()`
#' @inheritParams calculate.KH
standard.gibbs <- function(reactants, moles, phases,temperature,pressure=1){
  standard.gibbs <- 4.184*CHNOSZ::subcrt(reactants,moles,phases,T=temperature-273.15,P=pressure*1.01325)$out$G
  return(standard.gibbs)
}
