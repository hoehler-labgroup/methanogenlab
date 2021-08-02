#'Calculates Gibbs free energy
#'
#' `gibbs.step()` calculates the Gibbs free energy at each step in the
#' model for a given Q and standard Gibbs free energy.
#' @param standard.gibbs The standard Gibbs free energy of formation of the reaction in J/mol
#' @param Q Current reaction quotient
#' @param temperature Temperature of the system, in Kelvin.
gibbs.step <- function(standard.gibbs, Q, temperature){
  Gibbs.step <- (standard.gibbs + 8.314 * temperature * log(Q)) / 1000
  return(Gibbs.step)
}
