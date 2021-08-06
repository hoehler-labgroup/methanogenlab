#' Uses Henry's Law to calculate dissolved gas concentrations
#'
#' `aqueous.step()` calculates the aqueous concentration at each step using the current pressure.
#' @param P.step the current pressure in the system.
#' @param K Henry's equilibrium constant of the gas.
#' @return Dissolved gas concentration, in molar.
aqueous.step <- function(P.step, K){
  conc.step <- P.step * K
  return(conc.step)
}
