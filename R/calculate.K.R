#' Calculates Henry's constant for a gas
#'
#' `calculate.Kgas()` uses the provided Bunsen coefficient and moles per volume of gas to calculate Henry's constant.
#' @param bunsen The Bunsen coefficient for the gas.
#' @param molar.vol.gas The moles per unit volume of gas.
#' @export
calculate.Kgas <- function(bunsen,T,P){
  molar.vol.gas <- 0.082057366*T/P
  Kgas <- bunsen / molar.vol.gas
  return(Kgas)
}

#' Calculates equilibrium constants using CHNOSZ
#'
#' `calculate.KH()` utilizes the subcrt function from Jeffrey Dick's CHNOSZ R package to calculate equilibrium constants.
#' @param reactants A vector of all components involved in the reaction, both reactants and products.
#' @param moles A vector of the molar coefficients for the reaction, with negative values indicating reactants and positive values indicating products.
#' @param phases A vector of the phases for all components in the reaction, either "aq", "l", or "g".
#' @param temperature Temperature of the system, in Kelvin.
#' @param pressure The pressure of the system, in atm.
#' @details
#' For the following chemical equilibrium:
#' \loadmathjax
#' \mjdeqn{CO_2(aq)+H_2O(l) \leftrightharpoons {HCO_3}^{2-}(aq)+H^+(aq)
#' }{ASCII representation}
#' an example can found below
#' @examples
#' calculate_KH(c("CO2","H2O","HCO3-","H+"),c(-1,-1,1,1),c("aq","l","aq","aq"), 273.15+37,1.7)
#'
#' @export
calculate.KH <- function(reactants, moles, phases, temperature, pressure){
  KH = 10**(as.numeric(subcrt(reactants,moles,phases,T=temperature-273.15,P=pressure*1.01325)$out$logK))
  return(KH)
}

