#' Calculates CO2 pressure in the head space of the closed system.
#'
#' `PCO2()` Uses the calculate_closed_system_pCO2 function from Sebastian Kopf's microbialkitchen R package to calculate pressure of CO2 gas.
#'
#' @param pH pH of the solution.
#' @param nDIC Moles of dissolved inorganic carbon.
#' @param VolumeSolution Volume of liquid in the closed system, in liters.
#' @param VolumeHeadspace Volume of gaseous head space in the closed system, in liters.
#' @param temperature Temperature of the system, in Kelvin.
#' @param K.CO2HCO3 Equilibrium constant for the dissociation of CO2(aq) to HCO3-(aq). 5.223196e-07 by default.
#' @param KHCO3CO3 Equilibrium constant for the dissociation of HCO3- (aq) to CO3-- (aq). 6.01886e-11 by default.
#' @return CO2 pressure in the head space, in atm.
#'
PCO2 <- function(pH, nDIC, VolumeSolution, VolumeHeadspace, temperature, K.CO2HCO3 = 5.223196e-07, K.HCO3CO3 = 6.01886e-11){

  PCO2 <- as.numeric(scale_metric(calculate_closed_system_pCO2(
    pH = pH,
    TIC = qty(nDIC, "mol"),
    V_liquid = qty(VolumeSolution, "L"),
    V_gas = qty(VolumeHeadspace, "L"),
    solubility = calculate_gas_solubility("CO2", qty((temperature-273.15), "C")),
    temperature = qty((temperature-273.15), "C"),
    pKa1 = -log10(K.CO2HCO3),
    pKa2 = -log10(K.HCO3CO3)
  ), prefix = "m"))
  return(PCO2/1013.25)
}

#' Calculates head space pressure
#'
#' `pressure.step()` calculates head space partial pressure from the total moles of a particular gas in the system using Henry's law and the Ideal Gas law.
#'
#' @param n.total.step The total number of moles of a gas in the system.
#' @param K Henry's constant of the gas.
#' @param VolumeSolution Volume of liquid in the closed system, in liters.
#' @param VolumeHeadspace Volume of gaseous head space in the closed system, in liters.
#' @param temperature Temperature of the system, in Kelvin.
#' @return Partial pressure of the head space.
#'
pressure.step <- function(n.total.step, K, VolumeSolution, VolumeHeadspace, temperature){
  P.step <- n.total.step / ((VolumeHeadspace / (0.082057366 * temperature)) + (K * VolumeSolution))
  return(P.step)
}
