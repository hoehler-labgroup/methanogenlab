#' Calculates solution alkalinity
#'
#' `alkalinity()` Uses the calculate_closed_system_alkalinity function from Sebastian Kopf's microbialkitchen R package to calculate media alkalinity
#'
#' @param pH.initial Initial pH
#' @param nDIC Moles of dissolved inorganic carbon
#' @param VolumeSolution Volume of liquid in the closed system, in liters.
#' @param VolumeHeadspace Volume of gaseous headspace in the closed system, in liters.
#' @param temperature Temperature of the system, in Kelvin.
#' @param K.CO2HCO3 Equilibrium constant for the dissociation of CO2(aq) to HCO3-(aq).
#' @param KHCO3CO3 Equilibrium constant for the dissociation of HCO3- (aq) to CO3-- (aq).
#' @return alkalinity of the solution, in millimolar.

alkalinity <- function(pH, nDIC, VolumeSolution, VolumeHeadspace, temperature, K.CO2HCO3, K.HCO3CO3){

  alkalinity <- as.numeric(scale_metric(calculate_closed_system_alkalinity(
    pH = pH,
    TIC = qty(nDIC, "mol"),
    V_liquid = qty(VolumeSolution, "L"),
    V_gas = qty(VolumeHeadspace, "L"),
    solubility = calculate_gas_solubility("CO2", qty((temperature-273.15), "C")),
    temperature = qty((temperature-273.15), "C"),
    pKa1 = -log10(K.CO2HCO3),
    pKa2 = -log10(K.HCO3CO3),
    pKw = 14,
    buffer = qty(0, "M"),
    buffer_pKa = 0
  ), prefix = "m"))

  return(alkalinity)
}
