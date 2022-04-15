#' Calculates solution alkalinity
#'
#' `alkalinity()` Uses the calculate_closed_system_alkalinity function from Sebastian Kopf's microbialkitchen R package to calculate media alkalinity
#'
#' @inheritParams methanogenesis.time
#' @param pH pH of the system.
#' @param nDIC Moles of dissolved inorganic carbon
#' @param K.CO2HCO3 Equilibrium constant for the dissociation of CO2(aq) to HCO3-(aq).
#' @param K.HCO3CO3 Equilibrium constant for the dissociation of HCO3- (aq) to CO3– (aq).
#' @return alkalinity of the solution, in millimolar.

alkalinity <- function(pH, nDIC, VolumeSolution, VolumeHeadspace, temperature, K.CO2HCO3, K.HCO3CO3){

  alkalinity <- as.numeric(microbialkitchen::scale_metric(microbialkitchen::calculate_closed_system_alkalinity(
    pH = pH,
    TIC = qty(nDIC, "mol"),
    V_liquid = qty(VolumeSolution, "L"),
    V_gas = qty(VolumeHeadspace, "L"),
    solubility = microbialkitchen::calculate_gas_solubility("CO2", qty((temperature-273.15), "C")),
    temperature = qty((temperature-273.15), "C"),
    pKa1 = -log10(K.CO2HCO3),
    pKa2 = -log10(K.HCO3CO3),
    pKw = 14,
    buffer = qty(0, "M"),
    buffer_pKa = 0
  ), prefix = "m"))

  return(alkalinity)
}
