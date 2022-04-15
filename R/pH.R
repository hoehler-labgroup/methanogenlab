#' Calculates pH of the solution
#'
#' `pH()` Uses the calculate_closed_system_pH function from Sebastian Kopf's microbialkitchen R package to calculate pH.
#'
#' @param alkalinity initial alkalinity of the system, in millimolar.
#' @inheritParams alkalinity
#' @return pH of the solution.

pH <- function(nDIC, VolumeSolution, VolumeHeadspace, temperature, alkalinity, K.CO2HCO3, K.HCO3CO3){

  systempH <- as.numeric(microbialkitchen::calculate_closed_system_pH(
    TIC = qty(nDIC, "mol"),
    V_liquid = qty(VolumeSolution, "L"),
    V_gas = qty(VolumeHeadspace, "L"),
    solubility = microbialkitchen::calculate_gas_solubility("CO2", qty((temperature-273.15), "C")),
    temperature = qty((temperature-273.15), "C"),
    pKa1 = -log10(K.CO2HCO3),
    pKa2 = -log10(K.HCO3CO3),
    pKw = 14,
    buffer = qty(0, "M"),
    buffer_pKa = 0,
    alkalinity = qty(alkalinity, "mM")
  ))

  return(systempH)
}
