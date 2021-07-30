
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

pressure.step <- function(n.total.step, K, VolumeHeadspace, temperature, VolumeSolution){ #calculates pressure of headspace in gas
  P.step <- n.total.step / ((VolumeHeadspace / (0.08206 * temperature)) + (K * VolumeSolution))
  return(P.step)
}
