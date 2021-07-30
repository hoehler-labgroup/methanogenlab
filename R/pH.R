
pH <- function(nDIC, VolumeSolution, VolumeHeadspace, temperature, alkalinity, K.CO2HCO3 = 5.223196e-07, K.HCO3CO3 = 6.01886e-11){
  
  systempH <- as.numeric(calculate_closed_system_pH(
    TIC = qty(nDIC, "mol"),
    V_liquid = qty(VolumeSolution, "L"),
    V_gas = qty(VolumeHeadspace, "L"),
    solubility = calculate_gas_solubility("CO2", qty((temperature-273.15), "C")),
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
