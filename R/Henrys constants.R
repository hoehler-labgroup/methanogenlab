#' Henry's constant for typical gas
#' @parm test
#' @export
calculate_Kgas <- function(bunsen, molar.vol.gas){
  Kgas <- bunsen / molar.vol.gas
  return(Kgas)
}

#' Henry's constant for CO2
#' @export
calculate_KH <- function(reactants, moles, phases, temperature, pressure){
  KH = 10^(as.numeric(subcrt(reactants,moles,phases,T= (temperature-273.15),P=pressure)[[2]][4:4]))
  return(KH)
}
