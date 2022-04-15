#' Biomass to carbon yield
#'
#' `carbon.yield()` converts the known organism-specific biomass yield per mol of product into carbon yield.
#' @param biomass.yield mass of dry biomass produced per mol of product.
#' @param carbon.fraction w/w percent C of biomass, expressed as a decimal.
#' @return carbon yield, in mol C/mol product
carbon.yield <- function(biomass.yield, carbon.fraction){
  carbon.biomass.per.mol = biomass.yield*carbon.fraction*1/12.0107
  return(carbon.biomass.per.mol)
}

#' Amount of biomass produced for every DIC consumed
#'
#' `biomass.coefficient()` calculates the amount of biomass produced for every mol of DIC consumed
#' @inheritParams carbon.yield
#' @return biomass produced, in mol/mol DIC
biomass.coefficient <- function(biomass.yield,carbon.fraction){
  energy <- 1
  biomass <- carbon.yield(biomass.yield,carbon.fraction)
  return(biomass/(energy+biomass))
}

#' Amount of CH4 produced for every DIC consumed
#'
#' `CH4.coefficient()` calculates the amount of CH4 produced for every mol of DIC consumed.
#' @inheritParams carbon.yield
#' @return CH4 produced, in mol/mol DIC
CH4.coefficient <- function(biomass.yield, carbon.fraction){
  rxn.coefficient <- 1
  biomass <- carbon.yield(biomass.yield,carbon.fraction)
  return(rxn.coefficient/(rxn.coefficient+biomass))
}

#' Amount of H2 consumed for every DIC consumed
#'
#' `H2.coefficient()` calculates the amount of H2 consumed for every mol of DIC consumed.
#' @inheritParams methanogenesis.time
#' @return H2 consumed, in mol/mol DIC
H2.coefficient <- function(biomass.yield, carbon.fraction){
  rxn.coefficient <- 4
  biomass.rxn.coefficient <- 2.1
  biomass.contribution <- biomass.rxn.coefficient*biomass.coefficient(biomass.yield,carbon.fraction)
  energy.contribution <- rxn.coefficient*CH4.coefficient(biomass.yield, carbon.fraction)
  return(biomass.contribution+energy.contribution)
}
