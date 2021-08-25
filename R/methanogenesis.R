#' Determines initial conditions from initial inputs
#'
#' `init()` sets up the initial environment to be used by the methanogenesis model.
#' @param CH4.initial Concentration of initial dissolved CH4, in molarity.
#' @param K.CH4 Henry's constant for CH4.
#' @param H2.initial Concentration of initial dissolved H2, in molarity.
#' @param K.H2 Henry's constant for H2.
#' @param DIC.initial Concentration of initial dissolved inorganic carbon, in molarity.
#' @param pH.initial initial pH.
#' @param K.CO2 Henry's constant for CH4.
#' @param standard.gibbs Standard Gibbs free energy for the reaction,
#' @param temperature Temperature of the system, in Kelvin.
#' @param VolumeSolution Volume of liquid in the closed system, in liters.
#' @param VolumeHeadspace Volume of gaseous headspace in the closed system, in liters.
#' @param K.CO2HCO3 Equilibrium constant for the dissociation of CO2(aq) to HCO3-(aq). 5.223196e-07 by default.
#' @param KHCO3CO3 Equilibrium constant for the dissociation of HCO3- (aq) to CO3-- (aq). 6.01886e-11 by default.
#' @return A data frame to be used for the methanogenesis function
init <- function(CH4.initial, K.CH4, H2.initial, K.H2,
                 DIC.initial, pH.initial, K.CO2, standard.gibbs,
                 temperature, VolumeSolution, VolumeHeadspace, K.CO2HCO3 = 5.223196e-07,K.HCO3CO3 = 6.01886e-11){

  # CH4.initial = 0.000001
  # K.CH4 = 0.001128969
  # H2.initial = 0.0005
  # K.H2 = 0.000666252
  # DIC.initial = 0.0032
  # pH.initial = 7.5
  # K.CO2 = 0.023464592
  # standard.gibbs = -191359.46584
  # temperature = 313.15
  # VolumeSolution = 0.08
  # VolumeHeadspace = 0.02
  # #init(0.000001, 0.001128969, 0.0005, 0.000666252, 0.0032, 7.5, 0.023464592, -45736.01, 313.15, 0.08, 0.02)


  #CH4 initials
  nCH4.solution <- CH4.initial * VolumeSolution
  PCH4.initial <- CH4.initial / K.CH4
  nCH4.headspace.initial <- (PCH4.initial * VolumeHeadspace) / (0.08206 * temperature)
  nCH4.total.initial <- nCH4.headspace.initial + nCH4.solution

  #H2 initials
  nH2.solution <- H2.initial * VolumeSolution
  PH2.initial <- H2.initial / K.H2
  nH2.headspace.initial <- (PH2.initial * VolumeHeadspace) / (0.08206 * temperature)
  nH2.total.initial <- nH2.headspace.initial + nH2.solution

  #DIC, CO2 initials
  nDIC <- DIC.initial * VolumeSolution
  alkalinity.initial <- alkalinity(pH.initial, nDIC, VolumeSolution, VolumeHeadspace, temperature)
  system.pH <- pH(nDIC, VolumeSolution, VolumeHeadspace, temperature, alkalinity.initial)
  PCO2.initial <- PCO2(system.pH, nDIC, VolumeSolution, VolumeHeadspace, temperature,K.CO2HCO3, K.HCO3CO3)
  CO2.initial <- aqueous.step(PCO2.initial, K.CO2)

  #Gibbs free energy initials
  Q.initial <- CH4.initial/((H2.initial^4)*CO2.initial)
  Gibbs.free.energy.initial <- gibbs.step(standard.gibbs,Q.initial,temperature)

  init_frame <- data.frame(CH4.initial, nCH4.solution, PCH4.initial, nCH4.headspace.initial, nCH4.total.initial,
                           H2.initial, nH2.solution, PH2.initial, nH2.headspace.initial, nH2.total.initial,
                           DIC.initial, nDIC, PCO2.initial, CO2.initial, alkalinity.initial,
                           pH.initial,Gibbs.free.energy.initial)

  return(init_frame)
}
#' Steps through DIC consumption during a methanogenesis reaction
#'
#' `methanogenesis()` calculates CH4 produced, H2 consumed, CO2 consumed, and Gibbs free energy changes as dissolved inorganic carbon is consumed.
#'
#' `init()` sets up the initial environment to be used by the methanogenesis model.
#' @param CH4.initial Concentration of initial dissolved CH4, in molarity.
#' @param K.CH4 Henry's constant for CH4.
#' @param H2.initial Concentration of initial dissolved H2, in molarity.
#' @param K.H2 Henry's constant for H2.
#' @param DIC.initial Concentration of initial dissolved inorganic carbon, in molarity.
#' @param pH.initial initial pH.
#' @param K.CO2 Henry's constant for CH4.
#' @param standard.gibbs Standard Gibbs free energy for the reaction,
#' @param temperature Temperature of the system, in Kelvin.
#' @param VolumeSolution Volume of liquid in the closed system, in liters.
#' @param VolumeHeadspace Volume of gaseous headspace in the closed system, in liters.
#' @param K.CO2HCO3 Equilibrium constant for the dissociation of CO2(aq) to HCO3-(aq). 5.223196e-07 by default.
#' @param KHCO3CO3 Equilibrium constant for the dissociation of HCO3- (aq) to CO3-- (aq). 6.01886e-11 by default.
#' @param delta.DIC step size, in millimolar. 0.1 mM by default.
#' @param biomass.yield mass of dry biomass produced per mol of product. 2.4 g/mol product.
#' @param carbon.fraction w/w percent C of biomass, expressed as a decimal. 0.44 by default.
#' @return A data frame of the model results
#' @examples
#' methanogenesis(CH4.initial = 1e-6,H2.initial = 5e-4,DIC.initial = 3.2e-3,pH.initial = 7.5,standard.gibbs = -191359.46584,temperature = 273.15+40,VolumeSolution = 80e-3,VolumeHeadspace = 20e-3,delta.DIC = 0.0001)
#'
#' @export
methanogenesis <- function(CH4.initial, K.CH4=0.00112896948941469, H2.initial, K.H2=0.000666251556662516,
                           DIC.initial, pH.initial, K.CO2=0.023464592, standard.gibbs=-191359.46584, temperature,
                           VolumeSolution, VolumeHeadspace, K.CO2HCO3 = 5.223196e-07, K.HCO3CO3 = 6.01886e-11,
                           delta.DIC=0.0001, biomass.yield=2.4,carbon.fraction=0.44){
  # CH4.initial = 0.000001
  # K.CH4 = 0.001128969
  # H2.initial = 0.0005
  # K.H2 = 0.000666252
  # DIC.initial = 0.0032
  # pH.initial = 7.5
  # K.CO2 = 0.023464592
  # standard.gibbs = -45736.01
  # temperature = 313.15
  # VolumeSolution = 0.08
  # VolumeHeadspace = 0.02
  # delta.DIC = 0.0001

  total.steps <- DIC.initial/delta.DIC #Runs until all DIC has been consumed and returns all the model output

  #make init data frame
  init <- init(CH4.initial, K.CH4, H2.initial, K.H2,
               DIC.initial, pH.initial, K.CO2, standard.gibbs,
               temperature, VolumeSolution, VolumeHeadspace)

  #make empty data frame
  columns <- c("DIC.consumed", "nDIC.consumed","CH4.produced", "nCH4.produced","H2.consumed", "nH2.consumed",
               "nCH4.total.step", "PCH4.step", "[CH4].step",  "nH2.total.step", "PH2.step", "[H2].step",
               "nDIC.total.step", "[DIC].step","PCO2.step", "[CO2].step","systempH.step","Gibbs.free.energy.step")

  main <- setNames(data.frame(matrix(ncol = length(columns), nrow = 1)), columns)

  #Set up initial conditions, index of 1

  main$DIC.consumed[1] <- 0
  main$CH4.produced[1] <- 0
  main$H2.consumed[1] <- 0

  main$nDIC.consumed[1] <- 0
  main$nCH4.produced[1] <- 0
  main$nH2.consumed[1] <- 0

  main$nDIC.total.step[1] <- init$nDIC
  main$nCH4.total.step[1] <- init$nCH4.total.initial
  main$nH2.total.step[1] <- init$nH2.total.initial

  main$PCH4.step[1] <- init$PCH4.initial
  main$PH2.step[1] <- init$PH2.initial

  main$`[DIC].step`[1] <- init$DIC.initial
  main$`[CH4].step`[1] <- init$CH4.initial
  main$`[H2].step`[1] <- init$H2.initial

  main$systempH.step[1] <- pH.initial

  main$PCO2.step[1] <- init$PCO2.initial
  main$`[CO2].step`[1] <- init$CO2.initial

  main$Gibbs.free.energy.step[1] <- init$Gibbs.free.energy.initial

  CH4.per.step <- CH4.coefficient(biomass.yield,carbon.fraction)
  H2.per.step <- H2.coefficient(biomass.yield,carbon.fraction)

  main$`[H2]/[CO2] ratio`[1] <- 0
  main$`[H2]/[DIC] ratio`[1] <- 0

  for (i in 2:total.steps){

    main <- rbind(main, NA)
    main$DIC.consumed[i] <- main$DIC.consumed[i-1]+delta.DIC
    main$CH4.produced[i] <- main$CH4.produced[i-1]+CH4.per.step*delta.DIC
    main$H2.consumed[i] <- main$H2.consumed[i-1]+H2.per.step*delta.DIC

    main$nDIC.consumed[i] <- main$DIC.consumed[i]*VolumeSolution
    main$nCH4.produced[i] <- main$CH4.produced[i]*VolumeSolution
    main$nH2.consumed[i] <- main$H2.consumed[i]*VolumeSolution

    main$nDIC.total.step[i] <- main$nDIC.total.step[1]-main$nDIC.consumed[i]
    main$nCH4.total.step[i] <- main$nCH4.total.step[1]+main$nCH4.produced[i]
    main$nH2.total.step[i] <- main$nH2.total.step[1]-main$nH2.consumed[i]

    main$PCH4.step[i] <- pressure.step(main$nCH4.total.step[i], K.CH4, VolumeSolution, VolumeHeadspace, temperature)
    main$PH2.step[i] <- pressure.step(main$nH2.total.step[i], K.H2, VolumeSolution, VolumeHeadspace, temperature)

    main$`[DIC].step`[i] <- main$nDIC.total.step[i]/VolumeSolution
    main$`[CH4].step`[i] <- aqueous.step(main$PCH4.step[i], K.CH4)
    main$`[H2].step`[i] <- aqueous.step(main$PH2.step[i], K.H2)

    if (main$`[DIC].step`[i]<0|main$`[H2].step`[i]<0){
      main <- main[1:(nrow(main)-1),]
      break

    }




    main$systempH.step[i] <- pH(main$nDIC.total.step[i], VolumeSolution, VolumeHeadspace, temperature, init$alkalinity.initial)
    main$PCO2.step[i] <- PCO2(main$systempH.step[i], main$nDIC.total.step[i], VolumeSolution, VolumeHeadspace, temperature,K.CO2HCO3, K.HCO3CO3)

    main$`[CO2].step`[i] <- aqueous.step(main$PCO2.step[i], K.CO2)

    main$`[H2]/[CO2] ratio`[i] <- (main$`[H2].step`[i]/H2.per.step)/main$`[CO2].step`[i]
    main$`[H2]/[DIC] ratio`[i] <- (main$`[H2].step`[i]/H2.per.step)/main$`[DIC].step`[i]

    Q.step <- main$`[CH4].step`[i] / (main$`[H2].step`[i]^4 * main$`[CO2].step`[i])
    main$Gibbs.free.energy.step[i] <- gibbs.step(standard.gibbs, Q.step, temperature)

  }

  return(main)
}
