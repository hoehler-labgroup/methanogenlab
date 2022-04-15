#' Determines initial conditions from initial inputs
#'
#' `init()` sets up the initial environment to be used by the methanogenesis model.
#' @inheritParams methanogenesis.time
#' @return A data frame to be used for the methanogenesis function.
init <- function(CH4.initial, K.CH4, H2.initial, K.H2,
                 DIC.initial, pH.initial, K.CO2,
                 temperature, VolumeSolution, VolumeHeadspace,inoculum.cell.number,
                 biomass.yield,carbon.fraction,cell.weight,K.CO2HCO3,K.HCO3CO3){

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
  alkalinity.initial <- alkalinity(pH.initial, nDIC, VolumeSolution, VolumeHeadspace, temperature,K.CO2HCO3,K.HCO3CO3)
  system.pH <- pH(nDIC, VolumeSolution, VolumeHeadspace, temperature, alkalinity.initial,K.CO2HCO3, K.HCO3CO3)
  PCO2.initial <- PCO2(system.pH, nDIC, VolumeSolution, VolumeHeadspace, temperature,K.CO2HCO3, K.HCO3CO3)
  CO2.initial <- aqueous.step(PCO2.initial, K.CO2)

  #yield coefficients
  CH4.per.step <- CH4.coefficient(biomass.yield,carbon.fraction)
  H2.per.step <- H2.coefficient(biomass.yield,carbon.fraction)

  #ratios
  H2.CO2.ratio.initial<- (H2.initial/H2.per.step)/CO2.initial
  H2.DIC.ratio.initial <- (H2.initial/H2.per.step)/DIC.initial

  #Gibbs free energy initials
  standard.gibbs <- standard.gibbs(c("H2","CO2","CH4","H2O"),c(-4,-1,1,2),c("aq","aq","aq","l"),temperature)
  Q.initial <- CH4.initial/((H2.initial^4)*CO2.initial)
  Gibbs.free.energy.initial <- gibbs.step(standard.gibbs,Q.initial,temperature)

  #biomass
  g.biomass.per.nDIC <- biomass.coefficient(biomass.yield,carbon.fraction)*12.0107/carbon.fraction
  initial.cell.number <- inoculum.cell.number
  initial.g.biomass <- initial.cell.number*cell.weight
  initial.cell.per.mL <- initial.cell.number/(VolumeSolution*1e3)

  init_frame <- data.frame(CH4.initial, nCH4.solution, PCH4.initial, nCH4.headspace.initial, nCH4.total.initial,
                           H2.initial, nH2.solution, PH2.initial, nH2.headspace.initial, nH2.total.initial,
                           DIC.initial, nDIC, PCO2.initial, CO2.initial, alkalinity.initial,
                           pH.initial,standard.gibbs,Gibbs.free.energy.initial,CH4.per.step,H2.per.step,H2.CO2.ratio.initial,H2.DIC.ratio.initial,g.biomass.per.nDIC,
                           initial.cell.number,initial.g.biomass,initial.cell.per.mL)
  return(init_frame)
}

#' Steps through time during a methanogenesis reaction
#'
#' `methanogenesis.time()` calculates CH4 produced, H2 consumed, CO2 consumed, pH, and Gibbs free energy changes as dissolved inorganic carbon is consumed over time.
#'
#' @param CH4.initial Concentration of initial dissolved CH4, in molarity.
#' @param K.CH4 Henry's constant for CH4. NA by default (calculated by CHNOSZ).
#' @param H2.initial Concentration of initial dissolved H2, in molarity.
#' @param K.H2 Henry's constant for H2. NA by default (calculated by CHNOSZ).
#' @param is.H2.limiting Boolean. Determines whether to use the Monod equation for H2. NA by default (limiting contribution determined by the proportion of current H2 and Ks.H2).
#' @param Ks.H2 Half-saturation constant for H2 (microorganism dependent).
#' @param DIC.initial Concentration of initial dissolved inorganic carbon, in molarity.
#' @param pH.initial initial pH.
#' @param K.CO2 Henry's constant for CH4. NA by default (calculated by CHNOSZ).
#' @param is.CO2.limiting Boolean. Determines whether to use the Monod equation for CO2. NA by default (limiting contribution determined by the proportion of current CO2 and Ks.CO2).
#' @param Ks.CO2 Half-saturation constant for CO2 (microorganism dependent).
#' @param umax Maximum growth rate for the microorganism, \ifelse{html}{\out{hr<sup>-1</sup>}}{\eqn{hr^-1}}
#' @param temperature Temperature of the system, in Kelvin.
#' @param VolumeSolution Volume of liquid in the closed system, in liters.
#' @param VolumeHeadspace Volume of gaseous head space in the closed system, in liters.
#' @param K.CO2HCO3 Equilibrium constant for the dissociation of CO2(aq) to HCO3-(aq). NA by default (calculated by CHNOSZ).
#' @param K.HCO3CO3 Equilibrium constant for the dissociation of HCO3- (aq) to CO3-- (aq). NA by default (calculated by CHNOSZ).
#' @param time.step Step size in time, in hours. 0.1 hrs by default.
#' @param total.time Length of time the model will run, in hours. Note: the model may break before the total.time is reached if H2 or DIC runs out.
#' @param inoculum.cell.number Initial number of cells. 1e6 by default.
#' @param biomass.yield mass of dry biomass produced per mol of product. 2.4 g/mol product by default.
#' @param carbon.fraction w/w percent C of biomass, expressed as a decimal. 0.44 by default.
#' @param cell.weight individual cell mass, in grams. 30e-15 by default.
#' @return A data frame of the model results
#' @examples
#' methanogenesis(CH4.initial = 1e-6,H2.initial = 5e-4,DIC.initial = 3.2e-3,pH.initial = 7.5,temperature = 273.15+40,VolumeSolution = 80e-3,VolumeHeadspace = 20e-3,delta.DIC = 0.0001)
#'
#' @export
methanogenesis.time <- function(CH4.initial, K.CH4=NA, H2.initial, K.H2=NA,is.H2.limiting=NA,Ks.H2,
                           DIC.initial, pH.initial, K.CO2=NA,is.CO2.limiting=NA,Ks.CO2,umax, temperature,
                           VolumeSolution, VolumeHeadspace, K.CO2HCO3 = NA, K.HCO3CO3 = NA,
                           time.step=0.1, total.time=5, inoculum.cell.number = 1e6,biomass.yield=2.4,carbon.fraction=0.44,cell.weight=30e-15){

  #Calculates Henry's constants if they aren't already provided
  if (is.na(K.CH4)){
    K.CH4 <- calculate.KH(c("CH4","CH4"),c(-1,1),c("g","aq"),temperature = temperature,pressure = 1)
    }

  if (is.na(K.H2)){
    K.H2 <- calculate.KH(c("H2","H2"),c(-1,1),c("g","aq"),temperature = temperature,pressure = 1)
  }

  if (is.na(K.CO2)){
    K.CO2 <- calculate.KH(c("CO2","CO2"),c(-1,1),c("g","aq"),temperature = temperature,pressure = 1)
  }

  if (is.na(K.CO2HCO3)){
    K.CO2HCO3 <- calculate.KH(c("CO2","H2O","HCO3-","H+"),c(-1,-1,1,1),c("aq","l","aq","aq"), temperature = temperature,pressure = 1)
  }

  if (is.na(K.HCO3CO3)){
    K.HCO3CO3 <- calculate.KH(c("HCO3-","CO3-2","H+"),c(-1,1,1),c("aq","aq","aq"), temperature = temperature,pressure = 1)
  }


  #make init data frame
  init <- init(CH4.initial, K.CH4, H2.initial, K.H2,
               DIC.initial, pH.initial, K.CO2,
               temperature, VolumeSolution, VolumeHeadspace,inoculum.cell.number,
               biomass.yield,carbon.fraction,cell.weight,K.CO2HCO3,K.HCO3CO3)

  #make empty data frame
  columns <- c("DIC.consumed", "nDIC.consumed","CH4.produced", "nCH4.produced","H2.consumed", "nH2.consumed",
               "nCH4.total.step", "PCH4.step", "[CH4].step",  "nH2.total.step", "PH2.step", "[H2].step",
               "nDIC.total.step", "[DIC].step","PCO2.step", "[CO2].step","systempH.step","Gibbs.free.energy.step","[H2]/[CO2] ratio","[H2]/[DIC] ratio")

  main <- setNames(data.frame(matrix(ncol = length(columns), nrow = 1)), columns)

  #Set up initial conditions, index of 1
  main$time[1] <- 0
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

  standard.gibbs <- init$standard.gibbs
  main$Gibbs.free.energy.step[1] <- init$Gibbs.free.energy.initial

  CH4.per.step <- init$CH4.per.step
  H2.per.step <- init$H2.per.step

  g.biomass.per.nDIC <- init$g.biomass.per.nDIC

  main$g.biomass.step[1] <- init$initial.g.biomass
  main$cell.number.step[1] <- init$initial.cell.number
  main$cell.per.mL[1] <- init$initial.cell.per.mL

  main$`[H2]/[CO2] ratio`[1] <- init$H2.CO2.ratio.initial
  main$`[H2]/[DIC] ratio`[1] <- init$H2.DIC.ratio.initial
  percent.change.list <- c("PCH4.step","[CH4].step","PH2.step","[H2].step","[DIC].step","PCO2.step","[CO2].step","systempH.step","Gibbs.free.energy.step","cell.per.mL")

  #Calculated total number of steps in the loop
  total.steps <- total.time/time.step

  for (i in 2:total.steps){

    main <- rbind(main, NA)

    #Number of new cells sets the change in DIC
    main$cell.number.step[i] <- growth(main$cell.number.step[i-1],main$`[CO2].step`[i-1],is.CO2.limiting,Ks.CO2,main$`[H2].step`[i-1],is.H2.limiting,Ks.H2,umax,time.step)+main$cell.number.step[i-1]
    main$g.biomass.step[i] <- main$cell.number.step[i]*cell.weight
    main$cell.per.mL[i] <- main$cell.number.step[i]/(VolumeSolution*1e3)

    delta.DIC <- (main$g.biomass.step[i]-main$g.biomass.step[i-1])*carbon.fraction

    #Change in DIC determines how much H2 is consumed and CH4 produced
    main$DIC.consumed[i] <- main$DIC.consumed[i-1]+delta.DIC
    main$CH4.produced[i] <- main$CH4.produced[i-1]+CH4.per.step*delta.DIC
    main$H2.consumed[i] <- main$H2.consumed[i-1]+H2.per.step*delta.DIC

    #Bulk chemistry
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

    #Break if DIC or H2 runs oout
    if (main$`[DIC].step`[i]<0|main$`[H2].step`[i]<0){
      main <- main[1:(nrow(main)-1),]
      break

    }

    main$systempH.step[i] <- pH(main$nDIC.total.step[i], VolumeSolution, VolumeHeadspace, temperature, init$alkalinity.initial,K.CO2HCO3, K.HCO3CO3)
    main$PCO2.step[i] <- PCO2(main$systempH.step[i], main$nDIC.total.step[i], VolumeSolution, VolumeHeadspace, temperature,K.CO2HCO3, K.HCO3CO3)

    main$`[CO2].step`[i] <- aqueous.step(main$PCO2.step[i], K.CO2)

    main$`[H2]/[CO2] ratio`[i] <- (main$`[H2].step`[i]/H2.per.step)/main$`[CO2].step`[i]
    main$`[H2]/[DIC] ratio`[i] <- (main$`[H2].step`[i]/H2.per.step)/main$`[DIC].step`[i]

    Q.step <- main$`[CH4].step`[i] / (main$`[H2].step`[i]^4 * main$`[CO2].step`[i])
    main$Gibbs.free.energy.step[i] <- gibbs.step(standard.gibbs,Q.step, temperature)

    main$`time`[i] <- main$`time`[i-1]+time.step

  }


  main$`log([H2]/[CO2] ratio)` <- log10(main$`[H2]/[CO2] ratio`)
  main$`log([H2]/[DIC] ratio)` <- log10(main$`[H2]/[DIC] ratio`)


  for (column in percent.change.list){
    percent.change <- abs(((main[[column]]-main[[column]][1])/main[[column]][1])*100)
    column.name <- sprintf("percent.change from initial %s",column)
    main <- cbind(main,percent.change)
    colnames(main)[colnames(main)=="percent.change"] <- column.name
  }

  return(main)
}

#' Steps through DIC consumption during a methanogenesis reaction
#'
#' `methanogenesis.DIC()` calculates CH4 produced, H2 consumed, CO2 consumed, pH, and Gibbs free energy changes as dissolved inorganic carbon is consumed.
#'
#' @param delta.DIC step size, in millimolar DIC. 0.1 mM by default.
#' @inheritParams methanogenesis.time
#' @return A data frame of the model results
#' @examples
#' methanogenesis(CH4.initial = 1e-6,H2.initial = 5e-4,DIC.initial = 3.2e-3,pH.initial = 7.5,temperature = 273.15+40,VolumeSolution = 80e-3,VolumeHeadspace = 20e-3,delta.DIC = 0.0001)
#'
#' @export
methanogenesis.DIC <- function(CH4.initial, K.CH4=NA, H2.initial, K.H2=NA,
                           DIC.initial, pH.initial, K.CO2=NA, temperature,
                           VolumeSolution, VolumeHeadspace, K.CO2HCO3 = NA, K.HCO3CO3 = NA,
                           delta.DIC=0.0001, inoculum.cell.number = 1e6,biomass.yield=2.4,carbon.fraction=0.44,cell.weight=30e-15){

  #Calculates Henry's constants if they aren't already provided
  if (is.na(K.CH4)){
    K.CH4 <- calculate.KH(c("CH4","CH4"),c(-1,1),c("g","aq"),temperature = temperature,pressure = 1)
  }

  if (is.na(K.H2)){
    K.H2 <- calculate.KH(c("H2","H2"),c(-1,1),c("g","aq"),temperature = temperature,pressure = 1)
  }

  if (is.na(K.CO2)){
    K.CO2 <- calculate.KH(c("CO2","CO2"),c(-1,1),c("g","aq"),temperature = temperature,pressure = 1)
  }

  if (is.na(K.CO2HCO3)){
    K.CO2HCO3 <- calculate.KH(c("CO2","H2O","HCO3-","H+"),c(-1,-1,1,1),c("aq","l","aq","aq"), temperature = temperature,pressure = 1)
  }

  if (is.na(K.HCO3CO3)){
    K.HCO3CO3 <- calculate.KH(c("HCO3-","CO3-2","H+"),c(-1,1,1),c("aq","aq","aq"), temperature = temperature,pressure = 1)
  }

  total.steps <- DIC.initial/delta.DIC #Runs until all DIC has been consumed and returns all the model output

  #make init data frame
  init <- init(CH4.initial, K.CH4, H2.initial, K.H2,
               DIC.initial, pH.initial, K.CO2,
               temperature, VolumeSolution, VolumeHeadspace,inoculum.cell.number,
               biomass.yield,carbon.fraction,cell.weight,K.CO2HCO3,K.HCO3CO3)

  #make empty data frame
  columns <- c("DIC.consumed", "nDIC.consumed","CH4.produced", "nCH4.produced","H2.consumed", "nH2.consumed",
               "nCH4.total.step", "PCH4.step", "[CH4].step",  "nH2.total.step", "PH2.step", "[H2].step",
               "nDIC.total.step", "[DIC].step","PCO2.step", "[CO2].step","systempH.step","Gibbs.free.energy.step","[H2]/[CO2] ratio","[H2]/[DIC] ratio")

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

  standard.gibbs <- init$standard.gibbs
  main$Gibbs.free.energy.step[1] <- init$Gibbs.free.energy.initial

  CH4.per.step <- init$CH4.per.step
  H2.per.step <- init$H2.per.step

  g.biomass.per.nDIC <- init$g.biomass.per.nDIC

  main$g.biomass.step[1] <- init$initial.g.biomass
  main$cell.number.step[1] <- init$initial.cell.number
  main$cell.per.mL[1] <- init$initial.cell.per.mL

  main$`[H2]/[CO2] ratio`[1] <- init$H2.CO2.ratio.initial
  main$`[H2]/[DIC] ratio`[1] <- init$H2.DIC.ratio.initial
  percent.change.list <- c("PCH4.step","[CH4].step","PH2.step","[H2].step","[DIC].step","PCO2.step","[CO2].step","systempH.step","Gibbs.free.energy.step","cell.per.mL")

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

    main$systempH.step[i] <- pH(main$nDIC.total.step[i], VolumeSolution, VolumeHeadspace, temperature, init$alkalinity.initial,K.CO2HCO3, K.HCO3CO3)
    main$PCO2.step[i] <- PCO2(main$systempH.step[i], main$nDIC.total.step[i], VolumeSolution, VolumeHeadspace, temperature,K.CO2HCO3, K.HCO3CO3)

    main$`[CO2].step`[i] <- aqueous.step(main$PCO2.step[i], K.CO2)

    main$`[H2]/[CO2] ratio`[i] <- (main$`[H2].step`[i]/H2.per.step)/main$`[CO2].step`[i]
    main$`[H2]/[DIC] ratio`[i] <- (main$`[H2].step`[i]/H2.per.step)/main$`[DIC].step`[i]

    Q.step <- main$`[CH4].step`[i] / (main$`[H2].step`[i]^4 * main$`[CO2].step`[i])
    main$Gibbs.free.energy.step[i] <- gibbs.step(standard.gibbs,Q.step, temperature)

    main$g.biomass.step[i] <- main$g.biomass.step[i-1]+delta.DIC*VolumeSolution*g.biomass.per.nDIC
    main$cell.number.step[i] <- main$g.biomass.step[i]/cell.weight
    main$cell.per.mL[i] <- main$cell.number.step[i]/(VolumeSolution*1e3)

  }


  main$`log([H2]/[CO2] ratio)` <- log10(main$`[H2]/[CO2] ratio`)
  main$`log([H2]/[DIC] ratio)` <- log10(main$`[H2]/[DIC] ratio`)


  for (column in percent.change.list){
    percent.change <- abs(((main[[column]]-main[[column]][1])/main[[column]][1])*100)
    column.name <- sprintf("percent.change from initial %s",column)
    main <- cbind(main,percent.change)
    colnames(main)[colnames(main)=="percent.change"] <- column.name
  }

  return(main)
}
