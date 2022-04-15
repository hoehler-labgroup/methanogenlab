#' Determines the number of new cells generated for a given time step size and
#'
#' @inheritParams methanogenesis.time
#' @param cells.initial Number of cells currently in the system.
#' @param CO2 Current dissolved CO2
#' @param H2 Current dissolved H2
#' @return Number of new cells made

growth <- function(cells.initial,CO2,is.CO2.limiting=NA,Ks.CO2,H2,is.H2.limiting=NA,Ks.H2,umax,time.step){

  epsilon <- 0.9
  print(paste('H2 mu', (H2/(Ks.H2+H2))))
  print(paste('CO2 mu', (CO2/(Ks.CO2+CO2))))
  # H2 and CO2 limiting
  if (is.na(is.CO2.limiting) && is.na(is.H2.limiting)){
    if((H2/(Ks.H2+H2)) <= epsilon && CO2/(Ks.CO2+CO2) <= epsilon ){
      cells.new <- (cells.initial*umax*((CO2/(Ks.CO2+CO2))*(H2/(Ks.H2+H2))))*time.step
      print('both limited')
    }
    #H2 limiting
    else if ((H2/(Ks.H2+H2)) <= epsilon){
      cells.new <- (cells.initial*umax*(H2/(Ks.H2+H2)))*time.step
      print('h2 limited')
    }
    #CO2 limiting
    else{
      cells.new <- (cells.initial*umax*(CO2/(Ks.CO2+CO2)))*time.step
      print('co2 limited')
    }
  }
  else {
    if(is.CO2.limiting==TRUE && is.H2.limiting==TRUE){
      cells.new <- (cells.initial*umax*((CO2/(Ks.CO2+CO2))*(H2/(Ks.H2+H2))))*time.step
      print('both limited')
    }
    #H2 limiting
    else if (is.H2.limiting==TRUE){
      cells.new <- (cells.initial*umax*(H2/(Ks.H2+H2)))*time.step
      print('h2 limited')
    }
    #CO2 limiting
    else{
      cells.new <- (cells.initial*umax*(CO2/(Ks.CO2+CO2)))*time.step
      print('co2 limited')
    }
  }
  return(cells.new)
}
