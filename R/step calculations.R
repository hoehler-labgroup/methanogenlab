pressure.step <- function(n.total.step, K, VolumeHeadspace, temperature, VolumeSolution){ #calculates pressure of headspace in gas
  P.step <- n.total.step / ((VolumeHeadspace / (0.08206 * temperature)) + (K * VolumeSolution))
  return(P.step)
}

aqueous.step <- function(P.step, K){ #calculates [gas]
  conc.step <- P.step * K
  return(conc.step)
}

gibbs.step <- function(standard.gibbs, Q, temperature){
  Gibbs.step <- (standard.gibbs + 8.314 * temperature * log(Q)) / 1000
  return(Gibbs.step)
}
