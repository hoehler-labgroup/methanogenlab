gibbs.step <- function(standard.gibbs, Q, temperature){
  Gibbs.step <- (standard.gibbs + 8.314 * temperature * log(Q)) / 1000
  return(Gibbs.step)
}
