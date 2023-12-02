#' build_up() build up method for calculating Air Exchange Rate. This function is the objective function whose
#' roots we want to find with respect to the Air Exchange Rate
#' @param aer : Air Exchange Rate ( 1 / hour)
#' @param n : number of individuals
#' @param Gp : Average age-adjusted CO2 generation rate (L / (min·person))
#' @param V : Room volume in meters cubed
#' @param C0 : Initial CO2 concentration (ppm)
#' @param C1 : Final CO2 concentration (ppm)
#' @param delta_t : Time difference (hours)
#' @param Cr : Replacement air CO2 concentration. Will set to a default value of ambient air (400 ppm)
#' @returns : Return the value of the equation given the input parameters
build_up = function(aer, n, Gp, V, C0, C1, delta_t, Cr = 400){
  temp = (6e4 * n * Gp / (V * aer)) + Cr

  value = exp(aer * delta_t) - ( (temp - C0) / (temp - C1) )

  return(value)

}
