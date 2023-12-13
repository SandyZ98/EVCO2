# CO2rate adjusted to room temperature and convert to ug/s
# Gmole = 0.0446*G (L/s) (/60 if L/min)
# Gadj (L/s) = 8.314 (J/mol/K) * (273.15+temp) (K) / 101.325 (kPa) * Gmole (mol/s)
# note: kPa = J/L
# Gadj (g/s) = 1.965 * Gadj (L/s)
# putting it all together: convert from L/s to ug/s at a particular temperature

# CO2rate ug/s to Gp L/min
#' Convert CO2 emission rates (in micrograms/s) to Gp (in L/min)
#'
#' @param CO2rate 
#' @param temp 
#'
#' @return Gp
#' @export
#'
#' @examples
CO2rate_to_Gp <- function(CO2rate, temp=25){
  return(CO2rate / 1.965 / .0046 * 60 / (8.314*(273.15+temp)/101.325) / 1000)
  
}

# Gp L/min to CO2rate ug/s
#' Convert Gp (in L/min) to CO2 emission rates (in micrograms/s)
#'
#' @param Gp 
#' @param temp 
#'
#' @return CO2rate
#' @export
#'
#' @examples
Gp_to_CO2rate <- function(Gp, temp=25){
  return(Gp * 1.965 * .0046 / 60* (8.314*(273.15+temp)/101.325) * 1000)
}

