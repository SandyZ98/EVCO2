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
#'
CO2rate_to_Gp <- function(CO2rate, temp=25){
  # use 1.98 g/L as the density of CO2 at temperatures near 20-25 C
  # volume = mass / density
  # Gp L/min = CO2 ug/s * 60 s / min * 1 g / 1e6 ug  / 1.98 g / L
  return(CO2rate * 60 / 1e6 / 1.98)

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
#'
Gp_to_CO2rate <- function(Gp, temp=25){
  # use 1.98 g/L as the density of CO2 at temperatures near 20-25 C
  # mass = volume * density
  # CO2 ug/s = Gp L/min * 1 min / 60 s * 1.98 g / L * 1e6 ug / 1 g
  return(Gp / 60 * 1.98 * 1e6)
}

