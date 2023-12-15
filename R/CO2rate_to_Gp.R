# CO2rate mg/s to Gp L/min
#' Convert CO2 emission rates (in milligrams/s) to Gp (in L/min)
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
  # Gp L/min = CO2 g/s * 60 s / min * 1 g / 1e3 mg  / 1.98 g / L
  return(CO2rate * 60 / 1e3 / 1.98)

}

# Gp L/min to CO2rate mg/s
#' Convert Gp (in L/min) to CO2 emission rates (in milligrams/s)
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
  # CO2 mg/s = Gp L/min * 1 min / 60 s * 1.98 g / L * 1e3 mg / 1 g
  return(Gp / 60 * 1.98 * 1e3)
}


#' Convert from ppm to mg/m^3
#'
#' @param ppm concentration of CO2 in ppm
#' @param temp temp in celsius
#'
#' @return
#' @export
#'
#' @examples
#'
ppm_to_mg <- function(ppm, temp){
  ppm * 28.9647 / 8.314 / (273.15+temp)  * 101.325
}


#' convert from mg/m^3 to ppm
#'
#' @param mg concentration of CO2 in mg/m^3
#' @param temp temperature in celsius
#'
#' @return
#' @export
#'
#' @examples
mg_to_ppm <- function(mg, temp) {
  mg / 28.9647  / 101.325 * 8.314 * (273.15+temp)
}
