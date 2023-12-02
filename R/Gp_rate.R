#' Gp_rate - Estimates Gp based on ppm of CO2 in expired breath and the respiratory rate of an individual in m^3 / hour
#' @param co2_breath : Concentration of CO2 in exhaled breath (ppm)
#' @param resp_rate : Respiratory rate of an individual in m^3 / hour. Recommended to get this value from the
#' EPA Exposure Factors handbook Chapter 6: Inhalaltion Rates table 6-23
#' @return : Returns the calculated Gp value

Gp_rate = function(resp_rate, co2_breath = 4e4){
  value = resp_rate * co2_breath / 6e4

  return(value)
}
