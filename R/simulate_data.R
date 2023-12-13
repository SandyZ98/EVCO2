library(Rcpp)

#' Simulate CO2 time series data
#'
#' @param persondata data.frame with time (in hours since start of measurement), n (number of people), and the following optional columns: age, gender, MET, CO2rate in L/s at 1 atm & 0 Celsius. If you have Gp instead of CO2rate, call Gp_to_CO2rate first
#' @param volume volume of room in m^3
#' @param ventilation_rate ventilation rate of room in m^3/s
#' @param envCO2 environmental CO2 concentration in ppm
#' @param startCO2 starting CO2 concentration in ppm
#' @param freq frequency that data should be "measured", in seconds
#' @param CO2var variance of noise added to simulated CO2 values at each step
#' @param temp temperature of room in Celsius
#' @param method method for simulating CO2 concentration. Options are 'Euler' or 'Exponential'
#' @param useCpp whether to use Rcpp for simulation. Default is TRUE
#'
#' @return
#' @export
#'
#' @examples
#' inputdata_example <- data.frame(
#'   time = c(0, 8, 8, 17, 24),
#'   n = c(0, 15, 2, 0, 0),
#'   age = c(0, 8, 50, 0, 0),
#'   gender = rep(NA, 5),
#'   MET = c(NA, 3, 1.8, NA, NA),
#'   CO2rate = rep(NA, 5)
#' ) 
#' 
#' exdata <- simulateData(persondata = inputdata_example, volume=500, ventilation_rate = .1, CO2var = 0.01)
#' exdata2 <- simulateData(persondata = inputdata_example, volume=500, ventilation_rate = .5, CO2var = 0.01)

simulateData <- function(persondata, volume, ventilation_rate, envCO2=400, startCO2=400, freq = 1, CO2var = .01, temp = 25, method='Euler', useCpp = TRUE){
  
  ## get emissions from persondata ##
  emissions <- persondata_to_emission(persondata, temp, freq)
  times <- emissions$times
  nadj_CO2rate <- emissions$nadj_CO2rate
  
  ### Simulate CO2 concentration ###
  # Equation: dC/dt = (N*CO2rate + (envCO2-CO2)*ventilation_rate)/volume
  
  # Convert from ppm to ug/m^3
  ppm_to_ug <- function(ppm, temp){
    ppm / 44.01 * 8314 * (273.15+temp) / 1000 / 101.325
  }
  
  # convert from ug/m^3 to ppm
  ug_to_ppm <- function(ug, temp) {
    ug * 44.01 * 1000 * 101.325 / 8314 / (273.15+temp)
  }
  
  startCO2 <- ppm_to_ug(startCO2, temp)
  envCO2 <- ppm_to_ug(envCO2, temp)
  
  errors = rnorm(length(times), 0, sqrt(CO2var))
  
  # dC/dt = (N*CO2rate + (envCO2-CO2)*ventilation_rate)/volume
  # Euler method:
  # C[i+1] = C[i] + dC/dt*freq + error 
  # First-order exponential integrator (used in Batterman):
  # C[i+1] = N*CO2rate/ventilation_rate*(1-exp(-ventilation_rate/volume*freq)) + (C[i]-envCO2)*exp(-ventilation_rate/volume*freq) + envCO2 + error
  
  cppFunction('NumericVector simulateCO2_euler(double freq, double startCO2, double envCO2, double ventilation_rate, double volume, NumericVector times, NumericVector nadj_CO2rate, NumericVector errors) {
               int n = times.size();
               NumericVector CO2(n);
               CO2[0] = startCO2;
               for(int i=0; i<(n-1); ++i) {
                 CO2[i+1] = CO2[i] + (nadj_CO2rate[i] + (envCO2-CO2[i])*ventilation_rate)/volume*freq + errors[i];
               }
               return CO2;
               }')
  cppFunction('NumericVector simulateCO2_exp(double freq, double startCO2, double envCO2, double ventilation_rate, double volume, NumericVector times, NumericVector nadj_CO2rate, NumericVector errors) {
               int n = times.size();
               NumericVector CO2(n);
               CO2[0] = startCO2;
               for(int i=0; i<(n-1); ++i) {
                 CO2[i+1] = nadj_CO2rate[i]/ventilation_rate*(1-exp(-ventilation_rate/volume*freq)) + (CO2[i]-envCO2)*exp(-ventilation_rate/volume*freq) + envCO2 + errors[i];
               }
               return CO2;
               }')
  simulateCO2_euler_R <- function(freq, startCO2, envCO2, ventilation_rate, volume, times, nadj_CO2rate, errors) {
    n = length(times)
    CO2 = rep(0, n)
    CO2[1] = startCO2
    for(i in 1:(n-1)) {
      CO2[i+1] = CO2[i] + (nadj_CO2rate[i] + (envCO2-CO2[i])*ventilation_rate)/volume*freq + errors[i]
    }
    return(CO2)
  }
  simulateCO2_exp_R <- function(freq, startCO2, envCO2, ventilation_rate, volume, times, nadj_CO2rate, errors) {
    n = length(times)
    CO2 = rep(0, n)
    CO2[1] = startCO2
    for(i in 1:(n-1)) {
      CO2[i+1] = nadj_CO2rate[i]/ventilation_rate*(1-exp(-ventilation_rate/volume*freq)) + (CO2[i]-envCO2)*exp(-ventilation_rate/volume*freq) + envCO2 + errors[i]
    }
    return(CO2)
  }
  
  if(useCpp){
  if(method=="Euler"){
    CO2 <- simulateCO2_euler(freq, startCO2, envCO2, ventilation_rate, volume, times, nadj_CO2rate, errors)
  } else if(method=='Exponential') {
    CO2 <- simulateCO2_exp(freq, startCO2, envCO2, ventilation_rate, volume, times, nadj_CO2rate, errors)
  } else {
    stop("method must be 'Euler' or 'Exponential'")
  }
  } else {
    if(method=="Euler"){
      CO2 <- simulateCO2_euler_R(freq, startCO2, envCO2, ventilation_rate, volume, times, nadj_CO2rate, errors)
    } else if(method=='Exponential') {
      CO2 <- simulateCO2_exp_R(freq, startCO2, envCO2, ventilation_rate, volume, times, nadj_CO2rate, errors)
    } else {
      stop("method must be 'Euler' or 'Exponential'")
    }
  }
  
  ret = data.frame('time' = times, 'CO2' = ug_to_ppm(CO2, temp))
  return(ret)
}

