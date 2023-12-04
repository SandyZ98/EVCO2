library(Rcpp)


# frequency in seconds
# CO2 in ppm 
# temp in celsius
# volume in m^3
# ventilation rate in m^3/s 
# input data has following colums:
# time, n, optional: age, gender, MET, CO2rate in L/s at 1 atm & 0 Celsius 
# method choices: Euler's method or Exponential integrator 
simulateData <- function(persondata, volume, ventilation_rate, envCO2=400, startCO2=400, freq = 1, CO2var = 1, temp = 25, method='Euler', useCpp = TRUE){
  
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

# ## Examples
# 
# inputdata_example <- data.frame(
#   time = c(0, 8, 8, 17, 24),
#   n = c(0, 15, 2, 0, 0),
#   age = c(0, 8, 50, 0, 0),
#   gender = rep(NA, 5),
#   MET = c(NA, 3, 1.8, NA, NA),
#   CO2rate = rep(NA, 5)
# ) 
# 
# exdata <- simulateData(persondata = inputdata_example, volume=500, ventilation_rate = .1, CO2var = 0.01)
# exdata2 <- simulateData(persondata = inputdata_example, volume=500, ventilation_rate = .5, CO2var = 0.01)
# 
# library(ggplot2)
# library(patchwork)
# plot1 <- ggplot(exdata, aes(x=time, y=CO2)) + geom_point(shape='.') + ggtitle(label='Volume=500 m3, Ventilation rate = .1 m3/s', subtitle='17 people from 8 am to 5 pm') + ylab('CO2 (ppm)') + xlab('Time (hrs)') + ylim(300, 1000)
# plot2 <- ggplot(exdata2, aes(x=time, y=CO2)) + geom_point(shape='.') + ggtitle(label='Volume=500 m3, Ventilation rate = .5 m3/s', subtitle='17 people from 8 am to 5 pm') + ylab('CO2 (ppm)') + xlab('Time (hrs)') + ylim(300, 1000)
# plot1 + plot2 + plot_annotation(title='Euler method')
# 
# # compare to using Exponential method
# exdata_exp <- simulateData(persondata = inputdata_example, volume=500, ventilation_rate = .1, CO2var = 0.01, method='Exponential')
# exdata2_exp <- simulateData(persondata = inputdata_example, volume=500, ventilation_rate = .5, CO2var = 0.01, method='Exponential')
# plot1_exp <- ggplot(exdata_exp, aes(x=time, y=CO2)) + geom_point(shape='.') + ggtitle(label='Volume=500 m3, Ventilation rate = .1 m3/s', subtitle='17 people from 8 am to 5 pm') + ylab('CO2 (ppm)') + xlab('Time (hrs)') + ylim(300, 1000)
# plot2_exp <- ggplot(exdata2_exp, aes(x=time, y=CO2)) + geom_point(shape='.') + ggtitle(label='Volume=500 m3, Ventilation rate = .5 m3/s', subtitle='17 people from 8 am to 5 pm') + ylab('CO2 (ppm)') + xlab('Time (hrs)') + ylim(300, 1000)
# plot1_exp + plot2_exp + plot_annotation(title='Exponential method')
# 
# 
# inputdata_example_staggered_entrance <- data.frame(
#   time = c(0, 8, 8, 10,10, 12,12, 17, 24),
#   n = c(0, 5, 2, 10, 2, 15, 2, 0, 0),
#   age = c(0, 8, 30, 8, 30, 8, 30, 0, 0),
#   gender = rep(NA, 9),
#   MET = c(NA, 3, 1.8, 3, 1.8, 3, 1.8, NA, NA),
#   CO2rate = rep(NA, 9)
# ) 
# exdata3 <- simulateData(persondata = inputdata_example_staggered_entrance, volume=500, ventilation_rate = .1, CO2var = 0.01)
# exdata4 <- simulateData(persondata = inputdata_example_staggered_entrance, volume=500, ventilation_rate = .5, CO2var = 0.01)
# 
# plot3 <- ggplot(exdata3, aes(x=time, y=CO2)) + geom_point(shape='.') + ggtitle(label='Volume=500 m3, Ventilation rate = .1 m3/s', subtitle='17 people; staggered entrances') + ylab('CO2 (ppm)') + xlab('Time (hrs)') + ylim(300, 1000)
# plot4 <- ggplot(exdata4, aes(x=time, y=CO2)) + geom_point(shape='.') + ggtitle(label='Volume=500 m3, Ventilation rate = .5 m3/s', subtitle='17 people; staggered entrances') + ylab('CO2 (ppm)') + xlab('Time (hrs)') + ylim(300, 1000)
# plot3 + plot4
# 
# 
# inputdata_example_moving <- data.frame(
#   time = c(0, 8, 8, 12, 12, 13, 13, 17, 24),
#   n = c(0, 15, 2, 0, 0, 15, 2, 0, 0),
#   age = c(0, 8, 30, 8, 30, 8, 30, 0, 0),
#   gender = rep(NA, 9),
#   MET = c(NA, 3, 1.8, 3, 1.8, 3, 2, NA, NA),
#   CO2rate = rep(NA, 9)
# ) 
# exdata5 <- simulateData(persondata = inputdata_example_moving, volume=500, ventilation_rate = .1, CO2var = 0.01)
# ggplot(exdata5, aes(x=time, y=CO2)) + geom_point(shape='.') + ggtitle(label='Volume=500 m3, Ventilation rate = .1 m3/s', subtitle='17 people; go in and out') + ylab('CO2 (ppm)') + xlab('Time (hrs)') + ylim(300, 1000)
