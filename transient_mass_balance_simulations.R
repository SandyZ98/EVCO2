library(EVCO2)
library(dplyr)
library(ggplot2)
library(tidyr)
library(bench)
library(patchwork)

## Simulation performance ##
# benchmarking data simulation
persondata <- data.frame(
  time = c(0, 8, 8, 12, 12, 13, 13, 17, 24),
  n = c(0, 15, 2, 0, 0, 15, 2, 0, 0),
  age = c(0, 8, 30, 8, 30, 8, 30, 0, 0),
  gender = rep(NA, 9),
  MET = c(NA, 3, 1.8, 3, 1.8, 3, 2, NA, NA),
  CO2rate = rep(NA, 9)
)

freq.sim=1
freq.est = 120
volume=500
envCO2=400
temp=25

Cppresults <- bench::mark(simulate_data(persondata = persondata,
                                       volume=volume,
                                       ventilation_rate = .1,
                                       CO2var = 1,
                                       method='Exponential',
                                       freq=freq.sim), iterations = 100)
Rresults <- bench::mark(simulate_data(persondata = persondata,
                                     volume=volume,
                                     ventilation_rate = .1,
                                     CO2var = 1,
                                     method='Exponential',
                                     freq=freq.sim, useCpp = FALSE), iterations = 100)
Rtimes <- as.numeric(Rresults$time[[1]])
Cpptimes <- as.numeric(Cppresults$time[[1]])
mean(Rtimes)
sd(Rtimes)
mean(Cpptimes)
sd(Cpptimes)

# example of what the simulated data looks like
sim_example <- simulate_data(persondata = persondata,
                            volume=volume,
                            ventilation_rate = .1,
                            CO2var = .01,
                            method='Exponential',
                            freq=freq.sim)
plot.01 <- ggplot(sim_example, aes(x=time, y=CO2)) + geom_line() + labs(x="Time (hours)", y="CO2 (ppm)") + theme_bw() + ggtitle('CO2 variance = 0.01') + geom_vline(xintercept=c(8, 12, 13, 17), col='red', lty=2)
sim_example <- simulate_data(persondata = persondata,
                            volume=volume,
                            ventilation_rate = .1,
                            CO2var = .1,
                            method='Exponential',
                            freq=freq.sim)
plot.1 <- ggplot(sim_example, aes(x=time, y=CO2)) + geom_line() + labs(x="Time (hours)", y="CO2 (ppm)") + theme_bw() + ggtitle('CO2 variance = 0.1') + geom_vline(xintercept=c(8, 12, 13, 17), col='red', lty=2)
sim_example <- simulate_data(persondata = persondata,
                            volume=volume,
                            ventilation_rate = .1,
                            CO2var = 1,
                            method='Exponential',
                            freq=freq.sim)
plot1 = ggplot(sim_example, aes(x=time, y=CO2)) + geom_line() + labs(x="Time (hours)", y="CO2 (ppm)") + theme_bw() + ggtitle('CO2 variance = 1') + geom_vline(xintercept=c(8, 12, 13, 17), col='red', lty=2)
sim_example <- simulate_data(persondata = persondata,
                            volume=volume,
                            ventilation_rate = .1,
                            CO2var = 2,
                            method='Exponential',
                            freq=freq.sim)
plot10 = ggplot(sim_example, aes(x=time, y=CO2)) + geom_line() + labs(x="Time (hours)", y="CO2 (ppm)") + theme_bw() + ggtitle('CO2 variance = 2') + geom_vline(xintercept=c(8, 12, 13, 17), col='red', lty=2)
plot.01+plot.1+plot1+plot10


## 1-D Newton Raphson Accuracy ##
persondata <- data.frame(
  time = c(0, 8, 8, 12, 12, 13, 13, 17, 24),
  n = c(0, 15, 2, 0, 0, 15, 2, 0, 0),
  age = c(0, 8, 30, 8, 30, 8, 30, 0, 0),
  gender = rep(NA, 9),
  MET = c(NA, 3, 1.8, 3, 1.8, 3, 2, NA, NA),
  CO2rate = rep(NA, 9)
)

freq.sim=1
freq.est = 120
volume=500
envCO2=400
temp=25

# this function simulates a dataset then runs 1-D Newton Raphson on it
testing_NR <- function(Q, var, init.Q){
  sim <- simulate_data(persondata = persondata,
                      volume=volume,
                      ventilation_rate = Q,
                      CO2var = var,
                      method='Exponential',
                      freq=freq.sim)
  index = seq(1, length(sim$CO2), freq.est)
  CO2 = sim$CO2[index]
  out <- transient_mass_balance(freq=freq.est,
                                CO2=CO2,
                                volume=volume,
                                envCO2known=envCO2,
                                init.Q=init.Q,
                                temp=temp,
                                persondata=persondata,
                                method='NR',
                                max.iter=1000,
                                tol=1e-10)
  return(c(Q.est = out$ventilation, iter = out$iter, convergence = out$convergence, Q.true = Q, var = var, Q.init = init.Q))
}

Q = rep(seq(0.01, 1, length=10), 5) # do 5 simulations on each Q
var = c(.01, .1, 1, 2)
test_cases = merge(Q, var)


out_results = matrix(NA, nrow = nrow(test_cases), ncol = 6)
colnames(out_results) = c("Q.est", "iter", "convergence", "Q.true", "var", "Q.init")
for(i in 1:nrow(test_cases)){
  #print(paste("Running test case", i, "of", nrow(test_cases), "test cases."))
  out_results[i,] = testing_NR(test_cases[i, 1], test_cases[i, 2], 1)
}

out_results <- as.data.frame(out_results)
out_results <- out_results %>% mutate(rel.error = abs(Q.est-Q.true)/Q.true)
summary <- out_results %>% group_by(Q.true, var) %>% summarize(mean.Q.est = mean(Q.est), mean.iter = mean(iter), mean.convergence = mean(convergence), mean.rel.error = mean(rel.error), sd.rel.error = sd(rel.error))

ggplot(summary, aes(x=Q.true, y=mean.rel.error)) +     geom_point() +
  geom_line() +
  geom_ribbon(aes(ymin=mean.rel.error-sd.rel.error, ymax=mean.rel.error+sd.rel.error), alpha=.2) +
  facet_wrap(~paste('Simulated CO2 variance:',var))+
  labs(x="True Ventilation Rate (m3/s)", y="Relative error in estimated ventilation rate") +
  theme_bw() +
  ylim(-0.01, .125)


## Multi-dimensional Accuracy ##
persondata <- data.frame(
  time = c(0, 8, 8, 17, 24),
  n = c(0, 15, 2,  0, 0),
  age = c(0, 8, 30,  0, 0),
  gender = rep(NA, 5),
  MET = c(NA, 3, 1.8,  NA, NA),
  CO2rate = rep(NA, 5)
)
freq.sim=1
freq.est = 120
volume=500
envCO2=400
temp=25

testing_all <- function(Q, envCO2, var, init.Q, init.envCO2){
  sim <- simulate_data(persondata = persondata,
                      volume=volume,
                      ventilation_rate = Q,
                      CO2var = var,
                      method='Exponential',
                      freq=freq.sim,
                      envCO2=envCO2)
  index = seq(1, length(sim$CO2), freq.est)
  CO2 = sim$CO2[index]
  outNewton <- tryCatch({
    transient_mass_balance(freq=freq.est,
                           CO2=CO2,
                           volume=volume,
                           envCO2.init=init.envCO2,
                           init.Q=init.Q,
                           temp=temp,
                           persondata=persondata,
                           method='Newton',
                           max.iter=1000,
                           tol=1e-6)
  },
  error= function(e) return(list('ventilation' = NA, 'envCO2' = NA, 'iter' = NA, 'convergence' = NA)))
  outBFGS <- transient_mass_balance(freq=freq.est,
                                    CO2=CO2,
                                    volume=volume,
                                    envCO2.init=init.envCO2,
                                    init.Q=init.Q,
                                    temp=temp,
                                    persondata=persondata,
                                    method='L-BFGS-B',
                                    max.iter=1000,
                                    tol=1e-6)
  outNM <- transient_mass_balance(freq=freq.est,
                                  CO2=CO2,
                                  volume=volume,
                                  envCO2.init=init.envCO2,
                                  init.Q=init.Q,
                                  temp=temp,
                                  persondata=persondata,
                                  method='Nelder-Mead',
                                  max.iter=1000,
                                  tol=1e-6)
  return(c(Newton.Q.est = outNewton$ventilation,
           BFGS.Q.est = outBFGS$ventilation,
           NM.Q.est = outNM$ventilation,
           Newton.envCO2.est = outNewton$envCO2,
           BFGS.envCO2.est = outBFGS$envCO2,
           NM.envCO2.est = outNM$envCO2,
           Newton.iter = outNewton$iter,
           BFGS.iter = outBFGS$iter,
           NM.iter = outNM$iter,
           Newton.convergence = outNewton$convergence,
           BFGS.convergence = outBFGS$convergence,
           NM.convergence = outNM$convergence,
           Q.true = Q,
           envCO2.true = envCO2,
           var = var,
           Q.init = init.Q,
           envCO2.init = init.envCO2))
}

Q = rep(seq(0.01, .5, length=5))
envCO2 = rep(c(385, 400, 415))
Q.init = seq(.1, 1, length=10)
envCO2.init = c(375, 400, 425)
var = .01
test_cases = matrix(NA, nrow = 0, ncol = 4)
for(i in 1:length(Q)){
  for(j in 1:length(envCO2)){
    for(k in 1:length(Q.init)){
      for(l in 1:length(envCO2.init)){
        params <- c(Q[i], envCO2[j], Q.init[k], envCO2.init[l])
        test_cases <- rbind(test_cases, params)
      }
    }
  }
}


out_results = matrix(NA, nrow = nrow(test_cases), ncol = 17)
colnames(out_results) = c('Newton.Q.est', 'BFGS.Q.est', 'NM.Q.est', 'Newton.envCO2.est', 'BFGS.envCO2.est', 'NM.envCO2.est', 'Newton.iter', 'BFGS.iter', 'NM.iter', 'Newton.convergence', 'BFGS.convergence', 'NM.convergence', 'Q.true', 'envCO2.true', 'var', 'Q.init', 'envCO2.init')
for(i in 1:nrow(test_cases)){
  #print(paste("Running test case", i, "of", nrow(test_cases), "test cases."))
  out_results[i,] = testing_all(test_cases[i, 1], test_cases[i, 2], var=.01, init.Q=test_cases[i, 3], init.envCO2=test_cases[i, 4])
  #print(testing_Newton(test_cases[i, 1], test_cases[i, 2], var=.01, init.Q=.5, init.envCO2=425))
}

out_results <- as.data.frame(out_results)
# pivot longer so that there is a column for method (Newton, BFGS, or NM), and columns for Q.est, envCO2.est, convergence, and iter
out_results <- out_results %>% pivot_longer(cols = c(Newton.Q.est, BFGS.Q.est, NM.Q.est,
                                      Newton.envCO2.est, BFGS.envCO2.est, NM.envCO2.est,
                                      Newton.convergence, BFGS.convergence, NM.convergence,
                                      Newton.iter, BFGS.iter, NM.iter
                                      ), names_to = c("method", ".value"), names_pattern = "(\\w+)\\.(\\w+)")
out_results <- out_results %>% mutate(rel.error.Q = abs(Q-Q.true)/Q.true) %>% mutate(rel.error.env = abs(envCO2-envCO2.true)/envCO2.true)


summary <- out_results %>% group_by(Q.true, envCO2.true, method) %>%
  summarize(mean.Q = mean(na.omit(Q)), mean.iter = mean(na.omit(iter)),
            mean.convergence = mean(na.omit(convergence)), mean.rel.error.Q = mean(na.omit(rel.error.Q)),
            sd.rel.error.Q = sd(na.omit(rel.error.Q)), mean.rel.error.env = mean(na.omit(rel.error.env)),
            sd.rel.error.env = sd(na.omit(rel.error.env)),NAs = sum(is.na(iter)),
            sd.iter = sd(na.omit(iter)), sd.convergence = sd(na.omit(convergence)))


summary <- summary %>% pivot_longer(cols = c(mean.rel.error.Q, mean.rel.error.env), names_to = "parameter", values_to = "mean.rel.error")
summary <- summary %>% pivot_longer(cols = c(sd.rel.error.Q, sd.rel.error.env), names_to = "parameter2", values_to = "sd.rel.error")
summary <- summary %>% mutate(parameter = factor(parameter, levels = c("mean.rel.error.Q", "mean.rel.error.env"), labels=c("Ventilation rate", "Environmental CO2")))
summary <- summary %>% mutate(parameter2 = factor(parameter2, levels = c("sd.rel.error.Q", "sd.rel.error.env"), labels=c("Ventilation rate", "Environmental CO2")))
summary <- summary %>% filter(parameter == parameter2)

summary$method <- factor(summary$method, levels=c("Newton", "BFGS", "NM"), labels=c("Newton", "L-BFGS-B", "Nelder Mead"))

ggplot(summary, aes(x=Q.true, y=mean.rel.error, color=parameter))+     geom_point() +
  geom_line() +
  geom_ribbon(aes(ymin=mean.rel.error-sd.rel.error, ymax=mean.rel.error+sd.rel.error), alpha=.2) +
  facet_grid(rows= vars(method), cols=vars(paste('Simulated environmental CO2:', envCO2.true)))+
  labs(x="True Ventilation Rate (m3/s)", y="Relative error in estimated parameters") + theme_bw()  + labs(color="Parameter")

ggplot(summary, aes(x=Q.true, y=mean.iter, color = method))+     geom_point() +
  geom_line() +
  geom_ribbon(aes(ymin=mean.iter-sd.iter, ymax=mean.iter+sd.iter), alpha=.2) +
  facet_wrap(vars(paste('Simulated environmental CO2:', envCO2.true)))+
  labs(x="True Ventilation Rate (m3/s)", y="Number of iterations until convergence") + theme_bw()  + labs(color="Method")

summary %>% group_by(Q.true, envCO2.true, method) %>% select(NAs) %>% filter(NAs > 0) %>% print(n=30)
sum(summary$NAs)/2



### Real data
times = example_data$`Time(dd/mm/yyyy)`
times = as.POSIXct(times, format = "%d/%m/%Y %I:%M:%S %p")
CO2 = example_data$`Carbon dioxide(ppm)`
ind0 <- which(times == "2023-11-13 18:01:35 EST") # people leave for the day
ind1 <- which(times == "2023-11-14 8:01:36 EST") # people come in
ind2 <- which(times == "2023-11-14 18:01:36 EST") # people leave
ind3 <- which(times == "2023-11-15 8:01:35 EST") # people come in
ind4 <- which(times == "2023-11-15 18:01:35 EST") # people leave

hr0 <- as.numeric(difftime(times[ind0], times[1], units = "hours"))
hr1 <- as.numeric(difftime(times[ind1], times[1], units = "hours"))
hr2 <- as.numeric(difftime(times[ind2], times[1], units = "hours"))
hr3 <- as.numeric(difftime(times[ind3], times[1], units = "hours"))
hr4 <- as.numeric(difftime(times[ind4], times[1], units = "hours"))
end <- as.numeric(difftime(times[length(times)], times[1], units = "hours"))
Gp = Gp_rate(16/24)
# transient mass balance method uses CO2rate in L/s at 0C, so we need to convert:
# first convert L to mol, then use ideal gas law to convert to 0C
# divide by 60 to convert from L/min to L/s
CO2rate = Gp * .0446 * (8.314*273.15/101.325) / 60

# create persondata dataframe
ex_persondata <- data.frame(
  time = c(0, hr0, hr1, hr2, hr3, hr4, end+.02), # include small offset at end
  n = c(5, 0, 5, 0, 5, 0, 0),
  CO2rate = c(CO2rate, 0, CO2rate, 0, CO2rate, 0, 0)
)
# volume is based on measurements of the room and should be in m^3
volume = .0254^3*321*378*114
# environmental CO2 concentration in ppm
envCO2 = 400
# room temperature in degrees Celsius
temp = 25
# freq is time step between measurements
# note, measurements are equally spaced
freq = as.numeric(difftime(times[2], times[1], units = "secs"))
NRreal <- transient_mass_balance(freq=freq,
                       CO2=CO2,
                       envCO2known=envCO2,
                       volume=volume,
                       init.Q = 1,
                       temp=temp,
                       persondata=ex_persondata,
                       method='NR')
NRrealach <- NRreal$ACH
NRrealf <- NRreal$f
NRrealiter <- NRreal$iter
# try using Newton method
Newtonreal <- transient_mass_balance(freq=freq,
                       CO2=CO2,
                       envCO2.init=400,
                       volume=volume,
                       init.Q = 3,
                       temp=temp,
                       persondata=ex_persondata,
                       method='Newton',
                       tol=1e-6)
Newtonrealach <- Newtonreal$ACH
NewtonrealCO2 <- Newtonreal$envCO2
Newtonrealf <- Newtonreal$f
Newtonrealiter <- Newtonreal$iter
# Try using the L-BFGS-B method
BFGSreal <- transient_mass_balance(freq=freq,
                       CO2=CO2,
                       envCO2known=NULL,
                       envCO2.init=400,
                       volume=volume,
                       init.Q = 1,
                       temp=temp,
                       persondata=ex_persondata,
                       method='L-BFGS-B')
BFGSrealach <- BFGSreal$ACH
BFGSrealCO2 <- BFGSreal$envCO2
BFGSrealf <- BFGSreal$f
BFGSrealiter <- BFGSreal$iter

NMreal <- transient_mass_balance(freq=freq,
                       CO2=CO2,
                       envCO2known=NULL,
                       envCO2.init=400,
                       volume=volume,
                       init.Q = 1,
                       temp=temp,
                       persondata=ex_persondata,
                       method='Nelder-Mead')
NMrealach <- NMreal$ACH
NMrealCO2 <- NMreal$envCO2
NMrealf <- NMreal$f
NMrealiter <- NMreal$iter

BFGSreal_estE <- transient_mass_balance(freq=freq,
                       CO2=CO2,
                       envCO2.init=400,
                       volume=volume,
                       init.Q = 1,
                       temp=temp,
                       ELB = 0,
                       EUB = 100,
                       critpoints=c(ind0, ind1, ind2, ind3, ind4),
                       method='L-BFGS-B')
BFGSreal_estEach <- BFGSreal_estE$ACH
BFGSreal_estECO2 <- BFGSreal_estE$envCO2
BFGSreal_estEf <- BFGSreal_estE$f
BFGSreal_estEiter <- BFGSreal_estE$iter

# Try using the Nelder-Mead method
NMreal_estE <- transient_mass_balance(freq=freq,
                       CO2=CO2,
                       envCO2.init=400,
                       volume=volume,
                       init.Q = 1,
                       temp=temp,
                       critpoints=c(ind0, ind1, ind2, ind3, ind4),
                       method='Nelder-Mead')
NMreal_estEach <- NMreal_estE$ACH
NMreal_estECO2 <- NMreal_estE$envCO2
NMreal_estEf <- NMreal_estE$f
NMreal_estEiter <- NMreal_estE$iter

LBFGSreal_estEfreq <- transient_mass_balance(freq=freq,
                       CO2=CO2,
                       envCO2.init=400,
                       volume=volume,
                       init.Q = 1,
                       temp=temp,
                       critpoints=c(seq(2, length(CO2)-1), by=3),
                       method='L-BFGS-B')
LBFGSreal_estEfreqach <- LBFGSreal_estEfreq$ACH
LBFGSreal_estEfreqCO2 <- LBFGSreal_estEfreq$envCO2
LBFGSreal_estEfreqf <- LBFGSreal_estEfreq$f
LBFGSreal_estEfreqiter <- LBFGSreal_estEfreq$iter

NMreal_estEfreq <- transient_mass_balance(freq=freq,
                       CO2=CO2,
                       envCO2.init=400,
                       volume=volume,
                       init.Q = 1,
                       temp=temp,
                       critpoints=c(seq(2, length(CO2)-1), by=3),
                       method='Nelder-Mead')
NMreal_estEfreqach <- NMreal_estEfreq$ACH
NMreal_estEfreqCO2 <- NMreal_estEfreq$envCO2
NMreal_estEfreqf <- NMreal_estEfreq$f
NMreal_estEfreqiter <- NMreal_estEfreq$iter

# make results into table
realdata <- data.frame(
  method = c("Newton Raphson", "Newton", "L-BFGS-B", "Nelder-Mead", "L-BFGS-B", "Nelder-Mead", "Nelder-Mead"),
  estE = c('No', 'No', 'No', 'No', 'At 9am & 5pm', 'At 9 am & 5 pm', 'Every 6 minutes'),
  estCO2 = c('No', 'Yes', 'Yes', 'Yes', 'Yes', 'Yes', 'Yes'),
  ACH = c(NRrealach, Newtonrealach, BFGSrealach, NMrealach, BFGSreal_estEach, NMreal_estEach, NMreal_estEfreqach),
  envCO2 = c(400, NewtonrealCO2, BFGSrealCO2, NMrealCO2,  BFGSreal_estECO2, NMreal_estECO2, NMreal_estEfreqCO2),
  SumofSquares = c(NRrealf, Newtonrealf, BFGSrealf, NMrealf, BFGSreal_estEf, NMreal_estEf, NMreal_estEfreqf),
  iterations = c(NRrealiter, Newtonrealiter, BFGSrealiter, NMrealiter, BFGSreal_estEiter, NMreal_estEiter, NMreal_estEfreqiter)
)

library(knitr)
kable(realdata, col.names = c('Method', 'Estimate Emission Rates', 'Estimate CO2', 'Ventilation (ACH)', 'envCO2 (ppm)', 'Sum of Squares', 'Iterations'))
