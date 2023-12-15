## Simulation performance ##
# benchmarking data simulation
library(bench)
Cppresults <- bench::mark(simulateData(persondata = persondata,
                                       volume=volume,
                                       ventilation_rate = .1,
                                       CO2var = 1,
                                       method='Exponential',
                                       freq=freq), iterations = 100)
Rresults <- bench::mark(simulateData(persondata = persondata,
                                     volume=volume,
                                     ventilation_rate = .1,
                                     CO2var = 1,
                                     method='Exponential',
                                     freq=freq, useCpp = FALSE), iterations = 100)
Rtimes <- as.numeric(Rresults$time[[1]])
Cpptimes <- as.numeric(Cppresults$time[[1]])
mean(Rtimes)
sd(Rtimes)
mean(Cpptimes)
sd(Cpptimes)


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
  sim <- simulateData(persondata = persondata,
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
  return(c(Q.est = out$est_ventilation, iter = out$iter, convergence = out$convergence, Q.true = Q, var = var, Q.init = init.Q))
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
  ylim(-0.1, 1)


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
  sim <- simulateData(persondata = persondata,
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
  return(c(Q.est.Newton = outNewton$ventilation,
           Q.est.BFGS = outBFGS$ventilation,
           Q.est.NM = outNM$ventilation,
           envCO2.est.Newton = outNewton$envCO2,
           envCO2.est.BFGS = outBFGS$envCO2,
           envCO2.est.NM = outNM$envCO2,
           iter.Newton = outNewton$iter,
           iter.BFGS = outBFGS$iter,
           iter.NM = outNM$iter,
           convergence.Newton = outNewton$convergence,
           convergence.BFGS = outBFGS$convergence,
           convergence.NM = outNM$convergence,
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
colnames(out_results) = c('Q.est.Newton', 'Q.est.BFGS', 'Q.est.NM',
                          'envCO2.est.Newton', 'envCO2.est.BFGS', 'envCO2.est.NM',
                          'iter.Newton', 'iter.BFGS', 'iter.NM',
                          'convergence.Newton', 'convergence.BFGS', 'convergence.NM',
                          'Q.true', 'envCO2.true', 'var', 'Q.init', 'envCO2.init')
for(i in 1:nrow(test_cases)){
  #print(paste("Running test case", i, "of", nrow(test_cases), "test cases."))
  out_results[i,] = testing_all(test_cases[i, 1], test_cases[i, 2], var=.01, init.Q=test_cases[i, 3], init.envCO2=test_cases[i, 4])
  #print(testing_Newton(test_cases[i, 1], test_cases[i, 2], var=.01, init.Q=.5, init.envCO2=425))
}

out_results <- as.data.frame(out_results)
# next line is not working quite right
out_results <- out_results %>% pivot_longer(cols = c(Q.est.Newton, Q.est.BFGS, Q.est.NM, envCO2.est.Newton, envCO2.est.BFGS, envCO2.est.NM, iter.Newton, iter.BFGS, iter.NM, convergence.Newton, convergence.BFGS, convergence.NM), names_to = c('method', '.value'), names_pattern = '(.*)\\.(.*)')
out_results <- out_results %>% mutate(rel.error.Q = abs(Q.est-Q.true)/Q.true) %>% mutate(rel.error.env = abs(envCO2.est-envCO2.true)/envCO2.true)


summary <- out_results %>% group_by(Q.true, envCO2.true) %>% summarize(mean.Q.est = mean(na.omit(Q.est)), mean.iter = mean(na.omit(iter)), mean.convergence = mean(na.omit(convergence)), mean.rel.error.Q = mean(na.omit(rel.error.Q)), sd.rel.error.Q = sd(na.omit(rel.error.Q)), mean.rel.error.env = mean(na.omit(rel.error.env)), sd.rel.error.env = sd(na.omit(rel.error.env)),NAs = sum(is.na(iter)))

library(tidyr)
summary <- summary %>% pivot_longer(cols = c(mean.rel.error.Q, mean.rel.error.env), names_to = "parameter", values_to = "mean.rel.error")
summary <- summary %>% pivot_longer(cols = c(sd.rel.error.Q, sd.rel.error.env), names_to = "parameter2", values_to = "sd.rel.error")
summary <- summary %>% mutate(parameter = factor(parameter, levels = c("mean.rel.error.Q", "mean.rel.error.env"), labels=c("Ventilation rate", "Environmental CO2")))
summary <- summary %>% mutate(parameter2 = factor(parameter2, levels = c("sd.rel.error.Q", "sd.rel.error.env"), labels=c("Ventilation rate", "Environmental CO2")))
summary <- summary %>% filter(parameter == parameter2)

ggplot(summary, aes(x=Q.true, y=mean.rel.error, color=parameter))+     geom_point() +
  geom_line() +
  geom_ribbon(aes(ymin=mean.rel.error-sd.rel.error, ymax=mean.rel.error+sd.rel.error), alpha=.2) + facet_wrap(~paste('Simulated environmental CO2:', envCO2.true))+
  labs(x="True Ventilation Rate (m3/s)", y="Relative error in estimated parameters") + theme_bw()  + labs(color="Parameter")

