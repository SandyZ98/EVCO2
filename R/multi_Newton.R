
multi_Newton <- function(persondata, max.iter, tol, verbose, record.steps, critpoints,
                         nadj_CO2rate, CO2, init.Q, envCO2known, E.init, envCO2.init, freq, temp){
  # convert from ug/m^3 to ppm
  ug_to_ppm <- function(ug, temp) {
    ug * 44.01 * 1000 * 101.325 / 8314 / (273.15+temp)
  }

  if(is.null(envCO2known)){
    estimate_envCO2 <- TRUE
  }else{
    estimate_envCO2 <- FALSE
  }
  if(is.null(persondata)){
    if(is.null(critpoints)) {
      stop('If you want to estimate emission rates, must provide a vector of indices corresponding to times where people enter/exit the room')
    }
    estimate_E <- TRUE
    nadj_CO2rate <- NULL
  }else estimate_E = FALSE

  convergence = 0
  Gradient_summand = function(){ # calculates the summands of the gradient
    if(estimate_E){
      dSSdE <- (Chat_C/Q*(1-q))
    } else {
      dSSdE <- 0
    }
    dSSdQ <- (Chat_C*derivChat_C)
    if(estimate_envCO2){
      dSSdenvCO2 <- (Chat_C*(1-q))
    } else {
      dSSdenvCO2 <- 0
    }
    return(list(dSSdE, dSSdQ, dSSdenvCO2))
  }
  Hessian_summand = function(){ # calculates the summands of the hessian
    if(estimate_E){
      d2SSdE2 <- ((1/Q*(1-q))  ^2)
      d2SSdEdQ <- (
        Chat_C*(-1/Q^2*(1-q) + 1/Q*(freq/volume)*q) +
          1/Q*(1-q)*derivChat_C)
    } else {
      d2SSdE2 <- 0
      d2SSdEdQ <- 0
    }
    d2SSdQ2 <- (derivChat_C^2 + Chat_C*(
      2*E/(Q^3)*(1-q) -
        q * (freq/volume)^2 * (E/Q - Ci_1 + envCO2) +
        2*E/Q^2 * q *(-freq/volume)))
    if(estimate_envCO2){
      d2SSdenvCO22 <- ((1-q)^2)
      d2SSdQdenvCO2 <- (
        Chat_C*(q*freq/volume) +
          (1-q)*derivChat_C)
    } else {
      d2SSdenvCO22 <- 0
      d2SSdQdenvCO2 <- 0
    }
    if(estimate_E & estimate_envCO2){
      d2SSdEdenvCO2 <- (1/Q*(1-q)^2)
    } else {
      d2SSdEdenvCO2 <- 0
    }
    return(list(d2SSdE2, d2SSdEdQ, d2SSdQ2, d2SSdEdenvCO2, d2SSdQdenvCO2, d2SSdenvCO22))
  }

  gradient = c(0) # initialize gradient as vector of length 1 (only estimating Q)
  hessiandim = 1 # initialize hessian as matrix of size 1x1 (only estimating Q)
  Eindices = list()

  if(estimate_E){
    # indices corresponding to the intervals over which we will estimate E separately
    Eindices[[1]] = 1:critpoints[1]
    for(i in 2:length(critpoints)){
      Eindices[[i]]= critpoints[i-1]:critpoints[i]
    }
    if(critpoints[length(critpoints)] < length(CO2)) { # add the last interval, if needed
      Eindices[[length(critpoints)+1]] = critpoints[length(critpoints)]:(length(CO2)-1)
    }
    E = rep(0, length(CO2)-1)
    if(!is.null(E.init)){
      if(length(E.init) != length(Eindices)){
        stop('Length of E.init must be equal to the number of intervals over which you are estimating E: number of critical points +1')
      }
    } else {
      E.init = (1:length(Eindices))/10 # choose something small-ish but non-zero
    }
    for(i in 1:length(Eindices)){
      E[Eindices[[i]]] = E.init[i]
    }
    gradient = c(gradient, rep(0, length(Eindices))) # add to gradient for each E we need to estimate
    hessiandim = hessiandim + length(Eindices) # add to hessian for each E we need to estimate
  } else {
    Eindices[[1]] = 1:(length(CO2)-1)
    E = nadj_CO2rate[1:(length(CO2)-1)]
  }

  if(estimate_envCO2){
    if(!is.null(envCO2.init)){
      envCO2 = envCO2.init
    } else {
      envCO2 = 222 # 400 ppm in ug
    }
    gradient = c(gradient, 0) # add to gradient for envCO2
    hessiandim = hessiandim + 1 # add to hessian for envCO2
  }
  Ci_1 = CO2[1:(length(CO2)-1)]
  Ci = CO2[2:length(CO2)]
  Q = init.Q
  hessian = matrix(0, nrow=hessiandim, ncol=hessiandim)
  q = exp(-Q/volume*freq)
  Chat_C = E/Q*(1-q)+(Ci_1 - envCO2)*q + envCO2 - Ci
  f0 = sum(Chat_C^2)
  if (record.steps) {
    steps <- data.frame(iter=numeric(max.iter),
                        Q=numeric(max.iter), # estimates for ventilation rate
                        envCO2=numeric(max.iter), # estimates for envCO2
                        f=numeric(max.iter), # sum of squares
                        sumdelta = numeric(max.iter) # step size (sum)
    )
    for(Enum in 1:length(Eindices)){
      steps[paste0('E', Enum)] = numeric(max.iter)
    }
  }

  for(iter in 1:max.iter){
    #q = exp(-Q/volume*freq)
    #Chat_C = E/Q*(1-q)+(Ci_1 - envCO2)*q + envCO2 - Ci
    #derivChat_C =(-E/Q^2)*(1-q) + E/Q*(freq/volume)*q + (Ci_1 - envCO2)*(-freq/volume)*q
    derivChat_C = -E/(Q^2) * (1-q) + q * (freq / volume) * (E/Q - Ci_1 + envCO2)

    gradientres <- Gradient_summand()
    dSSdE = gradientres[[1]]
    dSSdQ = gradientres[[2]]
    dSSdenvCO2 = gradientres[[3]]
    hessianres <- Hessian_summand()
    d2SSdE2 = hessianres[[1]]
    d2SSdEdQ = hessianres[[2]]
    d2SSdQ2 = hessianres[[3]]
    d2SSdEdenvCO2 = hessianres[[4]]
    d2SSdQdenvCO2 = hessianres[[5]]
    d2SSdenvCO22 = hessianres[[6]]
    gradient[1] = 2*sum(dSSdQ)
    hessian[1,1] = 2*sum(d2SSdQ2)
    # gradient[1] = 2*mean(dSSdQ)
    # hessian[1,1] = 2*mean(d2SSdQ2)
    if(estimate_E){
      for(i in 1:length(Eindices)){
        index = Eindices[[i]]
        gradient[i+1] = 2*sum(dSSdE[index])
        hessian[i+1, i+1] = 2*d2SSdE2*length(index)
        hessian[1, i+1] = 2*sum(d2SSdEdQ[index])
        # gradient[i+1] = 2*sum(dSSdE[index])/(length(CO2) - 1)
        # hessian[i+1, i+1] = 2*d2SSdE2*(length(index)/(length(CO2) - 1))
        # hessian[1, i+1] = 2*sum(d2SSdEdQ[index])/(length(CO2) - 1)
        hessian[i+1, 1] = hessian[1, i+1]
        if(estimate_envCO2){
          hessian[i+1, nrow(hessian)] = 2*d2SSdEdenvCO2*length(index)
          hessian[nrow(hessian), i+1] = hessian[i+1, nrow(hessian)]
        }
      }
    }
    if(estimate_envCO2){
      #gradient[length(gradient)] = 2*mean(dSSdenvCO2)
      gradient[length(gradient)] = 2*sum(dSSdenvCO2)
      hessian[nrow(hessian), nrow(hessian)] = 2*d2SSdenvCO22*(length(CO2)-1)
      hessian[1, nrow(hessian)] = 2*sum(d2SSdQdenvCO2)
      # hessian[nrow(hessian), nrow(hessian)] = 2*d2SSdenvCO22
      # hessian[1, nrow(hessian)] = 2*mean(d2SSdQdenvCO2)
      hessian[nrow(hessian), 1] = hessian[1, nrow(hessian)]
    }
    #delta =  solve(hessian + diag(1e-2, nrow=nrow(hessian)))%*%gradient
    delta =  solve(hessian)%*%gradient

    # update parameters
    #Qnew = max(Q - delta[1], 1e-6) # constrain to be positive
    Qnew = Q - delta[1]
    if(estimate_E){
      Enew = rep(0, length(E))
      for(i in 1:length(Eindices)){
        index = Eindices[[i]]
        Enew[index] = E[index] - delta[i+1]
        #Enew[index] = pmax(E[index] - delta[i+1], 0) # constrain to be positive
      }
    }
    if(estimate_envCO2){
      #envCO2new = max(envCO2 - delta[length(delta)], 0) # constrain to be positive
      envCO2new = envCO2 - delta[length(delta)]
    }

    # constraints : envCO2 between 375 and 425, all E >=0, ventilation between 0 and 1000
    # if((envCO2new <= 425) & (envCO2new >= 375) &
    #    (sum(Enew>=0) == length(Enew)) & (Qnew > 0) & (Qnew < 1000)){
    #   conditions_met = TRUE
    # } else {
    #   conditions_met = FALSE
    # }
    #
    # while(!conditions_met){
    #   delta = delta/2
    #   Qnew = Q - delta[1]
    #   if(estimate_E){
    #     for(i in 1:length(Eindices)){
    #       index = Eindices[[i]]
    #       Enew[index] = E[index] - delta[i+1]
    #     }
    #   }
    #   if(estimate_envCO2){
    #       envCO2new = envCO2 - delta[length(delta)]
    #   }
    #   if(envCO2new <= 425 & envCO2new >= 375 &
    #      sum(Enew>=0) == length(Enew) & Qnew > 0 & Qnew < 1000){
    #     conditions_met = TRUE
    #   } else conditions_met = FALSE
    # }
    Q = Qnew
    if(estimate_E) E = Enew
    if(estimate_envCO2) envCO2 = envCO2new
    q=exp(-Q/volume*freq)
    Chat_C = E/Q*(1-q)+(Ci_1 - envCO2)*q + envCO2 - Ci
    f1 = Chat_C%*%Chat_C
    if (verbose) {
      Es = rep(0, length(Eindices))
      for(i in 1:length(Eindices)){
        index = Eindices[[i]]
        Es[i] = paste("E", i, ":", E[index[1]])
      }
      print(paste("Iter:",iter,
                  "Q:", Q,
                  "envCO2:", ug_to_ppm(envCO2, temp),
                  "f:",f1,
                  "sum_delta:", sum(delta)))
      print(Es)
    }
    if (record.steps) {
      Etemp = rep(0, length(Eindices))
      for(i in 1:length(Eindices)){
        index = Eindices[[i]]
        Etemp[i] = E[index[1]]
      }
      steps[iter,] <- c(iter, Q, ug_to_ppm(envCO2, temp), f1, sum(delta), Etemp)

    }
    #if(abs(f1-f0)<=tol*(abs(f1)+abs(f0))){
    if(sqrt(gradient %*% gradient)<=tol){
      convergence = 1
      break
    }
    f0 = f1
  }
  if (record.steps) {
    steps=steps[1:iter,]

    return(list(ventilation = Q,
                E = E,
                envCO2 = ug_to_ppm(envCO2, temp),
                iter=iter,
                convergence=convergence,
                steps=steps))
  } else {
    return( list(ventilation = Q,
                 E = E,
                 envCO2 = ug_to_ppm(envCO2, temp),
                 iter=iter,
                 convergence=convergence))
  }
}
