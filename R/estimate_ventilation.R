
#' Estimate Ventilation Rate
#'
#' This function estimates ventilation rates in indoor environments based on CO2 time-series data.
#'
#' @param freq The frequency of CO2 measurements (e.g., measurements per hour).
#' @param CO2 A numeric vector representing the CO2 concentrations over time.
#' @param volume The volume of the room in cubic meters.
#' @param envCO2known Known environmental CO2 concentration in ug/m^3.
#' @param init.Q Initial estimate of ventilation rate.
#' @param temp The temperature in degrees Celsius (default is 25).
#' @param persondata A data frame containing information about individuals in the room, used for Newton-Raphson optimization.
#' @param critpoints Critical points for multi-dimensional Newton method.
#' @param method The optimization method to use ('NR' for Newton-Raphson, 'Newton' for multi-dimensional Newton).
#' @param max.iter Maximum number of iterations for optimization.
#' @param tol Tolerance for convergence in optimization.
#' @param E.init Initial estimate of emissions rate for multi-dimensional Newton method.
#' @param envCO2.init Initial estimate of environmental CO2 for multi-dimensional Newton method.
#' @param verbose If TRUE, print iteration details during optimization.
#' @param record.steps If TRUE, record optimization steps.
#'
#' @return A list containing the estimated ventilation rate, number of iterations, and convergence status.
#'
#' @examples
#' # Example using Newton-Raphson method
#' estimate_ventilation(freq = 1, CO2 = c(400, 600, 800), volume = 100,
#'                      envCO2known = 500, init.Q = 1, temp = 25,
#'                      persondata = data.frame(time = c(0, 1, 2), activity = c(1, 2, 1)),
#'                      method = 'NR')
#'
#' # Example using multi-dimensional Newton method
#' estimate_ventilation(freq = 1, CO2 = c(400, 600, 800), volume = 100,
#'                      envCO2known = 500, init.Q = 1, temp = 25,
#'                      persondata = data.frame(time = c(0, 1, 2), activity = c(1, 2, 1)),
#'                      method = 'Newton')
#'
#' @export
#'
estimate_ventilation <- function(freq,
                                 CO2,
                                 volume,
                                 envCO2known=NULL,
                                 init.Q,
                                 temp=25,
                                 persondata=NULL,
                                 critpoints=NULL,
                                 method='NR',
                                 max.iter=1000,
                                 tol=1e-8,
                                 E.init = NULL,
                                 envCO2.init=NULL,
                                 verbose=FALSE, record.steps=FALSE) {

  # format data
  CO2 = as.numeric(CO2)

  # get emissions from persondata
  if (is.null(persondata) & method=="NR") {
    stop("persondata must be provided for 1D Newton-Raphson optimization")
  }

  if(!is.null(persondata)){

    emissions <- persondata_to_emission(persondata, temp, freq)
    times <- emissions$times
    nadj_CO2rate <- emissions$nadj_CO2rate

    if(length(CO2) > length(times)) {
      warning(sprintf("More CO2 measurements than expected from persondata. Make sure the persondata starts at 0 and covers the same time period as the CO2 data. Using the last row of persondata for the last %d CO2 measurements", length(CO2) - length(times)))
      nadj_CO2rate <- c(nadj_CO2rate, rep(nadj_CO2rate[length(nadj_CO2rate)], length(CO2) - length(times)))
    } else if(length(CO2) < length(times)) {
      warning("More persondata than CO2 measurements. Make sure the persondata starts at 0 and covers the same time period as the CO2 data. Only using the provided persondata up until the last CO2 measurement")
    }
    nadj_CO2rate <- nadj_CO2rate[1:length(CO2)]
  }

  # Convert from ppm to ug/m^3
  ppm_to_ug <- function(ppm, temp){
    ppm / 44.01 * 8314 * (273.15+temp) / 1000 / 101.325
  }

  # convert from ug/m^3 to ppm
  ug_to_ppm <- function(ug, temp) {
    ug * 44.01 * 1000 * 101.325 / 8314 / (273.15+temp)
  }


  if(!is.null(envCO2known)){
    envCO2known <- ppm_to_ug(envCO2known, temp)
  }else{
    envCO2.init <- ppm_to_ug(envCO2.init, temp)
  }
  CO2 <- ppm_to_ug(CO2, temp)

  # # Gradient descent
  #   Q = init.Q
  #   convergence = 0
  #   for(iter in 1:max.iter){
  #     ChatminusC = (nadj_CO2rate[1:(length(CO2)-1)]/Q*(1-exp(-Q/volume*freq))+(CO2[1:(length(CO2)-1)] - envCO2)*exp(-Q/volume*freq) + envCO2 - CO2[2:length(CO2)])
  #     f = sum(ChatminusC^2)
  #     #df = 2*(ChatminusC %*% ( -nadj_CO2rate[1:(length(CO2)-1)]/(Q^2)*(1-exp(-Q/volume*freq)) + (nadj_CO2rate[1:(length(CO2)-1)]/Q - CO2[1:(length(CO2)-1)] + envCO2)*exp(-Q/volume*freq)*freq/volume))
  #     #df = 2*sum(ChatminusC *(-nadj_CO2rate[1:(length(CO2)-1)]/(Q^2) * (1-exp(-Q/volume * freq)) + exp(-Q / volume * freq) * (freq / volume) * (nadj_CO2rate[1:(length(CO2)-1)] / Q - CO2[1:(length(CO2)-1)] + envCO2)))
  #     df = 2*sum(ChatminusC *(-nadj_CO2rate[1:(length(CO2)-1)]/(Q^2) * (1-exp(-Q/volume * freq)) + exp(-Q / volume * freq) * (freq / volume) * (nadj_CO2rate[1:(length(CO2)-1)] / Q - CO2[1:(length(CO2)-1)] + envCO2)))
  #     Q2 = Q - 1e-4*f/df
  #     # print(paste("f:",f))
  #     # print(paste("df:",df))
  #     # print(paste("Q2:",Q2))
  #     if(Q2 < 0){
  #       warning("Large jump resulted in negative Q. Returning previous Q")
  #       break()
  #     }
  #     if(abs(Q2-Q) < tol){
  #       convergence=1
  #       break()
  #     }
  #     Q = Q2
  #   }
  # }

  # Newton-Raphson
  if(method == 'NR'){
    if (verbose) {
      print("Running Newton-Raphson")
    }

    if (record.steps) {
      steps <- data.frame(iter=numeric(max.iter),
                          Q=numeric(max.iter), # estimates for ventilation rate
                          f=numeric(max.iter), # sum of squares
                          df=numeric(max.iter), # derivative of SS
                          d2f=numeric(max.iter)) # 2nd derivative of SS
    }

    Q = init.Q
    convergence = 0
    E = nadj_CO2rate[1:(length(CO2)-1)]
    Ci_1 = CO2[1:(length(CO2)-1)]
    Ci = CO2[2:length(CO2)]
    envCO2 = envCO2known
    for(iter in 1:max.iter){
      expQ = exp(-Q/volume*freq)
      Chat_C = E/Q*(1-expQ)+(Ci_1 - envCO2)*expQ + envCO2 - Ci

      derivChat_C = -E/(Q^2) * (1-expQ) + expQ * (freq / volume) * (E/Q - Ci_1 + envCO2)

      f = 2*sum(Chat_C*derivChat_C)
      df = 2*sum(derivChat_C^2 + Chat_C*(
        2*E/(Q^3)*(1-expQ) -
          expQ * (freq/volume)^2 * (E/Q - Ci_1 + envCO2) +
          2*E/Q^2 * expQ *(-freq/volume)
      ))
      Q2 = Q - f/df

      if (verbose) {
        print(paste("Iter:",iter,
                    "Q:", Q2,
                    "df:",f,
                    "d2f:",df))
      }

      if (record.steps) {
        steps[iter,] <- c(iter=iter,
                          Q=Q,
                          f=Chat_C %*% Chat_C,
                          df=f,
                          d2f=df)
      }

      if(abs(Q-Q2) < tol){
        convergence=1
        break()
      }
      Q=Q2
    }
    if (record.steps) {
      steps=steps[1:iter,]

      return(list(est_ventilation = Q,
                  iter=iter,
                  convergence=convergence,
                  steps=steps))
    } else {
      return(list(est_ventilation = Q,
                  iter=iter,
                  convergence=convergence))
    }
  }

  if(method == 'Newton'){
    if (verbose) {
      print("Running Multi-dimensional Newton Method")
    }
    return(multi_Newton(persondata=persondata,
                        max.iter=max.iter, tol=tol,
                        verbose=verbose, record.steps=record.steps,
                        critpoints=critpoints, nadj_CO2rate=nadj_CO2rate,
                        CO2=CO2, init.Q=init.Q, envCO2known=envCO2known,
                        E.init = E.init, envCO2.init = envCO2.init, freq=freq, temp=temp))
  }
}
