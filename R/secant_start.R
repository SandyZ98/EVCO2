
#' secant_starter() Find appropriate starting values for the secant method based on parameter to estimate.
#' @param aer : Air Exchange Rate ( 1 / hour)
#' @param n : number of individuals
#' @param Gp : Average age-adjusted CO2 generation rate (L / (min·person))
#' @param V : Room volume in meters cubed
#' @param C0 : Initial CO2 concentration (ppm)
#' @param C1 : Final CO2 concentration (ppm)
#' @param delta_t : Time difference (hours)
#' @param Cr : Replacement Air CO2 concentration (ppm). Will set to a default value of ambient air (400 ppm)
#' @returns A list with the following attributes:
#'        *start : start value for the secant method
#'        *stop : stop value for the secant method
#'        *asymptote : value of x-axis asymptote. Mainly used for debugging. If estimating parameters in the denominator
#'        of the temp variable in build_up(), Volume or Air Exchange Rate, x-axis value must be selected to be smaller
#'        than the asymptote. If estimating values in the numerator of temp, number of individuals (n) or Gp, x-axis
#'        value must be selected to be larger than the asymptote.
secant_starter = function(aer = NA, n = NA, Gp = NA, V = NA, C0, C1, delta_t, Cr = 400){

  #The first thing we need to do is:
  #1. Make sure we only have one variable to estimate. If not, stop the function and return error message
  #2. If we only have one variable to estimate, determine which variable that is.

  #To check 1. we can sum the values of is.na(variable) for n, Gp, V, and aer. If this value is not equal to 1,
  #then we will stop the function and return an error message. First, we will store these values in boolean
  #variables, as we can use these variables to check for 2.
  e_n = is.na(n)
  e_Gp = is.na(Gp)
  e_V = is.na(V)
  e_aer = is.na(aer)

  num_NA = sum(e_n, e_Gp, e_V, e_aer)

  if (num_NA != 1){
    stop("One dimensional parameter estimation requires that only one parameter is unknown / only one parameter
         is being estimated. Ensure that of the arguments n, Gp, V, and aer, that only one is not given an
         argument value. If a value must be assigned, assign NA as the value for the argument to be estimated")
  }


  #Now, we need to calculate the asymptote. The calculation changes depending on which variable is missing. Thus,
  #We have 4 if statements, each with a different calculation.

  asymptote = 0

  if(e_n){
    asymptote = (V * aer * (C1 - Cr)) / (6e4 * Gp)
  }else if(e_Gp){
    asymptote = (V * aer * (C1 - Cr)) / (6e4 * n)
  }else if(e_V){
    asymptote = (6e4 * n * Gp) / (aer*(C1 - Cr))
  }else if(e_aer){
    asymptote = (6e4 * n * Gp) / (V*(C1 - Cr))
  }


  #Initialize value_start and value_stop as 0. We will try moving away from asymptote by step size until we calculate
  #negative values for both value_start and value_stop. We will initialize step size as either half of the value of
  #asymptote or 1, whichever value is smaller. Depending on the parameter to be estimated, we either want step_size
  #to be positive or negative.
  step_size = min((asymptote / 2), 1)
  #For V and aer, we need to subtract from the asymptote. So assign negative value if estimating either parameter
  if (e_V | e_aer){
    step_size = -1 * step_size
  }
  value_start = 0
  value_stop = 0


  #We only assign value_stop a value if we have found 2 valid points. If build_up(value_start) < 0, then
  #we guarantee that asymptote + / - (step_size)/2 is also < 0. Thus, we have our two points.
  #While we haven't found values. Completely separate into 4 different while loops to speed up computation
  #time. It would be slower to do our if statement check each loop iteration, rather than having our if
  #statments first. Only downside is the code takes up much more space, as there are 4 while loops instead of one.

  if(e_n){
    while (value_stop == 0){

      #assign values
      value_start = asymptote + step_size
      n = value_start

      #calculate bulid_up(value_start)
      f_start = build_up(aer, n, Gp, V, C0, C1, delta_t, Cr)

      #if less than 0
      if (f_start < 0){
        #assign value_stop the following value
        value_stop = asymptote + (step_size / 2)
      }else{#else,
        #decrease step size by half
        step_size = step_size / 2
      }

    }
  }else if(e_Gp){
    while (value_stop == 0){

      #assign values
      value_start = asymptote + step_size
      Gp = value_start

      #calculate bulid_up(value_start)
      f_start = build_up(aer, n, Gp, V, C0, C1, delta_t, Cr)

      #if less than 0
      if (f_start < 0){
        #assign value_stop the following value
        value_stop = asymptote + (step_size / 2)
      }else{#else,
        #decrease step size by half
        step_size = step_size / 2
      }

    }
  }else if(e_V){
    while (value_stop == 0){

      #assign values
      value_start = asymptote + step_size
      V = value_start

      #calculate bulid_up(value_start)
      f_start = build_up(aer, n, Gp, V, C0, C1, delta_t, Cr)

      #if less than 0
      if (f_start < 0){
        #assign value_stop the following value
        value_stop = asymptote + (step_size / 2)
      }else{#else,
        #decrease step size by half
        step_size = step_size / 2
      }

    }
  }else if(e_aer){
    while (value_stop == 0){

      #assign values
      value_start = asymptote + step_size
      aer = value_start

      #calculate bulid_up(value_start)
      f_start = build_up(aer, n, Gp, V, C0, C1, delta_t, Cr)

      #if less than 0
      if (f_start < 0){
        #assign value_stop the following value
        value_stop = asymptote + (step_size / 2)
      }else{#else,
        #decrease step size by half
        step_size = step_size / 2
      }

    }
  }

  return(list(start = value_start,
              stop = value_stop,
              asymptote = asymptote))

}
