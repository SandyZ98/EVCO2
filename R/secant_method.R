
#' secant_method() build up method for estimating parameters. Estimates parameter value by solving for the
#' root of the build_up() function
#' @param aer : Air Exchange Rate ( 1 / hour)
#' @param n : number of individuals
#' @param Gp : Average age-adjusted CO2 generation rate (L / (min·person))
#' @param V : Room volume in meters cubed
#' @param C0 : Initial CO2 concentration (ppm)
#' @param C1 : Final CO2 concentration (ppm)
#' @param delta_t : Time difference (hours)
#' @param Cr : Replacement Air CO2 concentration (ppm). Will set to a default value of ambient air (400 ppm)
#' @param tol : absolute difference of function values between steps we use as a stopping condition
#' @param max_iter : maximum number of iterations
#' @returns : A list containing the following attributes:
#'    * root - Parameter value with build_up(parameter) close to zero
#'    * f_root - build_up(root)
#'    * iter - number of iterations to reach the solution
#'    * convergence - 0 if the root was found successfully, 1 if not found
secant_method = function(aer = NA, n = NA, Gp = NA, V = NA, C0, C1, delta_t, Cr = 400, tol=1e-10, max_iter=1000){

  #set convergence to 1
  convergence = 1

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

  #Next, we need to find appropriate starting conditions for the method.
  starter_vals = secant_starter(aer, n, Gp, V, C0, C1, delta_t, Cr)
  value_0 = starter_vals$start
  value_1 = starter_vals$stop

  #calculate our function values for build_up(). This depends on the parameter being estimated. Assign
  #appropriate

  if(e_n){
    n = value_0
  }else if(e_Gp){
    Gp = value_0
  }else if(e_V){
    V = value_0
  }else if(e_aer){
    aer = value_0
  }
  f0 = build_up(aer, n, Gp, V, C0, C1, delta_t, Cr)

  if(e_n){
    n = value_1
  }else if(e_Gp){
    Gp = value_1
  }else if(e_V){
    V = value_1
  }else if(e_aer){
    aer = value_1
  }
  f1 = build_up(aer, n, Gp, V, C0, C1, delta_t, Cr)

  #Calculate our change in aer
  delta_value = -f1 / (f1 - f0)*(value_1 - value_0)
  #interpolate to get our value_2 value
  value_2 = value_1 + delta_value

  #Separate loop for each parameter to be estimated. The only difference between each loop is that
  #a different parameter is updated with values of the new value_1, as we only want to update the
  #parameter that is being estimated. For each parameter:
  if(e_n){
    #Now we are ready to iterate! For each iteration
    for(iter in 1:max_iter){
      #check to see if we have reached convergence. If so, we will return value_2 as our root
      if(abs(delta_value)<tol){
        convergence = 0
        break
      }
      #update values for this iteration. value_0 = value_1, value_1 = value_2. Calculate value_2
      f0 = f1
      value_1 = value_2
      #Assign parameter to be estimated the new value
      n = value_1
      #Calculate new f1 value.
      f1 = build_up(aer, n, Gp, V, C0, C1, delta_t, Cr)
      #can use difference in function values for simple calculation of change in aer
      delta_f = f1 - f0
      #New change in aer = (f1 / change in f) * change in aer
      delta_value = -(f1/delta_f)*delta_value
      #Calculate value_2 via interpolation
      value_2 = value_1 + delta_value
    }
  }else if(e_Gp){
    #Now we are ready to iterate! For each iteration
    for(iter in 1:max_iter){
      #check to see if we have reached convergence. If so, we will return value_2 as our root
      if(abs(delta_value)<tol){
        convergence = 0
        break
      }
      #update values for this iteration. value_0 = value_1, value_1 = value_2. Calculate value_2
      f0 = f1
      value_1 = value_2
      #Assign parameter to be estimated the new value
      Gp = value_1
      #Calculate new f1 value.
      f1 = build_up(aer, n, Gp, V, C0, C1, delta_t, Cr)
      #can use difference in function values for simple calculation of change in aer
      delta_f = f1 - f0
      #New change in aer = (f1 / change in f) * change in aer
      delta_value = -(f1/delta_f)*delta_value
      #Calculate value_2 via interpolation
      value_2 = value_1 + delta_value
    }
  }else if(e_V){
    #Now we are ready to iterate! For each iteration
    for(iter in 1:max_iter){
      #check to see if we have reached convergence. If so, we will return value_2 as our root
      if(abs(delta_value)<tol){
        convergence = 0
        break
      }
      #update values for this iteration. value_0 = value_1, value_1 = value_2. Calculate value_2
      f0 = f1
      value_1 = value_2
      #Assign parameter to be estimated the new value
      V = value_1
      #Calculate new f1 value.
      f1 = build_up(aer, n, Gp, V, C0, C1, delta_t, Cr)
      #can use difference in function values for simple calculation of change in aer
      delta_f = f1 - f0
      #New change in aer = (f1 / change in f) * change in aer
      delta_value = -(f1/delta_f)*delta_value
      #Calculate value_2 via interpolation
      value_2 = value_1 + delta_value
    }
  }else if(e_aer){
    #Now we are ready to iterate! For each iteration
    for(iter in 1:max_iter){
      #check to see if we have reached convergence. If so, we will return value_2 as our root
      if(abs(delta_value)<tol){
        convergence = 0
        break
      }
      #update values for this iteration. value_0 = value_1, value_1 = value_2. Calculate value_2
      f0 = f1
      value_1 = value_2
      #Assign parameter to be estimated the new value
      aer = value_1
      #Calculate new f1 value.
      f1 = build_up(aer, n, Gp, V, C0, C1, delta_t, Cr)
      #can use difference in function values for simple calculation of change in aer
      delta_f = f1 - f0
      #New change in aer = (f1 / change in f) * change in aer
      delta_value = -(f1/delta_f)*delta_value
      #Calculate value_2 via interpolation
      value_2 = value_1 + delta_value
    }
  }

  #One last calculation. Assign value_2 to the correct parameter. Function call will
  #occur in the return list
  if(e_n){
    n = value_2
  }else if(e_Gp){
    Gp = value_2
  }else if(e_V){
    V = value_2
  }else if(e_aer){
    aer = value_2
  }

  return(list(root=value_2,
              f_root = build_up(aer, n, Gp, V, C0, C1, delta_t, Cr),
              iter=iter,
              convergence=convergence))

}
