#'exp_periods() Identify periods of exponential growth / decay from time series CO2 data
#'@param time_vec : Vector of time values. Must have POSIXct attribute. Must call in vectors
#'from dataframes using the $ operator so that the data keeps the POSIXct attribute
#'@param co2_vec : Vector of CO2 values. Must be in ppm, although there is no way for the
#'program to check this
#'@param quant_cutoff : Percentile cutoff for identifying exponential periods. Variance within
#'the CO2 values p index values to the left and right is calculated. High values of variance indicate
#'an exponential period. Percentile based cutoff is made to identify these high values, and therefore,
#'the exponential periods. A higher value for this cutoff is more strict, and may identify less exponential
#'periods. A lower values for this cutoff is less strict, and may identify more exponential periods.
#'@returns a data frame with the following columns:
#'        *Index Start : Index value for the start of the exponential period
#'        *Index Stop : Index value for the stop of the exponential period
#'        *Type : Type of Exponential period. Either Growth or Decay
#'        *Time Start : The POSIXct string time representation of the start time
#'        *Time Stop : The POSIXct string time representation of the stop time
#'        *Interval P Size : Gives the neighborhood interval size on both sides of a point. If p = 2,
#'        then the interval used for a point at index i goes from i - 2 to i + 2
        
exp_periods = function(time_vec, co2_vec, quant_cutoff = 0.9){
  #We have two vector inputs. One for our carbon dioxide data, one for our time data. Time data must be of the class POSIXct.
  #Quant_cutoff is the decimal representation of the quantile cutoff for variance. The base value is 0.9. 
  #We will first make sure that the input variables are of the correct format
  
  if(abs(quant_cutoff) > 1){
    stop("quant_cutoff must be a value between 0 and 1.")
  }
  
  #Check to make sure they are one dimensional
  if((!is.null(dim(time_vec)) | (!is.null(dim(co2_vec))) )){
    stop("Please make sure inputs are vectors. If calling from a data frame, please use the $ operator for 
    inputing a column as a vector, as this will keep the proper attributes attached to the data.")
  }
    
  #Check to make sure that time_vec is of class POSIXct
  if(inherits(time_vec, "POSIXct") == FALSE){
    stop("Time series data must be of the POSIXct type.If calling from a data frame, please use the $ operator for 
    inputing a column as a vector, as this will keep the proper attributes attached to the data.")
  }
  
  #Create data frame to hold our values
  data = data.frame("Time" = time_vec, "Carbon dioxide(ppm)" = co2_vec)
  #rename the column names so that they are correct
  colnames(data) = c("Time", "Carbon dioxide(ppm)")
  
  
  #Lets do a running p point forward average and variance
  #Need to find the lowest p value which causes smoothing
  #of the data. Can check for smoothness by doing u tests / t tests
  #between base data slope and p point slope distributions
  
  #Create data frames to hold slope, avg, and var
  p_slope = data.frame("base" = data$`Carbon dioxide(ppm)`)
  p_slope$'base'= c(0, (data$`Carbon dioxide(ppm)`[2:dim(data)[1]] -  data$`Carbon dioxide(ppm)`[1:dim(data)[1] - 1]) / 
                      as.numeric(difftime(data$`Time`[2:dim(data)[1]] , data$`Time`[1:dim(data)[1] - 1], units = 'secs')))
  p_avg = data.frame('base' = data$`Carbon dioxide(ppm)`)
  p_var = data.frame('base' = rep(0,dim(data)[1]))
  
  #Start the loop
  p = 1
  while (p > 0) {
    #initialize storage of data
    p_avg[, paste(p)] = p_slope$base 
    p_var[, paste(p)] = p_slope$base
    #for each index value, calculate the slope, avg, and var
    for (i in seq(1:dim(data)[1])){
      #get values
      #if i < p+1, just go from 1 to i + p
      #if i > dim(data)[1] - p, just go from i - p to dim(data)[1]
      if (i < p + 1){
        values = data$`Carbon dioxide(ppm)`[1:(i + p)]
      }else if (i > dim(data)[1] - p){
        values = data$`Carbon dioxide(ppm)`[(i - p):dim(data)[1]]
      } else {
        values = data$`Carbon dioxide(ppm)`[(i - p):(i + p)]
      }
      #assign values
      p_avg[i, paste(p)] = mean(values)
      p_var[i, paste(p)] = var(values)
      
    }
    #After having all the values, can calcualte slope
    p_slope[, paste(p)] = c(0, (p_avg[2:dim(data)[1], paste(p)] -  p_avg[1:dim(data)[1] - 1, paste(p)]) / 
                              as.numeric(difftime(data$`Time`[2:dim(data)[1]] , data$`Time`[1:dim(data)[1] - 1], units = 'secs')))
    #For p = 1, we want to skip
    if (p == 1){
      #iterate p
      p = p + 1
      next
    }
    #Our condition for stopping is if the distribution of slopes for p is:
    #Significantly different from base
    #Not significantly different from previous
    #In this case, we want the values from p - 1
    #Need to use KS test as this test is sensitive to changes in the shape of the distribution
    #The means of the distributions are not going to differ, but their shapes will
    test_base = ks.test(p_slope[, paste(p)], p_slope$base, exact = TRUE)
    test_prev = ks.test(p_slope[, paste(p)], p_slope[, paste(p - 1)], exact = TRUE)
    
    #If we meet our conditons
    if ((test_base$p.value < 0.05) & (test_prev$p.value > 0.05)){
      #save minimal p value
      p = p - 1
      break
    }
    
    p = p + 1
  }
  
  #Save relevant p value data to our main data frame
  data$'Interval Avg' = p_avg[, paste(p)]
  data$'Interval Var' = p_var[, paste(p)]
  data$'Interval Slope' = p_slope[, paste(p)]
  
  #save cutoff value based on argument input.
  cutoff = quantile(data$`Interval Var`, probs = quant_cutoff)
  
  #Create data frame to hold our exponential growth and decay periods
  useful_data = data.frame("Index Start" = c(), "Index Stop" = c(), "Type" = c(), "Time Start" = c(), "Time Stop" = c())
  
  #get indexes
  indexes = which(data$`Interval Var` >= cutoff)
  #Sequential indexes will be part of the same growth / decay process. Can assign these indexes the same number
  #Want to know when indexes[i] != indexes[i + 1]. Can use these index values as our index values for our points.
  #Need to call back from indexes to get the correct index from our data
  indexes = c(indexes[which(indexes[1: length(indexes) - 1] + 1 != indexes[2 : length(indexes)])], indexes[length(indexes)])
  
  #From each index, we need to go backwards in time and forwards in time until the slope of our
  #Interval Var changes from positive to negative or negative  to positive - need to find the 
  #adjacent local maxima and minima. Also want to make sure that when this happens, we are below
  #our 90% quantile cutoff
  
  
  for (ind in indexes){
    #Get Slope
    #Need to check if index == 1. If this is the case, we act as if index == 2. This will
    #prevent errors when index is one.
    if(ind == 1){
      ind = 2
    }
    slope = data$`Interval Slope`[ind]
    #If slope is positive
    if (slope > 0){
      #If slope is positive, type is growth
      type_val = "Growth"
      
      #Get the slope of all points preceding that are non-positive
      temp = which(data$`Interval Slope`[1:ind] <= 0)
      
      #Need to throw an exception if temp has no entries. This means that 
      #the exponential process starts at the beginning of the data
      if (length(temp) == 0){
        start = 1
      }else{
        #The last value of temp will give the index of the value at which the 
        #exponential process started
        start = temp[length(temp)]
      }
      
      #Now, we do the same, but in the forward direction
      temp = which(data$`Interval Slope`[ind:dim(data)[1]] <= 0)
      
      #Need to throw an exception if temp has no entries. This means that 
      #the exponential process starts at the end of the data
      if (length(temp) == 0){
        stop = dim(data)[1]
      }else{
        #This time, we want the first value that meets the critera. However, we need to add
        #a value of ind - 2 to get the correct index value
        stop = temp[1] + ind - 2
      }
      
      #Want to check. If start >= stop, then we do not want to return that period. simply continue to next iteration
      #of the loop
      if (start >= stop){
        next
      }
      
      #Create dummy dataframe to append values
      dummy = data.frame("Index Start" = c(start), "Index Stop" = c(stop), "Type" = c(type_val),
                         "Time Start" = data$`Time`[start], "Time Stop" = data$`Time`[stop])
      
      #Append to existing Data Frame
      useful_data = rbind(useful_data, dummy)
      
    }
    else{ 
      #othwerwise, slope is negative
      type_val = "Decay"
      
      #Get the slope of all points preceding that are non-negative
      temp = which(data$`Interval Slope`[1:ind] >= 0)
      
      #Need to throw an exception if temp has no entries. This means that 
      #the exponential process starts at the beginning of the data
      if (length(temp) == 0){
        start = 1
      }else{
        #The last value of temp will give the index of the value at which the 
        #exponential process started
        start = temp[length(temp)]
      }
      
      #Now, we do the same, but in the forward direction
      temp = which(data$`Interval Slope`[ind:dim(data)[1]] >= 0)
      
      
      #Need to throw an exception if temp has no entries. This means that 
      #the exponential process starts at the end of the data
      if (length(temp) == 0){
        stop = dim(data)[1]
      }else{
        #This time, we want the first value that meets the critera. However, we need to add
        #a value of ind - 2 to get the correct index value
        stop = temp[1] + ind - 2
      }
      
      #Want to check. If start >= stop, then we do not want to return that period. simply continue to next iteration
      #of the loop
      if (start >= stop){
        next
      }
      
      #Create dummy dataframe to append values
      dummy = data.frame("Index Start" = c(start), "Index Stop" = c(stop), "Type" = c(type_val),
                         "Time Start" = data$`Time`[start], "Time Stop" = data$`Time`[stop])
      
      #Append to existing Data Frame
      useful_data = rbind(useful_data, dummy)
    }
    
  }
  
  #Sometimes there are duplicates. We can remove duplicates using unique
  useful_data = unique(useful_data)
  
  #Lets also assign a new column that simply returns the neighborhood size used for calculating interval calculations
  useful_data$'Interval.P.Size' = rep(p, dim(useful_data)[1])
  
  return(useful_data)
  

}
