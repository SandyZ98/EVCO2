#' Analyze CO2 Data
#'
#' This function performs an analysis on CO2 data, calculating running averages, variances,
#' and identifying growth and decay periods.
#'
#' @param data The dataset containing the CO2 data.
#' @return A list containing analyzed data.
#' @export
#' @examples
#' analyze_co2_data("path/to/your/data.xlsx")

analyze_co2_data <- function(data) {
  data[,1] = as.POSIXct(data$`Time(dd/mm/yyyy)`, format = "%j/%m/%Y %I:%M:%S %p")

  #Lets do a running p point forward average and variance
  #Need to find the lowest p value which causes smoothing
  #of the data. Can check for smoothness by doing u tests / t tests
  #between base data slope and p point slope distributions

  #Create data frames to hold slope, avg, and var
  p_slope = data.frame("base" = data$`Carbon dioxide(ppm)`)
  p_slope$'base'= c(0, (data$`Carbon dioxide(ppm)`[2:dim(data)[1]] -  data$`Carbon dioxide(ppm)`[1:dim(data)[1] - 1]) /
                      as.numeric(difftime(data$`Time(dd/mm/yyyy)`[2:dim(data)[1]] , data$`Time(dd/mm/yyyy)`[1:dim(data)[1] - 1], units = 'secs')))
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
                              as.numeric(difftime(data$`Time(dd/mm/yyyy)`[2:dim(data)[1]] , data$`Time(dd/mm/yyyy)`[1:dim(data)[1] - 1], units = 'secs')))
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

  # Append p value data to the main data frame
  data$'Interval Avg' <- p_avg[, as.character(p)]
  data$'Interval Var' <- p_var[, as.character(p)]
  data$'Interval Slope' <- p_slope[, as.character(p)]
}



