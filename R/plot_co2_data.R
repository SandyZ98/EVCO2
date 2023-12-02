
#' Plot CO2 Data
#'
#' This function generates plots from the analyzed CO2 data.
#'
#' @param analyzed_data The analyzed data frame from `analyze_co2_data`.
#' @return A list of plots.
#' @export
#' @examples
#' analyzed_data <- analyze_co2_data("path/to/your/data.xlsx")
#' plot_co2_data(analyzed_data)

plot_co2_data <- function(analyzed_data) {
  plots <- list()

  # Ensure the input is as expected
  if ("data" %in% names(analyzed_data)) {
    data <- analyzed_data$data
  } else {
    stop("Invalid input: expected output from analyze_co2_data function")
  }

  # Plotting code
  par(mfrow = c(2, 2))

  plot1 <- plot(data$`Time(dd/mm/yyyy)`, data$`Interval Avg`, xlab = 'Time', ylab = 'Interval Average CO2 Concentration(ppm)',
                main = "Interval Average CO2 Concentration")
  plots[[1]] <- plot1

  plot2 <- plot(data$`Time(dd/mm/yyyy)`, data$`Carbon dioxide(ppm)`, xlab = 'Time', ylab = 'CO2 Concentration(ppm)',
                main = "CO2 Concentration")
  plots[[2]] <- plot2

  plot3 <- plot(data$`Time(dd/mm/yyyy)`, data$`Interval Var`, xlab = 'Time', ylab = 'Interval CO2 Variance',
                main = "Interval CO2 Variance")
  plots[[3]] <- plot3

  plot4 <- plot(data$`Time(dd/mm/yyyy)`, data$`Interval Slope`, xlab = 'Time', ylab = 'Interval CO2 Slope',
                main = "Interval CO2 Slope")
  plots[[4]] <- plot4

  return(plots)
}
