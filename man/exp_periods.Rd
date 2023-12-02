\name{exp_periods}
\alias{exp_periods}
\title{Identify Exponential Growth/Decay Periods in CO2 Time Series Data}
\description{
  Identifies periods of exponential growth or decay in time series data of CO2 concentrations. This function analyzes time-sequenced CO2 data to detect significant changes, categorizing them into growth or decay phases based on variance and trend analysis.
}
\usage{
exp_periods(time_vec, co2_vec, quant_cutoff = 0.9)
}
\arguments{
  \item{time_vec}{Vector of time values with POSIXct attribute.}
  \item{co2_vec}{Vector of CO2 values in ppm.}
  \item{quant_cutoff}{Percentile cutoff for identifying exponential periods, determining the variance threshold for significant changes in CO2 levels.}
}
\details{
  The function works by analyzing the variance within the CO2 values at each point in the time series, using a quantile-based cutoff to identify periods of high variance indicative of exponential changes. It is crucial that the input vectors maintain their appropriate formats and attributes, especially the POSIXct attribute for time data.
}
\value{
  Returns a data frame with the following columns:
  \item{Index Start}{Index value at the start of the exponential period.}
  \item{Index Stop}{Index value at the end of the exponential period.}
  \item{Type}{Type of exponential period, either 'Growth' or 'Decay'.}
  \item{Time Start}{POSIXct string time representation of the start time.}
  \item{Time Stop}{POSIXct string time representation of the stop time.}
}
\references{
  %% Add references to literature or websites related to time series analysis and CO2 data interpretation here %%
}
\author{
  %% Your Name / Your Institution %%
}
\note{
  Ensure that the time data is in the POSIXct format and that CO2 values are correctly inputted in ppm.
}
\seealso{
  \code{\link{other_related_function}}, \code{\link{another_function}}
}
\examples{
## Example of using exp_periods
## Not run:
# exp_periods(time_vector, co2_vector, quant_cutoff = 0.9)
## End(Not run)
}
\keyword{time series}
\keyword{environment}
