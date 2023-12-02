\name{secant_starter}
\alias{secant_starter}
\title{Estimate Starting Values for the Secant Method}
\description{
  This function calculates appropriate starting values for the secant method used in one-dimensional parameter estimation within a given context, such as air quality analysis. It determines these values based on several parameters, including air exchange rate, number of individuals, CO2 generation rate, room volume, and CO2 concentrations.
}
\usage{
secant_starter(aer = NA, n = NA, Gp = NA, V = NA, C0, C1, delta_t, Cr = 400)
}
\arguments{
  \item{aer}{Air Exchange Rate (1/hour), can be NA for estimation.}
  \item{n}{Number of individuals, can be NA for estimation.}
  \item{Gp}{Average age-adjusted CO2 generation rate (L/(min·person)), can be NA for estimation.}
  \item{V}{Room volume in cubic meters, can be NA for estimation.}
  \item{C0}{Initial CO2 concentration (ppm).}
  \item{C1}{Final CO2 concentration (ppm).}
  \item{delta_t}{Time difference in hours.}
  \item{Cr}{Replacement Air CO2 concentration (ppm), default is 400 ppm.}
}
\details{
  The function first ensures that only one of the parameters (n, Gp, V, aer) is set for estimation (marked as NA). It then calculates the asymptote based on the missing variable. The starting and stopping values for the secant method are determined by iteratively adjusting around the asymptote.
}
\value{
  The function returns a list with three components:
  \item{start}{The start value for the secant method.}
  \item{stop}{The stop value for the secant method.}
  \item{asymptote}{The x-axis asymptote value, useful for debugging and ensuring correct parameter estimation direction.}
}
\references{
  %% Add references to literature or websites related to the algorithm or context here %%
}
\author{
  %% Your Name / Your Institution %%
}
\note{
  Ensure only one of the parameters n, Gp, V, or aer is set as NA for estimation.
}
\seealso{
  \code{\link{other_function}}, \code{\link{another_relevant_function}}
}
\examples{
## Example of using secant_starter
## Not run:
# secant_starter(aer = 0.5, n = NA, Gp = 0.08, V = 50, C0 = 400, C1 = 800, delta_t = 1)
## End(Not run)
}
\keyword{math}
\keyword{optimization}
