\name{build_up}
\alias{build_up}
\title{Calculate Air Exchange Rate Impact on CO2 Concentration}
\description{
  Computes the impact of the air exchange rate (AER) on CO2 concentration in a given space. This function serves as the objective function to find roots with respect to the AER, playing a crucial role in air quality analysis involving CO2 buildup.
}
\usage{
build_up(aer, n, Gp, V, C0, C1, delta_t, Cr = 400)
}
\arguments{
  \item{aer}{Air Exchange Rate (1/hour).}
  \item{n}{Number of individuals.}
  \item{Gp}{Average age-adjusted CO2 generation rate (L/(min·person)).}
  \item{V}{Room volume in cubic meters.}
  \item{C0}{Initial CO2 concentration (ppm).}
  \item{C1}{Final CO2 concentration (ppm).}
  \item{delta_t}{Time difference in hours.}
  \item{Cr}{Replacement air CO2 concentration (ppm), default is 400 ppm.}
}
\details{
  The function calculates the buildup of CO2 in a room based on several factors including the air exchange rate, the number of individuals present, their CO2 generation rate, the volume of the room, and the initial and final CO2 concentrations. The equation used in this function is central to understanding how various factors influence the rate at which CO2 accumulates in a given space over time.
}
\value{
  Returns the calculated value of the CO2 concentration equation based on the input parameters.
}
\references{
  %% Add references to literature or websites related to the algorithm or context here %%
}
\author{
  %% Your Name / Your Institution %%
}
\note{
  The accuracy of the results depends on the precision of the input parameters.
}
\seealso{
  \code{\link{secant_method}}, \code{\link{secant_starter}}
}
\examples{
## Example of using build_up
## Not run:
# build_up(aer = 0.5, n = 5, Gp = 0.08, V = 50, C0 = 400, C1 = 800, delta_t = 1, Cr = 400)
## End(Not run)
}
\keyword{math}
\keyword{environment}
