\name{Gp_met}
\alias{Gp_met}
\title{Estimate Gp from Height, Weight, and Metabolic Equivalent}
\description{
  Calculates the average age-adjusted CO2 generation rate (Gp) based on an individual's height, weight, and their metabolic equivalent value (MET), which reflects their level of physical exertion.
}
\usage{
Gp_met(height, weight, met)
}
\arguments{
  \item{height}{Height of the individual in meters.}
  \item{weight}{Weight of the individual in kilograms.}
  \item{met}{Metabolic equivalent value for the level of exertion.}
}
\details{
  The function employs a formula that incorporates body surface area (calculated from height and weight) and the metabolic equivalent value to estimate the CO2 generation rate. This rate is crucial in understanding the individual's contribution to CO2 levels in a given environment.
}
\value{
  Returns the calculated Gp value, representing the CO2 generation rate in L/(min·person).
}
\references{
  %% Add references to literature or websites related to metabolic calculations and environmental studies here %%
}
\author{
  %% Your Name / Your Institution %%
}
\note{
  Accurate measurements of height and weight are essential for precise Gp estimation.
}
\seealso{
  \code{\link{another_related_function}}, \code{\link{a_complementary_function}}
}
\examples{
## Example of using Gp_met
## Not run:
# Gp_met(height = 1.75, weight = 70, met = 1.2)
## End(Not run)
}
\keyword{metabolism}
\keyword{environment}
