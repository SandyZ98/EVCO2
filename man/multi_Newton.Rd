\name{multi_Newton}
\alias{multi_Newton}
\title{Multi-dimensional Newton Method for Ventilation Rate Estimation}
\description{
  This function implements a multi-dimensional Newton method for estimating ventilation rates in indoor environments based on CO2 time-series data.
}
\usage{
  \code{multi_Newton(persondata, max.iter, tol, verbose, record.steps, critpoints, nadj\_CO2rate, CO2, init.Q, envCO2known, E.init, envCO2.init, freq, temp)}
}
\arguments{
  \item{persondata}{A data frame containing information about individuals in the room, used for Newton-Raphson optimization.}
  \item{max.iter}{Maximum number of iterations for optimization.}
  \item{tol}{Tolerance for convergence in optimization.}
  \item{verbose}{If TRUE, print iteration details during optimization.}
  \item{record.steps}{If TRUE, record optimization steps.}
  \item{critpoints}{Critical points for multi-dimensional Newton method.}
  \item{nadj\_CO2rate}{Numeric vector representing the age-adjusted CO2 generation rate.}
  \item{CO2}{Numeric vector representing the CO2 concentrations over time.}
  \item{init.Q}{Initial estimate of ventilation rate.}
  \item{envCO2known}{Known environmental CO2 concentration in ug/m\^3.}
  \item{E.init}{Initial estimate of emissions rate for multi-dimensional Newton method.}
  \item{envCO2.init}{Initial estimate of environmental CO2 for multi-dimensional Newton method.}
  \item{freq}{The frequency of CO2 measurements (e.g., measurements per hour).}
  \item{temp}{The temperature in degrees Celsius (default is 25).}
}
\value{
  A list containing the estimated ventilation rate, emissions rates, environmental CO2 concentration, number of iterations, and convergence status.
}
\examples{
  \dontrun{
    multi_Newton(persondata = data.frame(time = c(0, 1, 2), activity = c(1, 2, 1)),
                 max.iter = 1000,
                 tol = 1e-8,
                 verbose = TRUE,
                 record.steps = TRUE,
                 critpoints = c(10, 20, 30),
                 nadj\_CO2rate = c(0.5, 0.8, 1.2),
                 CO2 = c(400, 600, 800),
                 init.Q = 1,
                 envCO2known = 500,
                 E.init = c(0.1, 0.2, 0.3),
                 envCO2.init = 450,
                 freq = 1,
                 temp = 25)
  }
}

\keyword{ventilation}
\keyword{CO2}
\keyword{indoor environment}
\keyword{Newton method}
\keyword{optimization}
