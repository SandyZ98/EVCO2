\name{secant_method}
\alias{secant_method}
\title{Parameter Estimation using the Secant Method}
\description{
  Estimates a single parameter value by finding the root of the build_up() function using the secant method. This function is specifically designed for scenarios involving air quality analysis, where parameters such as air exchange rate, number of individuals, CO2 generation rate, and room volume play a crucial role.
}
\usage{
secant_method(aer = NA, n = NA, Gp = NA, V = NA, C0, C1, delta_t, Cr = 400, tol = 1e-10, max_iter = 1000)
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
  \item{tol}{Absolute difference of function values between steps used as a stopping condition.}
  \item{max_iter}{Maximum number of iterations.}
}
\details{
  The function first checks that only one of the parameters (aer, n, Gp, V) is set for estimation. It then uses the secant_starter function to find initial values for the iteration process. The secant method is applied iteratively to estimate the parameter by minimizing the build_up() function value.
}
\value{
  The function returns a list with the following components:
  \item{root}{The estimated parameter value at which build_up() is close to zero.}
  \item{f_root}{The value of build_up() at the estimated parameter value.}
  \item{iter}{The number of iterations taken to reach the solution.}
  \item{convergence}{Indicator of successful convergence: 0 if root was found, 1 otherwise.}
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
  \code{\link{secant_starter}}, \code{\link{another_relevant_function}}
}
\examples{
## Example of using secant_method
## Not run:
# secant_method(aer = 0.5, n = NA, Gp = 0.08, V = 50, C0 = 400, C1 = 800, delta_t = 1, tol = 1e-10, max_iter = 1000)
## End(Not run)
}
\keyword{math}
\keyword{optimization}
