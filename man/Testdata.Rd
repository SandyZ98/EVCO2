\name{Testdata}
\alias{Testdata}
\docType{data}
\title{
Real Sensor data
}
\description{

Aranet sensor data and Room Dimensions
}
\usage{data("Testdata")}
\format{
  A data frame with 1652 observations on the following 5 variables.
    \Time{\code{Time(dd/mm/yyyy)}}{a numeric vector}
    \Carbondioxide{\code{Carbon dioxide(ppm)}}{a numeric vector}
    \Temperature{\code{Temperature(°C)}}{a numeric vector}
    \humidity{\code{Relative humidity(%)}}{a numeric vector}
    \pressure{\code{Atmospheric pressure(hPa)}}{a numeric vector}

}

\examples{
data(Testdata)
## maybe str(Testdata) ; plot(Testdata) ...
}
\keyword{datasets}
