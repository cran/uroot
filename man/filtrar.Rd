\name{filtrar}
 \alias{filtrar}
  \title{Filter frequencies}

  \description{
      Filter selected frequencies in seasonal time series.
   }

   \usage{
      filtrar (vari, s, t0, filtra, plot)
   }

   \arguments{
       \item{vari}{a ts object.}
       \item{s}{periodicity of the time series.}
       \item{t0}{year and season of the first observation.}
       \item{filtra}{a vector indicating the frequencies to filter.}
       \item{plot}{logical. If TRUE a plot of the filtered series is showed.}
   }

   \details{
      \code{filtra} is a 0-1 vector of length s/2+1. The frequencies are ordered in the vector as follows: c (zero, pi, pi/2) for quarterly series and c (zero, pi, pi/2, 2pi/3, pi/3, 5pi/6, pi/6) for monthly series.

A zero in the element "i" of the vector indicates that the frequency associated to this element must not be filtered. The frequency related to elements equals 1 are filtered.
   }

   \value{ A list containing: Data of the filtered series, coefficients and lags of the polynomial that configure the filter. }

   \author{Javier López-de-Lacalle \email{javlacalle@yahoo.es} and Ignacio Díaz-Emparanza \email{etpdihei@bs.ehu.es}.}

   \examples{
      ## Filter seasonal frequencies.
      data(AirPassengers)
      filtrar (vari=AirPassengers, s=12, t0 = c(1949,1), filtra = c(0,1,1,1,1,1) , plot = TRUE)
   }
   \keyword{hplot}
