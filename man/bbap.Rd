\name{bbap}
 \alias{bbap}
  \title{Buys-Ballot anual path plot}

  \description{
     Plot anual paths.
   }

   \usage{
      bbap (vari, s, t0, yearsp)
   }

   \arguments{
       \item{vari}{a ts object.}
       \item{s}{periodicity of the series.}
       \item{t0}{a vector of two elements indicating the year and season of the first observation.}
       \item{yearsp}{vector. Years to plot.}
   }

   \seealso{
      \code{\link{quarterg}, \link{bbmp}, \link{bb3D}, \link{bbcn}}.
   }

   \value{
      anual path plot.
   }

   \author{Javier López-de-Lacalle \email{jedlobej@bs.ehu.es} and Ignacio Díaz-Emparanza \email{etpdihei@bs.ehu.es}.}

   \examples{
     ## Plot the following anual paths: 1949, 1950, 1953, 1957
     data(AirPassengers)
     bbap(vari=AirPassengers, s=12, t0=c(1949,1), yearsp=c(1949, 1950, 1953, 1957))
   }
   \keyword{hplot}
