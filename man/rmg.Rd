\name{rmg}
 \alias{rmg}
  \title{Range mean plot}

  \description{
      Plot of range against mean.
   }
   \usage{
      rmg (vari, krmg)
   }

   \arguments{
       \item{vari}{a ts object.}
       \item{krmg}{the number of observations employed for calculating each point.}
   }

   \details{
     In order to ensure that enough observations are employed for calculating each point,
it may choose a choose a value of <code>krmg</code> based on the number of observations T, for example: sqrt(T).
   }

   \value{ Shows a plot of the range against the mean, correlation between the range and the mean.}

   \author{Javier López-de-Lacalle \email{javlacalle@yahoo.es} and Ignacio Díaz-Emparanza \email{etpdihei@bs.ehu.es}.}

   \examples{
     data(AirPassengers)
     rmg (vari=AirPassengers, sqrt(length(AirPassengers)))
   }
   \keyword{hplot}
