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
       \item{krmg}{the number of observations in the intervals in which the the range and mean are calculated.}
   }

   \details{
     In order to ensure that enough observations are employed for calculating each point,
it may choose a choose a value of \code{krmg} based on the number of observations N, for example: sqrt(N).
   }

   \value{ A graphic of the range against the mean, and the correlation between the range and the mean.}

   \author{Javier López-de-Lacalle \email{javlacalle@yahoo.es} and Ignacio Díaz-Emparanza \email{etpdihei@bs.ehu.es}}

   \examples{
     data(AirPassengers)
     rmg (vari=AirPassengers, sqrt(length(AirPassengers)))
   }
   \keyword{hplot}
