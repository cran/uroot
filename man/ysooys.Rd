\name{ysooys}
 \alias{ysooys}
  \title{Search for the date of an observation and vice versa}

  \description{
     This function provides the year and season of an observation and the position of an observation in a given date.
   }

   \usage{
      ysooys (yso, t0, N, s)
   }

   \arguments{
       \item{yso}{the date or position desired to convert into the date or position. See details.}
       \item{t0}{year and season of the first observation in the observed time series.}
       \item{N}{length of the observed time series.}
       \item{s}{periodicity of the observed time series.}
   }

   \details{
      If it is desired to search for the date of an observation, the argument \code{yso} must be a numeric vector of length
      one with the position in the sample of the observation.

      If it is desired to search for the position of an observation related to a given date (year and season), the argument
      \code{yso} must be a numeric vector of length two containing the year and season of the observation.
   }

   \value{
      A list object that consist of the date or position of an observation, \code{out}, and a matrix that relates each date
      with each position in the sample, \code{index}.
   }


   \author{Javier López-de-Lacalle \email{javlacalle@yahoo.es} and Ignacio Díaz-Emparanza \email{etpdihei@bs.ehu.es}}

   \examples{
     ## Given an observed monthly time series that starts in 1990.1 and is of length 120:
     # Which is the position of the data observed in 1991.11?
     obs <- ysooys(yso=c(1991,11), t0=c(1990,1), N=120, s=12)$out
     # Which is the date of the 23rd observation?
     ys <- ysooys(yso=23, t0=c(1990,1), N=120, s=12)$out
   }
   \keyword{utilities}


