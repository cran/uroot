\name{CH.test}
 \alias{CH.test}
  \title{Canova-Hansen test}

  \description{
     Compute the Canova-Hansen test statistic for examining the null hypothesis
of no unit roots at seasonal frequencies against the alternative of a unit
root at either a single or a set of seasonal frequencies.
   }

   \usage{
      CH.test (label, frec, f0, DetTr, showcat)
   }

   \arguments{
       \item{label}{a list object with information about the series. The list consist of the following elements. \code{vari}, a ts object with the data of the series, \code{s}, the periodicity of the series (4, 12 or 1 for quarterly, monthly or anual data), \code{t0}, a vector of two elements indicating the year and season of the first observation, \code{N}, number of observations.}
       \item{frec}{a vector indicating the frequencies to analyse.}
       \item{f0}{a 0-1 (No-Yes) vector of length one indicating wether a lag of the dependent variable is included or not in the auxiliar regression.}
       \item{DetTr}{a logical argument. If TRUE a linear trend is included in the auxiliar regression.}
       \item{showcat}{how the results are showed and stored. If TRUE a brief explanation reports the results. If FALSE a list object stores the results.}
   }

   \Details{
      Elements of \code{frec} contain a 0 if the season assigned to this element
is not considered and contains a 1 for the frequencies that must be analysed. The position of each grequency in the vector is as follows: c(pi/2, pi) for quarterly series and c(pi/6, 5pi/6, pi/3, 2pi/3, pi/2, pi) for monthly series.
   }

   \seealso{
       \code{\link{CHseas.test}}.
   }

   \value{CH test statistic.}

   \author{Javier López-de-Lacalle \email{jedlobej@bs.ehu.es} and Ignacio Díaz-Emparanza \email{etpdihei@bs.ehu.es}.}

   \references{F. Canova and B.E. Hansen (1995), Are seasonal patterns constant over time? A test for seasonal stability. \emph{Journal of Business and Economic Statistics}, \bold{13}, 237-252.
   }

   \examples{
     ## Analyse frequency pi, including a first order lag and but not a linear trend.
     data(AirPassengers)
     AirP <- list(vari=AirPassengers, s=12, t0=c(1949, 1), N=length(AirPassengers))
     CH.test(label=AirP, frec=c(0,0,0,0,0,1), f0=1, DetTr=FALSE, showcat=TRUE)
   }
   \keyword{ts}
