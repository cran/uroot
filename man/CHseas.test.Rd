\name{CHseas.test}
 \alias{CHseas.test}
  \title{Canova-Hansen test for seasonal stability}

  \description{
      This function computes the Canova-Hansen statistic for testing the null hypothesis that seasonal stability exists against the alternative that the analysed season is not regular. This version of the CH test is equivalent to a KPSS test for each seasonal path.
   }

   \usage{
      CHseas.test (label, lmax, seas, showcat)
   }

   \arguments{
       \item{label}{a list object with information about the series. The list consist of the following elements. \code{vari}, a ts object with the data of the series, \code{s}, the periodicity of the series (4, 12, or 1 for quarterly, monthly, or anual data), \code{t0}, a vector of two elements indicating the year and season of the first observation, \code{N}, number of observations.}
       \item{lmax}{lag truncation parameter.}
       \item{seas}{a vector of length one indicating the season to analyse.}
       \item{showcat}{how the results are showed and stored. If TRUE a brief explanation reports the results. If FALSE a list object stores the results.}
   }

   \Details{
      Lag truncation parameter \code{lmax} indicates the number of autocovariances considered different from zero for estimating the variance of the residuals.

      \code{seas} is a numeric object indicating the position of the season to analyse.
   }

   \seealso{
       \code{\link{CH.test}}.
   }

   \value{CH test statistic.}

   \author{Javier López-de-Lacalle \email{javlacalle@yahoo.es} and Ignacio Díaz-Emparanza \email{etpdihei@bs.ehu.es}}

   \references{F. Canova and B.E. Hansen (1995), Are seasonal patterns constant over time? A test for seasonal stability. <em>Journal of Business and Economic Statistics</em>, <b>13</b>, 237-252.

      D. Kwiatkowski, P.C.B. Phillips, P. Schmidt and Y. Shin (1992), Testing the null hypothesis of stationarity against the alternative of a unit root: How sure are we that economic time series have a unit root? \emph{Journal of Econometrics}, \bold{54}, 159-178.
   }

   \examples{
     ## Analyse the stability of the first month.
     data(AirPassengers)
     AirP <- list(vari=AirPassengers, s=12, t0=c(1949, 1), N=length(AirPassengers))
     CHseas.test(label=AirP, lmax=3, seas=1, showcat=TRUE)
   }
   \keyword{ts}
