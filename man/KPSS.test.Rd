\name{KPSS.test}
 \alias{KPSS.test}
  \title{Kwiatkowski-Phillips-Schmidt-Shin test}

  \description{
     This function computes the Kwiatkowski-Phillips-Schmidt-Shin test statistic for examining the null hypothesis
that a given series is level-stationary or stationary around a deterministic trend against the alternative that the series is first difference stationary.
   }
   \usage{
      KPSS.test (vari, l, showcat)
   }

   \arguments{
       \item{vari}{a ts object.}
       \item{l}{lag truncation parameter.}
       \item{showcat}{how the results are showed and stored. If TRUE a brief explanation reports the results. If FALSE a list object stores the results.
}
   }

   \details{
     Lag truncation parameter indicates the number of autocovariances considered different from zero for estimating the variance of the residuals. According to the source paper cited below, the lag truncation parameter may be chosen either as \eqn{integer[4(T/100)^{1/4}]} or \eqn{integer[12(T/100)^{1/4}]}, as well as \code{l} =0.

     Rejection of the null hypothesis implies that the long term frequency contains a unit root.
   }

   \value{Statistic for the level stationary null hypothesis, and statistic for the trend stationary null hypothesis.
   }

   \author{Javier López-de-Lacalle \email{javlacalle@yahoo.es} and Ignacio Díaz-Emparanza \email{etpdihei@bs.ehu.es}}

   \references{
      D. Kwiatkowski, P.C.B. Phillips, P. Schmidt and Y. Shin (1992), Testing the null hypothesis of stationarity against the alternative of a unit root: How sure are we that economic time series have a unit root? \emph{Journal of Econometrics}, \bold{54}, 159-178.
   }

   \examples{
     data(AirPassengers)
     KPSS.test(AirPassengers, l=4, showcat=TRUE)
   }
   \keyword{ts}
