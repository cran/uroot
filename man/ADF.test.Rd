\name{ADF.test}
 \alias{ADF.test}
  \title{Augmented Dickey-Fuller test}

  \description{
     This function computes the augmented Dickey-Fuller statistic for testing the null hypothesis that exists a unit root at the zero frequency.
   }

   \usage{
     ADF.test (label, compdet, selecP, Mvfic, VFEp, showcat)
   }

   \arguments{
       \item{label}{a list object with information about the series. The list consists of the following elements. \code{vari}, a ts object with the data of the series, \code{s}, the periodicity of the series (4, 12, or 1 for quarterly, monthly, or anual data), \code{t0}, a vector of two elements indicating the year and season of the first observation, \code{N}, number of observations.}
       \item{compdet}{a vector indicating the deterministic components included in the auxiliar regression.}
       \item{selecP}{method for selecting lags included in the auxiliar regression.}
       \item{Mvfic}{a matrix containing any particular dummy.}
       \item{VFEp}{a matrix containing partial seasonal dummies.}
       \item{showcat}{how the results are showed and stored. If TRUE a brief explanation reports the results. If FALSE a list object stores the results.}
   }

   \details{
      \code{compdet} is a numeric vector of length three. If an element is equal to 1 it indicates that a constant, linear trend, or seasonal dummies, respectively are included. Otherwise, the element must be set to zero.

      \code{Mvfic} and \code{VFEp} are two kinds of dummies. The first one is a generic dummy and any 0-1 column-matrix can be designed for it. The second one includes seasonal dummies for only some, not necessarily all of the seasons.

      Available lag selection methods are the following: \code{aiclb} is based on the AIC criterion
and Ljung-Box test, \code{biclb} is based on the BIC criterion and Ljung-Box test, \code{aiclut} follows a top-down strategy based on the AIC criterion, \code{biclut} follows a top-down strategy based on the BIC criterion, \code{signf} retains the significant lags. It is also possible to set the argument \code{selecP} equals to a vector, for instance, \code{SelecP=c(1,3,4)} for those lags to be included in the auxiliar regression.
   }

   \value{ADF test statistic for long run unit root hypothesis, coefficient and t-statistics of the deterministic components, and number of available observations.
   }

   \author{Javier López-de-Lacalle \email{javlacalle@yahoo.es} and Ignacio Díaz-Emparanza \email{etpdihei@bs.ehu.es}}

   \references{
   D.A. Dickey and W.A. Fuller (1981), Likelihood ratio statistics for autoregressive time series with a unit root. \emph{Econometrica}, \bold{49}, 1057-1071.

   W.A. Fuller (1976), Introduction to Statistical Time Series. Jonh Wiley, New York.
   }

   \examples{

     ## ADF.test with constant and seasonal dummies and without other dummies.
     data(AirPassengers)
     AirP <- list(vari=AirPassengers, s=12, t0=c(1949,1), N=length(AirPassengers))
     ADF.test(label=AirP, compdet=c(1,0,1), selecP="biclb", Mvfic=0, VFEp=0, showcat=TRUE)
   }
   \keyword{ts}


