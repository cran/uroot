\name{HEGY.test}
 \alias{HEGY.test}
  \title{Hylleberg-Engle-Granger-Yoo test}

  \description{
     This function computes the Hylleberg-Engle-Granger-Yoo statistic for testing the null hypothesis that exits a unit root at long run or seasonal frequencies.
   }

   \usage{
     HEGY.test (label, compdet, selecP, Mvfic, VFEp, showcat)
   }

   \arguments{
       \item{label}{a list object with information about the series. The list consist of the following elements. \code{vari}, a ts object with the data of the series, \code{s}, the periodicity of the series (4, 12, or 1 for quarterly, monthly, or anual data), \code{t0}, a vector of two elements indicating the year and season of the first observation, \code{N}, number of observations.}
       \item{compdet}{a vector indicating the deterministic components included in the auxiliar regression.}
       \item{selecP}{method for selecting lags included in the auxiliar regression.}
       \item{Mvfic}{a matrix containing any particular dummy.}
       \item{VFEp}{a matrix containing partial seasonal dummies.}
       \item{showcat}{how the results are showed and stored. If TRUE a brief explanation reports
                      the results. If FALSE a list object stores the results.}
   }

   \details{
      \code{compdet} is a numeric vector of length three. If an element is equal to 1 it indicates that a constant, linear trend or seasonal dummies respectively are included. Otherwise, the element must be set to zero.

      \code{Mvfic} and \code{VFEp} are two kinds of dummies. The first one is a generic dummy and any 0-1 column-matrix can be designed for it. The second one includes seasonal dummies for only some, not necessarily all of the seasons.

      Available lag selection methods are the following: \code{aiclb} is based on the AIC criterion
and Ljung-Box test, \code{biclb} is based on the BIC criterion and Ljung-Box test, \code{aiclut} follows a top-down strategy based on the AIC criterion, \code{biclut} follows a top-down strategy based on the BIC criterion, \code{signf} retains the significant lags.
     It is also possible to set the argument \code{selecP} equals to a vector, for instance, \code{SelecP=c(1,3,4)} for those lags to be included in the auxiliar regression. If this option is used from the GUI menu, only four lag orders can be introduced, one in each cell. }

   \value{HEGY t-statistics and F-statistics for unit root hypothesis, coefficient and t-statistics of the deterministic components, number of available observations.
   }

   \author{Javier López-de-Lacalle \email{javlacalle@yahoo.es} and Ignacio Díaz-Emparanza \email{etpdihei@bs.ehu.es}}

   \references{S. Hylleberg, R. Engle, C. Granger and B. Yoo (1990),
Seasonal integration and cointegration. \emph{Journal of Econometrics}, \bold{44}, 215-238.

J. Beaulieu and J. Miron (1993), Seasonal unit roots in aggregate U.S. data.
\emph{Journal of Econometrics}, \bold{54}, 305-328.

P.H. Franses (1990), Testing for seasonal unit roots in monthly data, Technical Report 9032, Econometric Institute.
   }

   \examples{
     ## HEGY test with constant and seasonal dummies and without other dummies.
     data(AirPassengers)
     AirP <- list(vari=AirPassengers, s=12, t0=c(1949, 1), N=length(AirPassengers))
     HEGY.test(label=AirP, compdet=c(1,0,1), selecP="biclb", Mvfic=0, VFEp=0, showcat=TRUE)
   }
   \keyword{ts}


