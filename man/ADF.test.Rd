\name{ADF.test}
 \alias{ADF.test}
  \title{Augmented Dickey-Fuller test}

  \description{
     Compute the augmented Dickey-Fuller test statistic for examining the null hypothesis
that exists a unit root at the zero frequency.
   }

   \usage{
     ADF.test (label, compdet, selecP, Mvfic, VFEp, showcat)
   }

   \arguments{
       \item{label}{a list object with information about the series. The list consist of the following elements. \code{vari}, a ts object with the data of the series, \code{s}, the periodicity of the series (4, 12 or 1 for quarterly, monthly or anual data)}, \code{t0}, a vector of two elements indicating the year and season of the first observation, \code{N}, number of observations.
       \item{compdet}{a vector indicating the deterministic components included in the auxiliar regression.}
       \item{selecP}{method for selecting lags included in the auxiliar regression.}
       \item{Mvfic}{a matrix containing any particular dummy.}
       \item{VFEp}{a matrix containing partial seasonal dummy.}
       \item{showcat}{how the results are showed and stored. If TRUE a brief explanation reports the results. If FALSE a list object stores the results.}
   }

   \details{
      \code{compdet} is a numeric vector of length three. If an element equals 1 it indicates that a constant, linear trend or seasonal dummys respectively are included. Otherwise, the element must be set to zero.

      \code{Mvfic} and \code{VFEp} are two kinds of dummys. The first one is a generic dummy,
while the second includes seasonal dummys for only some, not necessarily all of the, seasons.
Any 0-1 column-matrix can be designed.

      Lag selection method availables are the following: "aiclb" is based on the AIC criterion
and Ljung-Box test, "biclb" is based on the BIC criterion and Ljung-Box test, "aiclut" follows a top-down strategy based on the AIC criterion, "biclut" follows a top-down strategy based on the BIC criterion, "signf" retains the significant lags, a numeric vector includes the lags reported on it.
   }

%   \seealso{
%      \code{\link{selecPabic}, \link{slecPv2}, \link{selecPv3}, \link{selecPv4}, \link{selecPv5},
%             \link{selecPv6}, \link{vfic}, \link{VFEp}}.
%   }

   \value{ADF test statistic, coefficient and t-statistics of the deterministic components, and number of available observations.
   }

   \author{Javier López-de-Lacalle \email{jedlobej@bs.ehu.es} and Ignacio Díaz-Emparanza \email{etpdihei@bs.ehu.es}.}

   \references{D.A. Dickey and W.A. Fuller (1981), Likelihood ratio statistics for autoregressive time series with a unit root. \emph{Econometrica}, \bold{49}, 1057-1071.}

   \examples{

     ## ADF.test with constant and seasonal dummys and without other dummys.
     data(AirPassengers)
     AirP <- list(vari=AirPassengers, s=12, t0=c(1949,1), N=length(AirPassengers))
     ADF.test(label=AirP, compdet=c(1,0,1), selecP="biclb", Mvfic=0, VFEp=0, showcat=TRUE)
   }
   \keyword{ts}


