\name{quarterg}
 \alias{quarterg}
  \title{Buys-Ballot quarterly plot}

  \description{
     Plot of quarterly paths.
   }

   \usage{
      quarterg (vari, s, t0, plot)
   }

   \arguments{
       \item{vari}{a ts object.}
       \item{s}{periodicity of the series.}
       \item{t0}{a vector of two elements indicating the year and season of the first observation.}
       \item{plot}{logical. If TRUE a plot is showed. Otherwise, quarterly data are stored.}
   }

   \seealso{
      \code{\link{bbmp}, \link{bbap}, \link{bb3D}, \link{bbcn}}.
   }

   \value{
      quarterly path data.

      quarterly path plot if it is asked for.
   }

   \author{Javier López-de-Lacalle \email{jedlobej@bs.ehu.es} and Ignacio Díaz-Emparanza \email{etpdihei@bs.ehu.es}.}

   \examples{
      ## Plot in a single graph the path of the four quarters and
      ## store its data in the object \code{MR}.
      data(JohnsonJohnson)
      MR <- quarterg (vari=JohnsonJohnson, s=4, t0=c(1960,1), plot=TRUE)
   }
   \keyword{hplot}
