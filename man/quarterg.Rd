\name{quarterg}
 \alias{quarterg}
  \title{Buys-Ballot quarterly plot}

  \description{
     Graphic of quarterly paths.
   }

   \usage{
      quarterg (vari, t0, plot)
   }

   \arguments{
       \item{vari}{a ts object.}
       \item{t0}{a vector of two elements indicating the year and season of the first observation.}
       \item{plot}{logical. If TRUE a graphic is showed. Otherwise, quarterly data are stored.}
   }

   \seealso{
      \code{\link{bbmp}}, \code{\link{bbap}}, \code{\link{bb3D}}, and \code{\link{bbcn}}.
   }

   \value{
      quarterly path data.

      graphic of the quarterly path if it is asked for.
   }

   \author{Javier López-de-Lacalle \email{javlacalle@yahoo.es} and Ignacio Díaz-Emparanza \email{etpdihei@bs.ehu.es}}

   \examples{
      ## Plot in a single graph the path of the four quarters and
      ## store its data in the object \code{MR}.
      data(JohnsonJohnson)
      MR <- quarterg (vari=JohnsonJohnson, t0=c(1960,1), plot=TRUE)
   }
   \keyword{hplot}
