\name{bbmp}
 \alias{bbmp}
  \title{Buys-Ballot monthly plot}

  \description{
     Plot of monthly paths.
   }

   \usage{
      bbmp (vari, s, t0, mp, vers, plot)
   }

   \arguments{
       \item{vari}{a ts object.}
       \item{s}{periodicity of the series.}
       \item{t0}{a vector of two elements indicating the year and season of the first observation.}
       \item{mp}{months to plot.}
       \item{vers}{how the plot is showed.}
       \item{plot}{logical. If TRUE a plot is showed. Otherwise, quarterly data are stored.}
   }

   \details{
      If \code{vers} is the character string "Prot" a panel of four graphics is showed. This option detracted from months selected in the vector \code{mp}.

      If \code{vers} is the character string "R" all the months selected in the vector \code{mp} are showed in the same graphic.
   }

   \seealso{
      \code{\link{quarterg}, \link{bbap}, \link{bb3D}, \link{bbcn}}.
   }

   \value{
      monthly path data.

      monthly path plot if it is asked for.
   }

   \author{Javier López-de-Lacalle \email{jedlobej@bs.ehu.es} and Ignacio Díaz-Emparanza \email{etpdihei@bs.ehu.es}.}

   \examples{
     ## Plot in a single graph the path of January, February, November and December and
     ## store its data in the object \code{jfnd}.
     data(AirPassengers)
     jfnd <- bbmp (AirPassengers, s=12, t0=c(1949, 1), mp=c(1,2,11,12), vers="R", plot=TRUE)
   }
   \keyword{hplot}
