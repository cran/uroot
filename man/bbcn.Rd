\name{bbcn}
 \alias{bbcn}
  \title{Buys-Ballot contour}

  \description{
     Creates a grid of rectangles with colors corresponding to the values in
\code{MR} and add contour lines to the existing plot.
   }

   \usage{
      bbcn (MR, s, t0, color)
   }

   \arguments{
      \item{MR}{a matrix containing the values to be plotted (NAs are allowed).}
      \item{s}{periodicity of the series}
      \item{t0}{year and season of the first observation.}
      \item{color}{logical. If TRUE a coloured graphic is showed. Otherwise, grey scale.}
   }

   \details{
       \code{MR} is a matrix with the data of the series like the matrix provided by
\link{quarterg} in quarterly series and \link{bbmp} in monthly series.
   }

   \seealso{
      \code{\link{persp}, \link{contour}, \link{image}, \link{bb3D}}.
   }

   \value{ Perspective plot of surface over the month-year plane. }

   \author{Javier López-de-Lacalle \email{javlacalle@yahoo.es} and Ignacio Díaz-Emparanza \email{etpdihei@bs.ehu.es}.}

   \examples{
      data(AirPassengers)
      MR <- bbmp (vari=AirPassengers, s=12, t0=c(1949,1), mp=c(1:12), vers="Prot", plot=FALSE)
      bbcn (MR, s=12, t0=c(1949,1), color=TRUE)
   }
   \keyword{hplot}
