\name{bb3D}
 \alias{bb3D}
  \title{Buys-Ballot 3D plot}

  \description{
     Draw perspective plot of surface over the month-year plane.
   }

   \usage{
      bb3D (MR, s, t0, color, x, y)
   }

   \arguments{
       \item{MR}{a matrix containing the values to be plotted (NAs are allowed).}
       \item{s}{periodicity of the series.}
       \item{t0}{a vector of two elements indicating the year and season of the first observation.}
       \item{color}{logical. If TRUE a coloured graphic is showed. Otherwise, grey scale.}
       \item{x, y}{angles defining the viewing direction. \code{x} gives the azimuthal direction and \code{y} the colatitude.}
   }

   \details{
       \code{MR} is a matrix with the data of the series like the matrix provided by
\link{quarterg} in quarterly series and \link{bbmp} in monthly series.

       When the plot is showed a window is launched as well. Coordinates of the mouse inside this window represent angles defining the viewing direction. Clicking on it with the right button of the mouse the plot rotates.
   }

   \seealso{
      \code{\link{persp}, \link{contour}, \link{image} and \link{bbcn}.}
   }

   \value{ Perspective plot of surface over the month-year plane. }

   \author{Javier López-de-Lacalle \email{javlacalle@yahoo.es} and Ignacio Díaz-Emparanza \email{etpdihei@bs.ehu.es}.}

   \examples{
      data(AirPassengers)
      MR <- bbmp (vari=AirPassengers, s=12, t0 = c(1949,1), mp = c(1:12), vers = "Prot", plot = FALSE)
      bb3D (MR, s=12, t0=c(1949,1), color = TRUE, x = 30, y = 30)
   }
   \keyword{hplot}
