\name{urootgui}
 \alias{urootgui}
  \title{uroot R-GUI}

  \description{
     Command for launching a graphical user interface for the package uroot.
   }

   \usage{
     urootgui ()
   }

   \details{
    \code{tcltk}  and \code{BWidget} are needed. The "BWidget" is not included in the minimal Tcl/TK package of R. You can download it at \url{http://www.activestate.com/Products/ActiveTcl/} or at maintainer's homepage (see DESCRIPTION file).

    Once \code{BWidget} has been installed a \eqn{TCL_LIBRARY} environment variable must be set to the Tcl library path, /usr/lib/tcl; C:\\Tcl\\lib\\tcl8.4, and to a non-empty string, e.g. "Yes", a \eqn{MY_TCLTK} environment variable. If after this there are still problems try loading BWidget from an R workspace: \code{tclRequire("BWidget")}, after changing the working directory to the folder where BWidget was extracted, it is suggested /usr/lib/BWidget, C:\\Program Files\\Tcl\\lib.
   }

   \author{Javier López-de-Lacalle \email{jedlobej@bs.ehu.es} and Ignacio Díaz-Emparanza \email{etpdihei@bs.ehu.es}.}

   \examples{
     ## library(uroot)
     ## urootgui()
   }
   \keyword{misc}


