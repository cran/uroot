# Process this file with autoconf to produce a configure script.
# This file is based on 
# https://cran.r-project.org/doc/manuals/r-devel/R-exts.html#Configure-and-cleanup
# configure.ac from packages 'gpuRcuda' and 'cudabayesreg' and from R sources

AC_INIT([uroot], [2.0-5])
AC_CONFIG_AUX_DIR([tools])
AM_INIT_AUTOMAKE

AC_CONFIG_MACRO_DIR([tools/m4])

: ${R_HOME=`R RHOME`}
if test -z "${R_HOME}"; then
  echo "could not determine R_HOME"
  exit 1
fi
CXX=`${R_HOME}/bin/R CMD config CXX`

AC_MSG_CHECKING("Checking environment variable CUDA_HOME")
if test -z "${CUDA_HOME}"; then
    CUDA_HOME=`find /usr/local/ -maxdepth 1 -type d -name "cuda*" | sort -V | tail -1`
    AC_MSG_RESULT("CUDA_HOME not set; using highest version found ${CUDA_HOME}")
else
    AC_MSG_RESULT("using CUDA_HOME=${CUDA_HOME}")
fi

AC_MSG_CHECKING("Checking if CUDA_HOME was found or set")
HAS_CUDA=TRUE
if test -z "${CUDA_HOME}"; then
  HAS_CUDA=FALSE
fi
AC_MSG_RESULT("HAS_CUDA is ${HAS_CUDA}")

# path to R include files, ${R_HOME}/include may not be correct;
# e.g., in my system it leads to /usr/local/lib/R/include, which does not exist,
# this way returns the existing path /usr/share/R/include
AC_MSG_CHECKING("Checking R include path")
R_INCL=`"${R_HOME}/bin/R" CMD config --cppflags`
AC_MSG_RESULT("R_INCL is ${R_INCL}")

# Libtool
#LT_INIT

AC_PROG_CXX
AC_PROG_INSTALL
AC_SUBST(R_HOME)
AC_SUBST(CXX)
AC_SUBST(CUDA_HOME)
AC_SUBST(HAS_CUDA)
AC_SUBST(R_INCL)
AC_CONFIG_FILES([src/Makevars R/hegy-test.R])
AC_OUTPUT
