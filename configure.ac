# Process this file with autoconf to produce a configure script.
# This file is based on 
# https://cran.r-project.org/doc/manuals/r-devel/R-exts.html#Configure-and-cleanup
# and on configure.ac from packages 'gpuRcuda', 'cudabayesreg' and from R sources

AC_INIT([uroot], [2.0-6])
AC_CONFIG_AUX_DIR([tools])
AM_INIT_AUTOMAKE

AC_CONFIG_MACRO_DIR([tools])

#: ${R_HOME=`R RHOME`}
#if test -z "${R_HOME}"; then
#  echo "could not determine R_HOME"
#  exit 1
#fi

AC_MSG_CHECKING("Checking host/Arina")
#it would be more natural to detect the server using the variable $HOST, which is 
#set to guinness, but when used by the configure script it seems it changes to "cn0",
#which seems less unique for a node name than "guinness", 
#as an alternative, hostid is used
#
#if test -z "${HOST}"; then
#  #HOST is already defined in the Guinness server, but may not be set in other systems,
#  #so define it via hostname
#  HOST=`hostname`
#fi

HOSTID=`hostid`

# if the hostid is "800a3002" and the directory "/software" is present, then it is assumed
# that we are connected to the server Guinness in the cluster Arina
#if test "${HOST}" = "guinness"; then
#if test "${HOST}" = "cn0"; then
if test "${HOSTID}" = "800a3002"; then
  if test -d "/software"; then
    ARINA=TRUE
  else
    ARINA=FALSE
  fi
else
  ARINA=FALSE
fi
AC_MSG_RESULT("using ARINA=${ARINA}")

AC_MSG_CHECKING("Checking environment variable R_HOME")
#if test -d "/software"; then
if test "${ARINA}" = "TRUE"; then
  #in the cluster arina the latest versions of software are under the directory /software:
  #there are multiple versions installed that directory and "R RHOME" does not return the latest one,
  #so set R_HOME here, current result:
  #R_HOME=/software/R-3.2.2/lib64/R
  #
  #Note: the version of "sort" in the cluster arina does not provide the flag -V, --version-sort,
  #which would be covenient to make sure that the latest version is chosen
  R_HOME=`find /software/ -maxdepth 1 -type d -name "R*" | sort | tail -1`/lib64/R
else
  #no search is done here, R_HOME should be defined or created before running this script
  #this also works in arina but sets an older version of R outside /software;  
  #R_HOME=`/usr/local/R-devel/lib/R/bin/R RHOME`
  R_HOME=`R RHOME`
fi
if test -z "${R_HOME}"; then
  # "R RHOME" above may fail, for example if the R binary is not in the search PATH, 
  # in that case RHOME should be set manually as indicated in the follwing message
  echo "Could not determine R_HOME:" 
  echo "you may set R_HOME by typing \"export R_HOME=path_to_R\" on a bash shell," 
  echo "  where \"path_to_R\" is the output output returned by R.home() in an R console."
  exit 1
fi
AC_MSG_RESULT("using R_HOME=${R_HOME}")

AC_MSG_CHECKING("Checking environment variable CUDA_HOME")
#if test -d "/software"; then
if test "${ARINA}" = "TRUE"; then
  #in the cluster arina the latest versions of software are under the directory /software;
  #CUDA_HOME is defined, but it points to a directory outside /software, so overwrite here 
  #with the latest version found under directory /software (flag -V for sort not available)
  #
  #CUDA_HOME=/software/cuda
  CUDA_HOME=`find /software/ -maxdepth 1 -type d -name "cuda*" | sort | tail -1`
  AC_MSG_RESULT("using CUDA_HOME=${CUDA_HOME}")
else
  if test -z "${CUDA_HOME}"; then
    CUDA_HOME=`find /usr/local/ -maxdepth 1 -type d -name "cuda*" | sort -V | tail -1`
    AC_MSG_RESULT("CUDA_HOME not set; using highest version found: CUDA_HOME=${CUDA_HOME}.")
  fi
fi

AC_MSG_CHECKING("Checking environment variable HAS_CUDA")
if test -z "${CUDA_HOME}"; then
  HAS_CUDA=FALSE
  echo "Could not determine CUDA_HOME:"
  echo "If CUDA is available in your system, the environment variable CUDA_HOME shoud be defined"
  echo "otherwise, the package will still be installed without the functionalies based on CUDA"
  #echo "in order to set CUDA_HOME you may, for example, type the following on a bash shell:"
  #echo "export CUDA_HOME=path_to_cuda_in_your_system"
else
  HAS_CUDA=TRUE
fi
AC_MSG_RESULT("using HAS_CUDA=${HAS_CUDA}")

AC_MSG_CHECKING("Checking R include path")
# Set R include path
#if test -d "/software"; then
if test "${ARINA}" = "TRUE"; then
  #in the cluster arina
  #`"${R_HOME}/bin/R" CMD config --cppflags` returns
  #"R was not built as a library" and
  #`"${R_HOME}/bin/R" CMD config CPPFLAGS` returns -I/usr/local/include 
  #which seems wrong path, so set:
  R_INCL=-I${R_HOME}/include
else
  R_INCL=`"${R_HOME}/bin/R" CMD config --cppflags`
fi
AC_MSG_RESULT("using R_INCL=${R_INCL}")

if test "${HAS_CUDA}" = "FALSE"; then
  #COMPILE_CUDA_FALSE.mk will be run; an empty 'uroot.so' object will be created
  #https://cran.r-project.org/doc/manuals/r-devel/R-exts.html#Configure-and-cleanup
  #packages written with C++ need to pick up the details for the C++ compiler and switch the current language to C++
  CXX=`"${R_HOME}/bin/R" CMD config CXX`
  #CXX is g++, icpc
  #an empty uroot.so object is created, so flags are not critical
  #CXXFLAGS=`"${R_HOME}/bin/R" CMD config CXXFLAGS`
  AC_LANG(C++)
  
  #if test -d "/software"; then
  if test "${ARINA}" = "TRUE"; then
    #set manually, for some reason the R CMD config LDFLAGS returns a wrong -L path,
    #-L/usr/local/lib64, where libR.so is not available
    LDFLAGS=-L/software/R-2.15.2/lib64/R/lib/
  else
    LDFLAGS=`"${R_HOME}/bin/R" CMD config LDFLAGS`
  fi
fi

#if ${HAS_CUDA} = TRUE no need to set CXX, CXXFLAGS,... since 
#there are no predetermined flags returned by "R CMD config" 
#for the nvcc compiler

# Libtool
#LT_INIT

AC_PROG_CXX
AC_PROG_INSTALL
#AC_SUBST(HOST)
#AC_SUBST(ARINA)
#AC_SUBST(R_HOME)
AC_SUBST(R_INCL)
AC_SUBST(CXX)
AC_SUBST(LDFLAGS)
AC_SUBST(CUDA_HOME)
AC_SUBST(HAS_CUDA)
AC_CONFIG_FILES([src/Makevars R/hegy-test.R])
AC_OUTPUT
