#include "R.h"

/*
this file is compiled if a CUDA enabled GPU is not detected in the system;
in that case uroot.so is generated with the purpose of meeting the requirement 
of the NAMESPACE "useDynLib(uroot)"

Alternatively it can be checked whether "if" statements are possible in 
the NAMESPACE file; if so, NAMESPACE.in could be created containing: 
  if (@HAS_CUDA@)
    useDynLib(uroot)
and then NAMESPACE.in would be added to AC_CONFIG_FILES in configure.ac
*/
