#                                               -*- Autoconf -*-
# netdiffuseR configure.ac
# (with some code borrowed from RcppArmadillo configure.ac
# and ARTP2 configure.ac)
# 
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])

AC_INIT(netdiffuseR, m4_esyscmd_s([awk '/^Version:/ {print $2}' DESCRIPTION]))

## Set R_HOME, respecting an environment variable if one is set 
: ${R_HOME=$(R RHOME)}
if test -z "${R_HOME}"; then
    AC_MSG_ERROR([Could not determine R_HOME.])   
fi
## Use R to set CXX and CXXFLAGS
CXX=$(${R_HOME}/bin/R CMD config CXX)
CXXFLAGS=$("${R_HOME}/bin/R" CMD config CXXFLAGS)

## We are using C++
AC_LANG(C++)
AC_REQUIRE_CPP

## Default the OpenMP flag to the empty string.
## If and only if OpenMP is found, expand to $(SHLIB_OPENMP_CXXFLAGS)
openmp_flag=''

## Check for broken systems produced by a corporation based in Cupertino
AC_MSG_CHECKING([for macOS])
RSysinfoName=$("${R_HOME}/bin/Rscript" --vanilla -e 'cat(Sys.info()[["sysname"]])')
if test x"${RSysinfoName}" == x"Darwin"; then
   AC_MSG_RESULT([found])
   AC_MSG_WARN([OpenMP unavailable and turned off.])
   arma_openmp='-DARMA_DONT_USE_OPENMP'
else
   AC_MSG_RESULT([not found as on ${RSysinfoName}])
   ## Check for OpenMP
   AC_MSG_CHECKING([for OpenMP])
   ## if R has -fopenmp we should be good
   allldflags=$(${R_HOME}/bin/R CMD config --ldflags)
   hasOpenMP=$(echo ${allldflags} | grep -- -fopenmp)
   if test x"${hasOpenMP}" == x""; then
      AC_MSG_RESULT([missing])
      arma_openmp='-DARMA_DONT_USE_OPENMP'
   else   
      AC_MSG_RESULT([found])
      arma_openmp='-DARMA_USE_OPENMP'
      openmp_flag='$(SHLIB_OPENMP_CXXFLAGS)'
   fi
fi

AC_SUBST(openmp_flag)
AC_SUBST(arma_openmp)
AC_CONFIG_FILES([src/Makevars])
AC_OUTPUT

