# Process this file with autoconf to produce a configure script.

AC_INIT([ZPriMATE], [0.1], [foldenauer@thphys.uni-heidelberg.de])
AM_INIT_AUTOMAKE([-Wall -Werror foreign subdir-objects])

#Checking for the compiler
AC_PROG_CC(clang  gcc)
AC_PROG_CXX(clang++  g++)
AC_PROG_CPP

dnl If CXX is Clang, check that it can find and parse C++ standard library
dnl headers.
if test "$CXX" = "clang++" ; then
  AC_MSG_CHECKING([whether clang works])
  AC_LANG_PUSH([C++])
  dnl Note that space between 'include' and '(' is required.  There's a broken
  dnl regex in aclocal that otherwise will think that we call m4's include
  dnl builtin.
  AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[#include <limits>
#if __has_include (<cxxabi.h>)
#include <cxxabi.h>
#endif
#if __has_include (<unwind.h>)
#include <unwind.h>
#endif
]])],
[
  AC_MSG_RESULT([yes])
],
[
  AC_MSG_RESULT([no])
  AC_MSG_ERROR([Selected compiler could not find or parse C++ standard library headers.  Rerun with CC=c-compiler CXX=c++-compiler ./configure ...])
])
  AC_LANG_POP([C++])
fi



AC_CONFIG_HEADERS([config.h])
AC_CONFIG_FILES([
 Makefile
 src/Makefile
])
AC_CONFIG_SUBDIRS([Cuba-4.2])

AC_CONFIG_MACRO_DIRS([m4])


#
# Platform specific setup
#
#############################
# Check for which host we are on and setup a few things
# specifically based on the host
AC_CANONICAL_HOST
case $host_os in
  darwin* )
        # Do something specific for mac
        AC_MSG_RESULT([[We noticed you are using macOS. If you are having any trouble, please consult: http://www.ubuntu.com/][$host_os]]) 
        CXXFLAGS="$CXXFLAGS -stdlib=libc++ -std=c++11"
        AC_DEFINE([MAC_SYS], [1], [Variable state is 1 if platform is Mac Osx])
        AC_SUBST([ON_MAC],[1])
	MAC_SYS=TRUE
        ;;
  linux*)
        # Do something specific for linux
        AC_MSG_RESULT([[We are on linux: ][$host_os]])
	CXXFLAGS="$CXXFLAGS -std=c++11"
        AC_DEFINE([MAC_SYS], [0], [Variable state is 1 if platform is Mac Osx])
        ;;
   *BSD*)
        # Do something specific for BSD
        ;;
    *)
        #Default Case
        AC_MSG_ERROR([Your platform is not currently supported])
        ;;
esac

AM_CONDITIONAL([MAC_SYS],test x$MAC_SYS = xTRUE)



#Check for GSL 
#If the library is found then the tests will define the macros HAVE_LIBGSL, HAVE_LIBGSLCBLAS, HAVE_LIBM 
#and add the options -lgsl -lgslcblas -lm to the variable LIBS.
AC_CHECK_LIB([m],[cos])
AC_CHECK_LIB([gslcblas],[cblas_dgemm])
AC_CHECK_LIB([gsl],[gsl_blas_dgemm])

##Check for boost
#AX_BOOST_BASE([1.33.1])
#AX_BOOST_FILESYSTEM
#AX_BOOST_SYSTEM

#Check for Python
PY_PROG_PYTHON([2.7.0], ,AC_MSG_ERROR(You need Python >= 2.7.0 to run ZPriMATE. Use --with-python and refer to a valid python binary.) )
AM_PATH_PYTHON([2.7.0], ,AC_MSG_ERROR(You need Python >= 2.7.0 to run ZPriMATE. Use --with-python and refer to a valid python binary.) )

#Check for Python modules
PY_PYTHON_IMPORT_CHECK(numpy, [import numpy])
PY_PYTHON_IMPORT_CHECK(scipy,[import scipy])
PY_PYTHON_IMPORT_CHECK(scipy.optimize, [from scipy.optimize import minimize, brentq])
PY_PYTHON_IMPORT_CHECK(scipy.optimize, [from scipy.interpolate import InterpolatedUnivariateSpline])
PY_PYTHON_IMPORT_CHECK(matplotlib, [from matplotlib import pyplot as plt])




AC_OUTPUT
