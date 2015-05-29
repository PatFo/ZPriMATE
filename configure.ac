# Process this file with autoconf to produce a configure script.

AC_INIT([cscan], [0.1], [foldenauer@thphys.uni-heidelberg.de])
AM_INIT_AUTOMAKE([-Wall -Werror foreign])

AC_PROG_CXX
AC_PROG_CC

AC_CONFIG_HEADERS([config.h])
AC_CONFIG_FILES([
 Makefile
 src/Makefile
])
AC_CONFIG_SUBDIRS([Cuba-4.2])



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
PY_PROG_PYTHON([2.7.3], ,AC_MSG_ERROR(You need Python >= 2.7.3 to run Cscan. Use --with-python and refer to a valid python binary.) )
AM_PATH_PYTHON([2.7.3], ,AC_MSG_ERROR(You need Python >= 2.7.3 to run Cscan. Use --with-python and refer to a valid python binary.) )

#Check for Python modules
#PY_PYTHON_IMPORT_CHECK(scipy.optimize, [from scipy.optimize import minimize, brentq])
PY_PYTHON_IMPORT_CHECK(numpy, [import numpy])
PY_PYTHON_IMPORT_CHECK(matplotlib, [from matplotlib import pyplot as plt])




AC_OUTPUT