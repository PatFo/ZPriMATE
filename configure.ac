# Process this file with autoconf to produce a configure script.

AC_INIT([cscan], [0.1], [foldenauer@thphys.uni-heidelberg.de])
AM_INIT_AUTOMAKE([-Wall -Werror foreign])

AC_PROG_CXX
AC_PROG_CC
## Needed to build cubature convenience lib:
AC_PROG_RANLIB  

AC_CONFIG_HEADERS([config.h])
AC_CONFIG_FILES([
 Makefile
 src/Makefile
 src/cubature/Makefile
])
AC_CONFIG_SUBDIRS([Cuba-4.2])

AC_OUTPUT