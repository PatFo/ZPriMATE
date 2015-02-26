# cscan
Calculate PDF convoluted crossection from general Z' model.


DEPENDENCIES:
----------------------------

1)MSTW2008:

You need to get the MSTW2008 package from https://mstwpdf.hepforge.org/code/code.html. 
Untar the source files into a local directory and then untar the "Grids" directory in that same location.
Finally you have to set the variable MSTWDIR in the Makefile to point to your local MSTW source directory:

  MSTWDIR = /path/to/local/MSTW/directory

  
  
2)GSL - GNU Scientific Library:

If your system does not have the GNU Scientific Library package you can get it from https://www.gnu.org/software/gsl/





USAGE:
----------------------------

Hadron Cross Section:

Currently there are two different integration strategies implemented, which have to be specified when calcuating 
a hadronic cross section via the int_strategy key:

  1 = GSL QAG adaptive integration (cf https://www.gnu.org/software/gsl/manual/html_node/QAG-adaptive-integration.html#QAG-adaptive-integration)
  2 = GSL Vegas monte carlo (cf https://www.gnu.org/software/gsl/manual/html_node/VEGAS.html)