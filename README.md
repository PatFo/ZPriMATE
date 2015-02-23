# cscan
Calculate PDF convoluted crossection from general Z' model.


DEPENDENCIES:
----------------------------

1)MSTW2008:

You need to get the MSTW2008 package from https://mstwpdf.hepforge.org/code/code.html. 
Untar the source files into a local directory and then untar the "Grids" directory in that same location.
Finally you have to set variable MSTWDIR in the Makefile to point to your local MSTW source directory:

  MSTWDIR = /path/to/local/MSTW/directory

  
  
2)GSL - GNU Scientific Library:

If your system does not have the GNU Scientific Library package you can get it from https://www.gnu.org/software/gsl/