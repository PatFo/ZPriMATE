# cscan
Calculate PDF convoluted crossection from general Z' model.


DEPENDENCIES:
----------------------------

1)MSTW2008:

You need to get the MSTW2008 package from https://mstwpdf.hepforge.org/code/code.html. 
Untar the source files into a local directory and then untar the "Grids" directory in that same location.
Finally you have to add the local MSTW directory to your include path by adding to your ~/.bashrc the following line:

  export CPLUS_INCLUDE_PATH=$CPLUS_INCLUDE_PATH:/path/to/local/MSTW/directory

  
  
2)GSL - GNU Scientific Library:

If your system does not have the GNU Scientific Library package you can get it from https://www.gnu.org/software/gsl/