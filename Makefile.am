#Include the autoconf macros in the m4 directory
ACLOCAL_AMFLAGS = -I m4

SUBDIRS = src src/cubature  
dist_doc_DATA = README.md

#Specify path to Cuba
CUBADIR = Cuba-4.2

all: all-local

all-local:
	cd $(CUBADIR) && $(MAKE) $(AM_MAKEFLAGS) libcuba.a

clean-local:
	-rm -f src/*.o
	-rm -f src/*~
	-rm -f *~
