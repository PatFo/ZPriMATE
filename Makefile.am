CUBADIR = Cuba-4.2

SUBDIRS = src src/cubature  
dist_doc_DATA = README.md


all: all-local

all-local:
	cd $(CUBADIR) && $(MAKE) $(AM_MAKEFLAGS) libcuba.a

clean-local:
	-rm -f src/*.o
	-rm -f src/*~
	-rm -f *~
