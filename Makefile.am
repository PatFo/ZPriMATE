#Include the autoconf macros in the m4 directory
ACLOCAL_AMFLAGS = -I m4

SUBDIRS = src  
dist_doc_DATA = README.md

#Specify path to Cuba
CUBADIR = Cuba-4.2
binary = zprimate

all: all-local

#From the Cuba package build only the libcub.a; rest is not necessary
all-local:
	cd $(CUBADIR) && $(MAKE) $(AM_MAKEFLAGS) libcuba.a
	mkdir -p results
	
# Create a symbolic link to the executable in the srcdir
install-exec-local:
	ln -s $(abs_top_srcdir)/$(binary) $(bindir)

#Remove the installed symbolic link
uninstall-local:
	rm $(bindir)/$(binary)

clean-local:
	-rm -f src/*.o
	-rm -f src/*~
	-rm -r results
	-rm -f *~
