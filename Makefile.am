#Include the autoconf macros in the m4 directory
ACLOCAL_AMFLAGS = -I m4

# The '.' is necessary as a subdir
# this way .NOTPARALLEL actually works and src is only built AFTER src
SUBDIRS = . src  
dist_doc_DATA = README.md

#Specify path to Cuba
CUBADIR = Cuba-4.2
binary = zprimate



if MAC_SYS
.IGNORE: config.h.in missing aclocal.m4 configure
endif

.NOTPARALLEL: all-recursive

all: all-local

#From the Cuba package build only the libcub.a; rest is not necessary
all-local:
	cd $(CUBADIR) && $(MAKE) $(AM_MAKEFLAGS) libcuba.a
	mkdir -p results

# Create a symbolic link to the executable in the srcdir
# Check whether the bindir already is in the SYSTEMS PATH variable
 #If the it is not present add to PATH
#Else do nothing
install-exec-local:
	mkdir -p $(prefix) && mkdir -p $(bindir)
	ln -s $(abs_top_srcdir)/$(binary) $(bindir)
	@NOT_PRESENT=`echo ${PATH} | grep -o $(bindir)` ; \
	if test -z "$$NOT_PRESENT" ; then \
		if test -n "$(ON_MAC)"; then \
			echo 'IS_PRESENT=`echo $$PATH | grep -o $(bindir)`\nif test -z "$$IS_PRESENT"\nthen\n\techo "export PATH=\$$PATH:$(bindir)" >> ~/.bash_profile\n\tsource ~/.bash_profile\nelse\n\techo "Path is already set up."\nfi'  > setup.sh \
			&& echo "\n\nRun the command \n\n\tsource setup.sh\n\nto set up the ZPriMATE command line tool.\n\n"; \
		else \
			echo -e 'IS_PRESENT=`echo $$PATH | grep -o $(bindir)`\nif test -z "$$IS_PRESENT"\nthen\n\techo "export PATH=\$$PATH:$(bindir)" >> ~/.bashrc\n\tsource ~/.bashrc\nelse\n\techo "Path is already set up."\nfi'  > setup.sh \
			&& echo -e "\n\nRun the command \n\n\tsource setup.sh\n\nto set up the ZPriMATE command line tool.\n\n"; fi;\
	else \
		echo "$(bindir) already in PATH"; fi;

#Remove the installed symbolic link
uninstall-local:
	rm $(bindir)/$(binary)

clean-local:
	-rm -f $(CUBADIR)/libcuba.a
	-rm -f src/*.o
	-rm -f src/*~
	-rm -f *~
