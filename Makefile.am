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
	sed -i '1s:.*:#!$(PYTHON):' 'zprimate'
	sed -i '1s:.*:#!$(PYTHON):' 'bin/calc_limit'
	sed -i '1s:.*:#!$(PYTHON):' 'bin/plot_signal'

# Create a symbolic link to the executable in the srcdir
# Check whether the bindir already is in the SYSTEMS PATH variable
 #If the it is not present add to PATH
#Else do nothing
install-exec-local:
	mkdir -p $(prefix) && mkdir -p $(bindir)
	@if [ -a $(bindir)/$(binary) ] ; \
	then \
		rm $(bindir)/$(binary) ; \
	fi;
	ln -s $(abs_top_srcdir)/$(binary) $(bindir)
	@NOT_PRESENT=`echo ${PATH} | grep -o $(bindir)` ; \
	if test -z "$$NOT_PRESENT" ; then \
		if test -n "$(ON_MAC)"; then \
			echo 'IS_PRESENT=`echo $$PATH | grep -o $(bindir)`\nif test -z "$$IS_PRESENT"\nthen\n\tcp ~/.bash_profile ~/.bash_profile_backup\n \techo "export PATH=\$$PATH:$(bindir)" >> ~/.bash_profile\n\tsource ~/.bash_profile\nelse\n\techo "Path is already set up."\nfi'  > setup.sh \
			&& echo 'IS_PRESENT=`grep -o "$(bindir)" ~/.bash_profile`\nif test -n "$$IS_PRESENT"\nthen\n\tcp ~/.bash_profile ~/.bash_profile_backup\n\tgrep -v "$(bindir)" ~/.bash_profile > ~/.bash_profile2\n\tmv ~/.bash_profile2 ~/.bash_profile\n\techo "Removed $(bindir) from \$$PATH environement variable\nClose session and restart shell to make changes effective"\nelse\n\techo "Nothing to be done"\nfi' > unset.sh \
			&& chmod a+x unset.sh \
			&& echo "\n\nRun the command \n\n\tsource setup.sh\n\nto set up the ZPriMATE command line tool.\nThis will add the install path to your .bash_profile file.\n\n"; \
		else \
			echo -e 'IS_PRESENT=`echo $$PATH | grep -o $(bindir)`\nif test -z "$$IS_PRESENT"\nthen\n\tcp ~/.bashrc ~/.bashrc_backup\n\techo "export PATH=\$$PATH:$(bindir)" >> ~/.bashrc\n\tsource ~/.bashrc\nelse\n\techo "Path is already set up."\nfi'  > setup.sh \
			&& echo -e 'IS_PRESENT=`grep -o "$(bindir)" ~/.bashrc`\nif test -n "$$IS_PRESENT"\nthen\n\tcp ~/.bashrc ~/.bashrc_backup\n\tgrep -v "$(bindir)" ~/.bashrc > ~/.bashrc2\n\tmv ~/.bashrc2 ~/.bashrc\n\techo "Removed $(bindir) from \$$PATH environement variable\nClose session and restart shell to make changes effective"\nelse\n\techo "Nothing to be done"\nfi' > unset.sh \
			&& chmod a+x unset.sh \
			&& echo -e "\n\nRun the command \n\n\tsource setup.sh\n\nto set up the ZPriMATE command line tool.\nThis will add the install path to your .bashrc file.\n\n"; fi;\
	else \
		echo "$(bindir) already in PATH"; fi;

#Remove the installed symbolic link
uninstall-local:
	@if [ -a $(bindir)/$(binary) ] ; \
	then \
		rm $(bindir)/$(binary) ; \
	fi;
	@if [ -a unset.sh ] ; \
	then \
		echo -e "\n\nTo also remove the ZPriMATE path from your \$$PATH environement variable run\n\n\t./unset.sh\n"; \
	fi;

clean-local:
	-rm -f $(CUBADIR)/libcuba.a
	-rm -f src/*.o
	-rm -f src/*~
	-rm -f *~
