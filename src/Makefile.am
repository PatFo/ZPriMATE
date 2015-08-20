CUBADIR = $(top_srcdir)/Cuba-4.2
BUILDDIR = build
CUBATUREDIR = cubature


AM_CPPFLAGS = -I$(top_srcdir)/mstw -I$(CUBADIR)
# AM_LDFLAGS =  @BOOST_FILESYSTEM_LIB@
# AM_LDFLAGS +=  @BOOST_SYSTEM_LIB@

noinst_PROGRAMS = core
core_SOURCES = main.cpp \
	input.cpp  input.h \
	model.cpp  model.h \
	pheno.cpp  pheno.h \
	spectrum_analysis.cpp  spectrum_analysis.h \
	xsec.cpp  xsec.h \
	../mstw/mstwpdf.cc ../mstw/mstwpdf.h \
	../Cuba-4.2/cuba.h
core_LDADD = ../Cuba-4.2/libcuba.a 

# This is redundand libcuba should be built previous to src
# in case this isn't built correctly uncomment these lines
#$(CUBADIR)/libcuba.a:
#	cd $(CUBADIR) && $(MAKE) $(AM_MAKEFLAGS) libcuba.a
