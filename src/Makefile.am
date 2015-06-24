CUBADIR = Cuba-4.2
MSTWDIR = mstw
BUILDDIR = build
CUBATUREDIR = cubature

AM_LDFLAGS = -std=c++11 -pthread
AM_CPPFLAGS = -I$(top_srcdir)/$(MSTWDIR) -I/Library/Developer/CommandLineTools/usr/lib/c++/v1 -I$(top_srcdir)/$(CUBADIR) -std=c++11 -pthread
# AM_LDFLAGS =  @BOOST_FILESYSTEM_LIB@
# AM_LDFLAGS +=  @BOOST_SYSTEM_LIB@

noinst_PROGRAMS = core
core_SOURCES = main.cpp \
	input.cpp  input.h \
	model.cpp  model.h \
	pheno.cpp  pheno.h \
	spectrum_analysis.cpp  spectrum_analysis.h \
	xsec.cpp  xsec.h \
	$(top_srcdir)/$(MSTWDIR)/mstwpdf.cc $(top_srcdir)/$(MSTWDIR)/mstwpdf.h \
	$(top_srcdir)/$(CUBADIR)/cuba.h
core_LDADD = $(top_srcdir)/$(CUBADIR)/libcuba.a 