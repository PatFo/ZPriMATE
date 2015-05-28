CUBADIR = Cuba-4.2
CUBA = cuba
MSTWDIR = mstw
MSTW = mstwpdf
BUILDDIR = build
CUBATUREDIR = cubature
CUBATURE = hcubature


AM_CPPFLAGS = -I$(top_srcdir)/$(MSTWDIR) -I$(top_srcdir)/$(CUBADIR)  -I$(CUBATUREDIR)
AM_LDFLAGS =  -lm -lgsl -lgslcblas -lboost_system -lboost_filesystem

noinst_PROGRAMS = core
core_SOURCES = main.cpp \
	input.cpp  input.h \
	model.cpp  model.h \
	pheno.cpp  pheno.h \
	spectrum_analysis.cpp  spectrum_analysis.h \
	xsec.cpp  xsec.h \
	$(top_srcdir)/$(MSTWDIR)/$(MSTW).cc $(top_srcdir)/$(MSTWDIR)/$(MSTW).h \
	$(top_srcdir)/$(CUBADIR)/$(CUBA).h
core_LDADD = $(CUBATUREDIR)/libcubature.a  $(top_srcdir)/$(CUBADIR)/libcuba.a
SUBDIRS = $(CUBATUREDIR)


# # Build the CUBATURE object file in the local build directory
# $(CUBATURE).o : $(CUBATURE).c  $(CUBATUREDIR).h
# 	$(CC) $(CFLAGS) -o $@  $<	
	
# Build the MSTWPDF object file in the local build directory
#$(MSTW).o : $(top_srcdir)/$(MSTWDIR)/$(MSTW).cc  $(top_srcdir)/$(MSTWDIR)/$(MSTW).h
#	$(CC) $(CFLAGS) -o $@  $<