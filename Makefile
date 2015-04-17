# The variable MSTWDIR must point to the MSTW source directory
MSTWDIR = /remote/pi104a/foldenauer/local/MSTW
CUBADIR = /remote/pi104a/foldenauer/local/cubature-1.0.2
#Search path
VPATH = $(MSTWDIR):$(CUBADIR)

MSTW = mstwpdf
CUBA = hcubature

# Macros
CC = g++
RELEASE = -O2
DEBUG = -g
MODE = $(RELEASE)
CFLAGS = -Wall -c $(MODE) -I$(MSTWDIR)  -I$(CUBADIR) 
LFLAGS = -Wall $(MODE) -lm -lgsl -lgslcblas -lcuba
BUILDDIR = build
SRCDIR = src
SRCEXT = cpp
BINDIR = bin
TARGET = $(BINDIR)/cscan


SOURCES = $(shell find $(SRCDIR) -type f -name *.$(SRCEXT))
OBJECTS := $(patsubst $(SRCDIR)/%,$(BUILDDIR)/%,$(SOURCES:.$(SRCEXT)=.o))
OBJECTS := $(OBJECTS) $(BUILDDIR)/$(MSTW).o $(BUILDDIR)/$(CUBA).o


$(TARGET) : $(OBJECTS)
	$(CC) $(OBJECTS) -o $@ $(LFLAGS)

$(BUILDDIR)/%.o : $(SRCDIR)/%.cpp
	@mkdir -p $(BUILDDIR)
	$(CC) $(CFLAGS) -o $@ $< 

# Build the MSTWPDF object file in the local build directory
$(BUILDDIR)/$(MSTW).o : $(MSTW).cc $(MSTW).h
	@mkdir -p $(BUILDDIR)
	$(CC) $(CFLAGS) -o $@  $<
	
	
# Build the CUBATUTRE object file in the local build directory
$(BUILDDIR)/$(CUBA).o : $(CUBA).c  cubature.h
	@mkdir -p $(BUILDDIR)
	$(CC) $(CFLAGS) -o $@  $<


clean:
	$(RM) -r   *~  $(TARGET)  $(BINDIR)/*~   $(BUILDDIR)  $(SRCDIR)/*~ 

#tar:
#	tar cfv $(EXECUTABLE).tar main.cpp