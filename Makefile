# The variable MSTWDIR must point to the MSTW source directory
MSTWDIR = /remote/pi104a/foldenauer/local/MSTW
MSTW = mstwpdf

# Macros
CC = g++
RELEASE = -O2
DEBUG = -g
MODE = $(RELEASE)
CFLAGS = -Wall -c $(MODE) -I$(MSTWDIR) 
LFLAGS = -Wall $(MODE) -lm -lgsl -lgslcblas
BUILDDIR = build
SRCDIR = src
SRCEXT = cpp
BINDIR = bin
TARGET = $(BINDIR)/cscan


SOURCES = $(shell find $(SRCDIR) -type f -name *.$(SRCEXT))
OBJECTS := $(patsubst $(SRCDIR)/%,$(BUILDDIR)/%,$(SOURCES:.$(SRCEXT)=.o))
OBJECTS := $(OBJECTS) $(BUILDDIR)/$(MSTW).o


$(TARGET) : $(OBJECTS)
	$(CC) $(OBJECTS) -o $@ $(LFLAGS)

$(BUILDDIR)/%.o : $(SRCDIR)/%.cpp
	@mkdir -p $(BUILDDIR)
	$(CC) $(CFLAGS) -o $@ $< 

# Build the MSTWPDF object file in the local build directory
$(BUILDDIR)/$(MSTW).o : $(MSTWDIR)/$(MSTW).cc  $(MSTWDIR)/$(MSTW).h
	@mkdir -p $(BUILDDIR)
	$(CC) $(CFLAGS) -o $@  $<


clean:
	$(RM) -r   *~  $(TARGET)  $(BINDIR)/*~   $(BUILDDIR)  $(SRCDIR)/*~ 

#tar:
#	tar cfv $(EXECUTABLE).tar main.cpp