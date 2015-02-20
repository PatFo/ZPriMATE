CC = g++
DEBUG = -g
CFLAGS = -Wall -c $(DEBUG)
LFLAGS = -Wall $(DEBUG) -lm
BUILDDIR = build
SRCDIR = src
SRCEXT = cpp
BINDIR = bin
TARGET = $(BINDIR)/cscan

SOURCES = $(shell find $(SRCDIR) -type f -name *.$(SRCEXT))
OBJECTS = $(patsubst $(SRCDIR)/%,$(BUILDDIR)/%,$(SOURCES:.$(SRCEXT)=.o))



$(TARGET) : $(OBJECTS)
	$(CC) $(LFLAGS) $(OBJECTS) -o $@

$(BUILDDIR)/%.o : $(SRCDIR)/%.cpp
	@mkdir -p $(BUILDDIR)
	$(CC) $(CFLAGS) -c -o  $@ $<


clean:
	$(RM) -r   *~  $(TARGET)  $(BINDIR)/*~   $(BUILDDIR)  $(SRCDIR)/*~ 

#tar:
#	tar cfv $(EXECUTABLE).tar main.cpp