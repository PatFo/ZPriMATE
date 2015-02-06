CC = g++
DEBUG = -g
CFLAGS = -Wall -c $(DEBUG)
LFLAGS = -Wall $(DEBUG) -lm
BUILDDIR = build
SRCDIR = src
SRCEXT = cpp
TARGET = bin/cscan

SOURCES = $(shell find $(SRCDIR) -type f -name *.$(SRCEXT))
OBJECTS = $(patsubst $(SRCDIR)/%,$(BUILDDIR)/%,$(SOURCES:.$(SRCEXT)=.o))


$(TARGET) : $(OBJECTS)
	$(CC) $(LFLAGS) $(OBJS) -o $@

$(BUILDDIR)/%.o : $(SRCDIR)/%.cpp
	@mkdir -p $(BUILDDIR)
	$(CC) $(CFLAGS) -c -o  $@ $<


clean:
	$(RM) -r  $(BUILDDIR) $(SRCDIR)/*~  *~ $(TARGET)

#tar:
#	tar cfv $(EXECUTABLE).tar main.cpp