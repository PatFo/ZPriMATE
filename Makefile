# The variable MSTWDIR must point to the MSTW source directory
MSTWDIR = mstw
CUBADIR = cubature

#Extend Search path to include the dependencies
VPATH = $(MSTWDIR):$(CUBADIR)

MSTW = mstwpdf
CUBA = hcubature

# Macros
CC = g++
RELEASE = -O2
DEBUG = -g
MODE = $(RELEASE)
CFLAGS = -Wall -c $(MODE) -I$(MSTWDIR)  -I$(CUBADIR) 
LFLAGS = -Wall $(MODE) -lm -lgsl -lgslcblas -lboost_system -lboost_filesystem
SRCEXT = cpp

#Directories
BINDIR = bin
BUILDDIR = build
EXDIR = example
SRCDIR = src
TSTDIR = test


#EXECUTABLES
TEST_EXEC   = $(TSTDIR)/test 
TARGET_EXEC = $(BINDIR)/core



# Calculate SOURCES, OBJECTS, TARGET and TEST
#Define the 'main' files
TARGET_SOURCE = $(SRCDIR)/main.cpp
TEST_SOURCE  = $(SRCDIR)/test.cpp

#Get common source files
SOURCES = $(shell find $(SRCDIR) -type f -name *.$(SRCEXT))
COMMON_SOURCES = $(filter-out $(TARGET_SOURCE) $(TEST_SOURCE), $(SOURCES))

#Get Common object files
COMMON_OBJECTS = $(patsubst $(SRCDIR)/%,$(BUILDDIR)/%,$(COMMON_SOURCES:.$(SRCEXT)=.o))
COMMON_OBJECTS := $(COMMON_OBJECTS) $(BUILDDIR)/$(MSTW).o $(BUILDDIR)/$(CUBA).o

#Get executable object files
TARGET_OBJECT = $(patsubst $(SRCDIR)/%,$(BUILDDIR)/%,$(TARGET_SOURCE:.$(SRCEXT)=.o)) 
TEST_OBJECT = $(patsubst $(SRCDIR)/%,$(BUILDDIR)/%,$(TEST_SOURCE:.$(SRCEXT)=.o)) 



#MAKE RULES
.PHONY: all cscan test

all: cscan test

cscan: $(TARGET_EXEC)

test: $(TEST_EXEC)


#Rule to make program executable
$(TARGET_EXEC): $(COMMON_OBJECTS) $(TARGET_OBJECT)
	$(CC) $^ -o $@ $(LFLAGS)

#Rule to make test
$(TEST_EXEC): $(COMMON_OBJECTS) $(TEST_OBJECT)
	$(CC) $^ -o $@ $(LFLAGS)


#Rule to make project object files
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
	$(RM) -r   *~  $(TARGET_EXEC) $(TEST_EXEC)   $(BUILDDIR)  $(BINDIR)/*~  $(TSTDIR)/*~ $(SRCDIR)/*~ $(EXDIR)/*~

#tar:
#	tar cfv $(EXECUTABLE).tar main.cpp