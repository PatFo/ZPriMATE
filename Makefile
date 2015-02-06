CC = g++
DEBUG = -g
CFLAGS = -Wall -c $(DEBUG)
LFLAGS = -Wall $(DEBUG) -lm
OBJS = main.o model.o pheno.o

EXECUTABLE = cscan

$(EXECUTABLE) : $(OBJS)
	$(CC) $(LFLAGS) $(OBJS) -o $@

main.o : main.cpp
	$(CC) $(CFLAGS) main.cpp
    
model.o : model.h model.cpp
	$(CC) $(CFLAGS) model.cpp

pheno.o : pheno.h pheno.cpp
	$(CC) $(CFLAGS) pheno.cpp

clean:
	\rm *.o *~ $(EXECUTABLE)

tar:
	tar cfv $(EXECUTABLE).tar main.cpp