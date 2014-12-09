
# MTD-code makefile

CC = gcc -O2 -fopenmp
CFLAGS = -c -Wall
SOURCES= src/migma.c
OBJECTS=$(SOURCES:.c=.o)
EXECUTABLE= migma 


.PHONY: all clean wipe depend

all: $(SOURCES) $(EXECUTABLE)

migma: src/migma.o
	$(CC) migma.o -o migma -lm -lz 

.c.o:
	$(CC) $(CFLAGS) $< 
	
	
clean:
	-\rm -f *~
	-\rm -f *.o
	-\rm -f *.d
	-\rm -f core

wipe: clean
	-\rm -f migma 
	





	
