
# MTD-code makefile

CC = gcc -O2 -fopenmp
CFLAGS = -c -Wall
# add -Wall to see most warning messages
SOURCES= train.c predict.c
OBJECTS=$(SOURCES:.c=.o)
EXECUTABLE= migma 
#train predict


.PHONY: all clean wipe depend

all: $(SOURCES) $(EXECUTABLE)

migma: migma.o
	$(CC) migma.o -o migma -lm -lz 

#train: train.o			
#	$(CC) train.o -o train -lm -lz 
#predict: predict.o			
#	$(CC) predict.o -o predict -lm -lz 

.c.o:
	$(CC) $(CFLAGS) $< 
	
	
clean:
	-\rm -f *~
	-\rm -f *.o
	-\rm -f *.d
	-\rm -f core

wipe: clean
	-\rm -f migma 

	#train predict
	





	
