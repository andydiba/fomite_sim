CC = gcc
CFLAGS = -Wall
TARGET = fomite
OBJFILES = fomite.o queue.o probability.o

$(TARGET): $(OBJFILES) 
	$(CC) $(CFLAGS) -o $(TARGET) $(OBJFILES) -lm

fomite.o: fomite.c queue.h
	$(CC) $(CFLAGS) -c fomite.c 

queue.o: queue.c queue.h
	$(CC) $(CFLAGS) -c queue.c

probability.o: probability.c probability.h
	$(CC) $(CFLAGS) -c probability.c 

clean:
	rm $(TARGET) $(OBJFILES)

