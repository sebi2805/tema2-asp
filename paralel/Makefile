CC=mpicc

CFLAGS=-O2

TARGET=program

OBJECTS=main.o pgm_IO.o helpers.o

all: $(TARGET)

$(TARGET): $(OBJECTS)
	$(CC) -o $@ $^

main.o: main.c pgm_IO.h helpers.h
	$(CC) $(CFLAGS) -c main.c

pgm_IO.o: pgm_IO.c pgm_IO.h
	$(CC) $(CFLAGS) -c pgm_IO.c

helpers.o: helpers.c helpers.h
	$(CC) $(CFLAGS) -c helpers.c

clean:
	rm -f $(TARGET) *.o

.PHONY: all clean
