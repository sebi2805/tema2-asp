# Definirea compilatorului
CC=mpicc

# Opțiuni de compilare, de exemplu -g pentru debug și -O2 pentru optimizare
CFLAGS=-O2

# Numele executabilului
TARGET=program

# Construiește ținta finală
all: $(TARGET)

# Link-ul final pentru executabil
$(TARGET): main2.o pgm_IO.o
	$(CC) -o $@ $^

# Regula pentru fișierele sursă
main.o: main2.c pgm_IO.h
	$(CC) $(CFLAGS) -c main.c

pgm_IO.o: pgm_IO.c pgm_IO.h
	$(CC) $(CFLAGS) -c pgm_IO.c

# Curățarea proiectului
clean:
	rm -f $(TARGET) *.o

.PHONY: all clean
