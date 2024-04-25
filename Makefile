CC=gcc
CFLAGS=-Wall -lm

# Numele executabilului
TARGET=reconstruct

# Toate fișierele sursă
SOURCES=main-seq2.c pgm_IO.c
HEADERS=pgm_IO.h

# Obiectele corespunzătoare surselor
OBJECTS=$(SOURCES:.c=.o)

# Regula implicită
all: $(TARGET)

# Regula pentru construirea executabilului
$(TARGET): $(OBJECTS)
	$(CC) $(OBJECTS) -o $@ $(CFLAGS)

# Regula pentru construirea fișierelor obiect
%.o: %.c
	$(CC) -c $< -o $@ $(CFLAGS)

# Regula pentru curățare
clean:
	rm -f $(OBJECTS) $(TARGET)

# Regula pentru rularea programului (opțional)
run: $(TARGET)
	./$(TARGET) image_640x480.pgm output.pgm 1000
