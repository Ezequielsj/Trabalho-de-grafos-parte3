CC = gcc
CFLAGS = -Wall -Wextra -O3 -march=native -g
LDFLAGS = -lm

SRC = main.c biblioteca_grafos.c fila.c minheap.c
OBJ = $(SRC:.c=.o)
TARGET = programa

all: $(TARGET)

$(TARGET): $(OBJ)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

rodar: all
	./$(TARGET) --tp3-todos

clean:
	rm -f $(OBJ) $(TARGET)

.PHONY: all rodar clean
