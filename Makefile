CC=gcc
CFLAGS=-std=gnu99 -Wall -Wextra
LDFLAGS=-lm -lgomp

# Files
TARGETS=tiny_ising
SOURCES=$(shell echo *.c)
OBJS=$(patsubst %.c, %.o, $(C_SOURCES))

# Rules
all: $(TARGETS)

tiny_ising: tiny_ising.o
	$(CC) $(CFLAGS) -o $@ $< $(LDFLAGS)

clean:
	rm -f $(TARGETS) *.o

.PHONY: clean all
