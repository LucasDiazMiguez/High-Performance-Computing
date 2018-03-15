include settings.inc

CC			= gcc
CFLAGS			= -std=gnu99 -Wall -Wextra
LDFLAGS			= -lm -lgomp

# Simulation Parameters
PARAMETERS		= -DQ=$(Q) -DL=$(L) -DSAMPLES=$(SAMPLES) \
			  -DTEMP_MIN=$(TEMP_MIN) -DTEMP_MAX=$(TEMP_MAX) -DDELTA_TEMP=$(DELTA_TEMP) \
			  -DTRAN=$(TRAN) -DTMAX=$(TMAX) -DDELTA_T=$(DELTA_T)
CPPFLAGS = $(PARAMETERS)

# Files
TARGETS		= tiny_ising
SOURCES		= $(shell echo *.c)
OBJS		= $(patsubst %.c, %.o, $(C_SOURCES))


# Rules
all: $(TARGETS)

tiny_ising: tiny_ising.o
	$(CC) $(CFLAGS) -o $@ $< $(LDFLAGS)

clean:
	rm -f $(TARGETS) *.o

.PHONY: clean all

