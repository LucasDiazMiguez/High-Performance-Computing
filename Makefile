CC			= gcc
CFLAGS			= -std=gnu99 -Wall -Wextra
LDFLAGS			= -lm -lgomp

# Default Values
L			= 128
SAMPLES			= 1
TEMP_MIN		= 2.1f
TEMP_MAX		= 2.5f
DELTA_TEMP		= 0.05f
TRAN			= 20
TMAX			= 800
DELTA_T			= 5
OFILE			= [CPU,$(Q),$(L),$(SAMPLES),$(TEMP_MIN),$(TEMP_MAX),$(DELTA_TEMP),$(TRAN),$(TMAX),$(DELTA_T)].dat

# Simulation Parameters
PARAMETERS		= -DQ=$(Q) -DL=$(L) -DSAMPLES=$(SAMPLES) \
			  -DTEMP_MIN=$(TEMP_MIN) -DTEMP_MAX=$(TEMP_MAX) -DDELTA_TEMP=$(DELTA_TEMP) \
			  -DTRAN=$(TRAN) -DTMAX=$(TMAX) -DDELTA_T=$(DELTA_T)
CPPFLAGS = $(PARAMETERS)

# Files
TARGETS		= tiny_ising tiny_ising_rb 
SOURCES		= $(shell echo *.c)
OBJS		= $(patsubst %.c, %.o, $(C_SOURCES))


# Rules
all: $(TARGETS)

tiny_ising: tiny_ising.o
	$(CC) $(CFLAGS) -o $@ $< $(LDFLAGS)

tiny_ising_rb: tiny_ising_rb.o
	$(CC) $(CFLAGS) -o $@ $< $(LDFLAGS)

clean:
	rm -f $(TARGETS) *.o

.PHONY: clean all

