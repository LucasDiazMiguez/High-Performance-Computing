# Binary file
BIN			= tiny_ising

# Flags
CFLAGS			= -std=gnu99 -Wall -Wextra
LDFLAGS			= -lm -lgomp

# Default Values
L			= 128
SAMPLES			= 1
TEMP_MIN		= 0.9f
TEMP_MAX		= 1.35f
DELTA_TEMP		= 0.05f
TRAN			= 20
TMAX			= 800
DELTA_T			= 5
OFILE			= [CPU,$(Q),$(L),$(SAMPLES),$(TEMP_MIN),$(TEMP_MAX),$(DELTA_TEMP),$(TRAN),$(TMAX),$(DELTA_T)].dat

# Simulation Parameters
PARAMETERS		= -DQ=$(Q) -DL=$(L) -DSAMPLES=$(SAMPLES) \
			  -DTEMP_MIN=$(TEMP_MIN) -DTEMP_MAX=$(TEMP_MAX) -DDELTA_TEMP=$(DELTA_TEMP) \
			  -DTRAN=$(TRAN) -DTMAX=$(TMAX) -DDELTA_T=$(DELTA_T)

# Compilers
CC			= gcc
LINKER			= gcc

# Files
C_SOURCES		= $(BIN).c
HEADERS			=
C_OBJS			= $(patsubst %.c, %.o, $(C_SOURCES))


# Rules
$(BIN): clean $(C_OBJS) $(HEADERS)
	$(LINKER) -o $(BIN) $(C_OBJS) $(LDFLAGS) $(INCLUDES) $(LIBS)

$(C_OBJS): $(C_SOURCES) $(HEADERS)
	$(CC) -c $(C_SOURCES) $(CFLAGS) $(INCLUDES) $(PARAMETERS)

run: $(BIN)
	./$(BIN) > $(OFILE) &

clean:
	rm -f $(BIN) *.o
