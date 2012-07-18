CC := gcc -DPAPI # -DSILO

# on pc168
CFLAGS := -O3
# on my mac:
# CFLAGS := -O3 -I/Users/weir/Installed/silo-4.8-bsd/include/

# on pc168:
LIBS := -lm -lsiloh5 -lpapi
# on my mac:
# LIBS := -lm -L/Users/weir/Installed/silo-4.8-bsd/lib/ -lsilo


OBJECTS := main.o evolve.o potential.o energy.o eos.o \
	transport.o initial.o output.o parameters.o silage.o \
	util.o papi.o

BINARY := run/hydro

FILES := $(OBJECTS) $(BINARY)

all: hydro

hydro: $(OBJECTS)
	$(CC) $(CFLAGS) $(OBJECTS) $(LIBS) -o $(BINARY)

clean:
	rm -f $(FILES) *\~ *\#
