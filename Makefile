CC := gcc

CFLAGS := -O3

LIBS := -lm

OBJECTS := main.o evolve.o potential.o energy.o eos.o \
	transport.o initial.o output.o parameters.o

BINARY := run/hydro

FILES := $(OBJECTS) $(BINARY)

all: hydro

hydro: $(OBJECTS)
	$(CC) $(CFLAGS) $(OBJECTS) $(LIBS) -o $(BINARY)

clean:
	rm -f $(FILES) *\~ *\#
