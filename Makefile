# CC := gcc -DPAPI # -DSILO
# CC := vtcc -vt:mpi -DMPI # -DFFT
CC := mpicc -DMPI -DFFT -DSILO

# on pc168
CFLAGS := -O3 -L/home/weir/local/lib -I/home/weir/local/include
# on afo
# CFLAGS := -O3 -I/home/weir/Installed/silo-4.8-bsd/include/
# on my mac:
# CFLAGS := -O3 -I/Users/weir/Installed/silo-4.8-bsd/include/

# on pc168:
LIBS :=   -lfftw3_mpi -lfftw3 -lm -lsiloh5 # -lpapi
# on afo
# LIBS := -L/home/weir/Installed/silo-4.8-bsd/lib/ -lsilo
# on my mac:
# LIBS := -lm -L/Users/weir/Installed/silo-4.8-bsd/lib/ -lsilo


OBJECTS := main.o evolve.o potential.o energy.o eos.o \
	transport.o initial.o output.o parameters.o silage.o \
	util.o papi.o arrangement.o fft.o

BINARY := run/hydro

FILES := $(OBJECTS) $(BINARY)

all: hydro

hydro: $(OBJECTS)
	$(CC) $(CFLAGS) $(OBJECTS) $(LIBS) -o $(BINARY)

clean:
	rm -f $(FILES) *\~ *\#
