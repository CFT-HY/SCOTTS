#### CC ####

# Everything turned on
CC := mpicc -DMPI -DFFT # -DEXPANSION # -DSILO # -DDUMPFFT

# Example for serial profiling
# CC := gcc -O0

# For generating VampirTrace logs
# CC := vtcc -vt:mpi -DMPI


#### CFLAGS ####

# on vuori
CFLAGS := -O3 -L/home/u1/weir/local/lib -I/home/u1/weir/local/include

# on pc168, need libraries
# CFLAGS := -O3 -L/home/weir/local/lib -I/home/weir/local/include

# on afo
# CFLAGS := -O3 -I/home/weir/Installed/silo-4.8-bsd/include/ -I/home/weir/local/include/

# on my mac
# CFLAGS := -O0 -I/Users/weir/Installed/silo-4.8-bsd/include/

#### LIBS ####

# on vuori
LIBS := -lfftw3_mpi -lfftw3 -lm

# on pc168:
# LIBS :=   -lfftw3_mpi -lfftw3 -lm -lsiloh5

# on afo
# LIBS := -L/home/weir/Installed/silo-4.8-bsd/lib/ -L/home/weir/local/lib/ -lsilo -lfftw3_mpi -lfftw3 -lm

# on my mac:
# LIBS := -L/Users/weir/Installed/silo-4.8-bsd/lib/ -lsilo -lfftw3_mpi -lfftw3 -lm


#### No user-serviceable parts below ####

OBJECTS := main.o evolve.o potential.o energy.o eos.o \
	transport.o initial.o output.o parameters.o silage.o \
	util.o papi.o arrangement.o fft.o alloc.o gw.o

BINARY := run/hydro

FILES := $(OBJECTS) $(BINARY)

all: hydro

hydro: $(OBJECTS)
	$(CC) $(CFLAGS) $(OBJECTS) $(LIBS) -o $(BINARY)

clean:
	rm -f $(FILES) *\~ *\#
