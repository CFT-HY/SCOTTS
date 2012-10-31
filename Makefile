#### CC ####

# Louhi
# CC := cc -DMPI -DFFT

# Everything turned on
# CC := mpicc -DMPI -DFFT # -DSILO # -DEXPANSION # -DSILO # -DDUMPFFT

# Example for serial profiling
# CC := gcc -O3 -DPAPI -lpapi

# For generating VampirTrace logs
# CC := vtcc -vt:mpi -DMPI

# Cosmos
CC := icc -DMPI -DFFT -DSILO


#### CFLAGS ####

# on louhi
# CFLAGS := -O3 -L/home/u1/weir/local/lib -I/home/u1/weir/local/include

# on vuori
# CFLAGS := -O3 -L/home/u1/weir/local/lib -I/home/u1/weir/local/include

# on pc168, need libraries
# CFLAGS := -O3 -L/home/weir/local/lib -I/home/weir/local/include

# on afo
# CFLAGS := -O3 -I/home/weir/local/include/

# on my mac
# CFLAGS := -O0 -I/Users/weir/Installed/silo-4.8-bsd/include/

# on Alcyone:
# CFLAGS := -O3 -ipo -m64 -align -fno-alias -fno-fnalias -falign-functions -unroll-aggressive -fp-model fast=2 -I/home/weir/local/include/

# on Cosmos:
CFLAGS := -g -O3 -xHost -align -ansi-alias -mcmodel=medium

#### LIBS ####

# on louhi
# LIBS := -lfftw3_mpi -lfftw3

# on vuori
# LIBS := -lfftw3_mpi -lfftw3 -lm

# on pc168:
# LIBS :=   -lfftw3_mpi -lfftw3 -lm -lsiloh5

# on afo
# LIBS := -L/home/weir/local/lib/ -lstdc++ -lhdf5 -lsiloh5 -lfftw3_mpi -lfftw3 -lm

# on my mac:
# LIBS := -L/Users/weir/Installed/silo-4.8-bsd/lib/ -lsilo -lfftw3_mpi -lfftw3 -lm

# on Alcyone
# LIBS := -L/home/weir/local/lib/ -lfftw3_mpi -lfftw3

# on Cosmos:
LIBS := -lmpi -lfftw3_mpi -lfftw3 -lsiloh5

#### No user-serviceable parts below ####

OBJECTS := main.o evolve.o potential.o energy.o eos.o \
	transport.o initial.o output.o parameters.o silage.o \
	util.o papi.o arrangement.o fft.o alloc.o gw.o velps.o \
	checkpoint.o

BINARY := run/hydro

FILES := $(OBJECTS) $(BINARY)

all: hydro

hydro: $(OBJECTS)
	$(CC) $(CFLAGS) $(OBJECTS) $(LIBS) -o $(BINARY)

clean:
	rm -f $(FILES) *\~ *\#
