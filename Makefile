#### CC ####

# Louhi
# CC := cc -DMPI -DFFT -DSILO

# Everything turned on
CC := mpicc -g -DMPI -DSILO -DFFT -DINITPS # -DSCALAR # -DFFT # -DEXPANSION

# Example for serial profiling
# CC := gcc -O3 -DPAPI -lpapi

# For generating VampirTrace logs
# CC := vtcc -vt:mpi -DMPI

# Cosmos
# CC := icc -DMPI -DFFT -DSILO

# Sisu
# CC := cc -DMPI -DSILO -DFFT # -DSCALAR

#  Taito
# CC := mpicc -DMPI -DSILO -DFFT # -DINITPS # -DDIVPS

#### CFLAGS ####

# on Louhi
# CFLAGS := -O3 -L/home/u1/weir/local/lib -I/home/u1/weir/local/include

# on Sisu
# CFLAGS := -L/homeappl/home/weir/local_sisu/lib \
#	-I/homeappl/home/weir/local_sisu/include \
#	 -O3 -opt-prefetch -unroll-aggressive -no-prec-div \
#	-fp-model fast=2 -align -fno-alias -fno-fnalias -ipo

# on Taito
# CFLAGS := -L/homeappl/home/weir/local_taito/lib \
#	-I/homeappl/home/weir/local_taito/include \
#	-O3 -opt-prefetch -unroll-aggressive -no-prec-div \
#	-fp-model fast=2 -align -fno-alias -fno-fnalias -ipo


# on vuori
# CFLAGS := -O3 -L/home/u1/weir/local/lib -I/home/u1/weir/local/include

# on pc168, need libraries
# CFLAGS := -O3 -L/home/weir/local/lib -I/home/weir/local/include

# on afo
CFLAGS := -O3 -I/home/weir/local/include/

# on my mac
# CFLAGS := -O0 -I/Users/weir/Installed/silo-4.8-bsd/include/

# on Alcyone:
# CFLAGS := -O3 -ipo -m64 -align -fno-alias -fno-fnalias \
#	-falign-functions -unroll-aggressive -fp-model fast=2 \
#	-I/home/weir/local/include/

# on Cosmos:
# CFLAGS := -g -O3 -xHost -align -ansi-alias -mcmodel=medium

#### LIBS ####

# on louhi
# LIBS := -lfftw3_mpi -lfftw3 -lsiloh5

# on sisu
# LIBS := -lfftw3_mpi -lfftw3 -lsiloh5 -lhdf5 -lstdc++ -lz

# on Taito
# LIBS := -lfftw3_mpi -lfftw3 -lsiloh5 -lhdf5 -lstdc++ -lz

# on vuori
# LIBS := -lfftw3_mpi -lfftw3 -lm

# on pc168:
# LIBS := -lfftw3_mpi -lfftw3 -lm -lsiloh5

# on afo
LIBS := -L/home/weir/local/lib/ \
	-lsiloh5 -lhdf5 -lfftw3_mpi -lfftw3 -lm -lstdc++

# on my mac:
# LIBS := -L/Users/weir/Installed/silo-4.8-bsd/lib/ \
#	-lsilo -lfftw3_mpi -lfftw3 -lm

# on Alcyone
# LIBS := -L/home/weir/local/lib/ \
#	-lfftw3_mpi -lfftw3 -lstdc++ -lhdf5 -lsiloh5 -lz

# on Cosmos:
# LIBS := -lmpi -lfftw3_mpi -lfftw3 -lsiloh5


#### No user-serviceable parts below ####

OBJECTS := main.o evolve.o potential.o energy.o eos.o \
	transport.o initial.o output.o parameters.o silage.o \
	util.o papi.o arrangement.o fft.o alloc.o gw.o velps.o \
	checkpoint.o uetc.o initps.o

BINARY := run/hydro # run/hydro-divps

FILES := $(OBJECTS) $(BINARY)

all: hydro

hydro: $(OBJECTS)
	$(CC) $(CFLAGS) $(OBJECTS) $(LIBS) -o $(BINARY)

clean:
	rm -f $(FILES) *\~ *\#
