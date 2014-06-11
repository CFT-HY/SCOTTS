#### CC ####

# afo
# CC := mpicc -g -DMPI -DFFT

# Sisu
# CC := cc -DMPI -DFFT

# Taito
CC := mpicc -DMPI -DFFT

#### CFLAGS ####

# on afo
# CFLAGS := -O3 -I/home/weir/local/include/

# on Sisu
# CFLAGS := -L/homeappl/home/weir/local_sisu/lib \
#	-I/homeappl/home/weir/local_sisu/include \
#	 -O3 -opt-prefetch -unroll-aggressive -no-prec-div \
#	-fp-model fast=2 -align -fno-alias -fno-fnalias -ipo

# on Taito 
CFLAGS := -L/homeappl/home/weir/local_taito/lib \
        -I/homeappl/home/weir/local_taito/include \
        -O3 -opt-prefetch -unroll-aggressive -no-prec-div \
        -fp-model fast=2 -align -fno-alias -fno-fnalias -ipo


#### LIBS ####

# on afo
# LIBS := -L/home/weir/local/lib/ \
#	-lsiloh5 -lhdf5 -lfftw3_mpi -lfftw3 -lm -lstdc++

# on sisu
# LIBS := -lfftw3_mpi -lfftw3 -lsiloh5 -lhdf5 -lstdc++ -lz

# on Taito
LIBS := -lfftw3_mpi -lfftw3 -lsiloh5 -lhdf5 -lstdc++ -lz


#### No user-serviceable parts below ####

OBJECTS := main.o evolve.o potential.o energy.o eos.o \
	transport.o initial.o output.o parameters.o silage.o \
	util.o papi.o arrangement.o fft.o alloc.o gw.o velps.o \
	checkpoint.o uetc.o initps.o

BINARY := run/hydro

FILES := $(OBJECTS) $(BINARY)

all: hydro

hydro: $(OBJECTS)
	$(CC) $(CFLAGS) $(OBJECTS) $(LIBS) -o $(BINARY)

clean:
	rm -f $(FILES) *\~ *\#

veryclean: clean
	rm -f *\~ *\#
