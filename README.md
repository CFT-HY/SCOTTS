# 3D hydro code

## Compilation

There are various makefiles in the `makefiles/` directory. By default,
most of them compile a parallelised (MPI) version of the code for the
platform identified in the makefile. They also all require FFTW to be
available, but where appropriate they can make use of system-provided
FFTW.

### Compilation flags:

* `-DSILO`: compile with support for writing 3D simulation data
  with Silo. This can include checkpointing (see `checkpoint.c`), or
  quantities relevant for visualisation (see `silage.c`). There is also
  routines which write only 2D slices in `silage_slice.c`. Once compiled
  with `-DSILO` the particular choice of target can be made with options
  in the parameter file.

* `-DFFT`: compile with FFT support for power spectra (and initial
  conditions).

* `-DINITPS`: compile with an initial power spectrum (requires
  `-DFFT`).

* `-DMPI`: compile with MPI support (necessary for parallelisation).