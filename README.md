# 3D hydro code

## The model

[Include links to papers here]

The potential is
\f[
V(\phi, T) = \frac{1}{2} \gamma(T^2 - T_0^2)\phi^2 - \frac{1}{3}\alpha T \phi^3 + \frac{1}{4}\lambda\phi^4
\f]

where, in the code:

* \f$ \gamma \f$ is `p.gamma`
* \f$ T_0 \f$ is `p.T0`
* \f$ \alpha \f$ is `p.alpha`
* \f$ \lambda \f$ is `p.lambda`

The internal energy of the fluid (`f.E`) is given by

\f[
\epsilon = 3 a T^4 + V(\phi, T) - T \frac{\partial V}{\partial T}
\f]

where:

* `f.E` is related to \f$ \epsilon \f$ by a gamma factor (`f.W`),
  thus: \f$ E = W \epsilon \f$
* \f$ V(\phi, T) \f$ is Vf()
* \f$ \frac{\partial V}{\partial T} \f$ is VTf()

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

* `-DSCALAR`: compile without the fluid.
