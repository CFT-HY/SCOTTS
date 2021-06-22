# 3D hydro code

Papers using/describing this code:
https://arxiv.org/abs/1511.04527
https://arxiv.org/abs/1704.05871
https://arxiv.org/abs/1906.00480

## The model

Under normal compliation the potential is
\f[ V_B(\phi, T) = \frac{1}{2} \gamma(T^2 - T_0^2)\phi^2
- \frac{1}{3}\alpha T \phi^3 +\frac{1}{4}\lambda\phi^4\text.  \f]

where, in the code:

* \f$ \gamma \f$ is `p.gamma`
* \f$ T_0 \f$ is `p.T0`
* \f$ \alpha \f$ is `p.alpha`
* \f$ \lambda \f$ is `p.lambda`

We renormalise so the potential energy is zero at zero temperature and
in the broken phase minimum \f$ \phi_b(T) \f$:
\f[ V(\phi, T) = V_B(\phi,T) -  V_B(\phi_b(0),0) \text,  \f]
and in the code \f$ V_B(\phi_b(0),0)  \f$ is `p.V0`.

The pressure is given by
\f[
p =  a T^4 - V(\phi, T)
\f]
and the internal energy is found through
\f[
\epsilon = T\frac{\partial p}{\partial T} - p
\f]
\f[
\epsilon = 3 a T^4 + V(\phi, T) - T \frac{\partial V}{\partial T}
\f]
with \f$ a = (\pi^2/90)g_* \f$ and \f$ g_* \f$ the effective degrees
of freedom.

In the code:
* `f.p` is the pressure
* The fluid energy density `f.E` is related to \f$ \epsilon \f$ by a
  gamma factor (`f.W`),  thus: \f$ E = W \epsilon \f$
* \f$ V(\phi, T) \f$ is Vf()
* \f$ \dfrac{\partial V}{\partial T} \f$ is VTf()
* \f$ a \f$ is `p.gdeg` and \f$ g_* \f$ is `p.gstar`


Also present in the code is the variable `f.kappa` which is given by
\f[
\kappa=1+\dfrac{PW}{E}
\f]

We follow Relativistic Numerical Hydrodynamics by
Wilson and Mathews to evolve the fluid. However, we define the
internal energy \f$\epsilon\f$ differently:
* \f$ \epsilon= \rho_\text{WM} (1 + \epsilon_\text{WM}) \f$

where the subscript WM indicates a variable given in Wilson and Mathews.

This then leads to:
* \f$ E= D_\text{WM} + E_\text{WM} \f$
* \f$ \dfrac{\kappa-1}{\Gamma_\text{WM}-1}=\dfrac{E_\text{WM}}{E} \f$

# BAG model

Alternatively we can compile with the `BAG` compiler flag to use the
bag model for the equation of state and the potential.

We have the potential
\f[
 V_A(\phi)=\frac{1}{2}\gamma \phi^2 - \frac{1}{3}\alpha \phi^3 +
  \frac{1}{4} \lambda \phi^4\text,
\f]
where, in the code:

* \f$ \gamma \f$ is `p.gamma`
* \f$ \alpha \f$ is `p.alpha`
* \f$ \lambda \f$ is `p.lambda`

We normalise the potential to be zero in the broken phase minima
\f$\phi_b \f$,
\f[
V_B(\phi)= V_A(\phi) - V_A(\phi_b)\text.
\f]
and in the code \f$ V_A(\phi_b) \f$ is `p.V0`.

In the bag model the effective degrees of freedom vary with \f$\phi\f$
and so \f$ a \rightarrow a(\phi) \f$.
\f[
p =  a(\phi) T^4 - V_B(\phi)
\f]
\f[
\epsilon = 3a(\phi) T^4 + V_B(\phi)
\f]
with
\f[
  a(\phi)= a_0 - \Delta V\left[3 \left(\frac{\phi}{\phi_b}\right)^2 -
    2 \left(\frac{\phi}{\phi_b}\right)^3\right]
\f]
and \f$ \Delta V = V_B(0) - V_B(\phi)\f$ the potential energy
difference between the two states. For our potential \f$ \Delta V\f$ is
`-p.V0`. We have \f$a_0= (\pi^2/90)g_*\f$.

We can group all terms that depend on \f$\phi\f$ into a temperature
dependent potential
\f[
V(\phi,T) = V_B(\phi) - (a(\phi) - a_0) T^4\text,
\f]
with
\f[
p =  a_0 T^4 - V(\phi,T)
\f]
\f[
\epsilon = 3a_0 T^4 + V(\phi, T) - T \frac{\partial V}{\partial T}\text.
\f]

In the code:
* `f.p` is the pressure
* The fluid energy density `f.E` is related to \f$ \epsilon \f$ by a
  gamma factor (`f.W`),  thus: \f$ E = W \epsilon \f$
* \f$ V(\phi, T) \f$ is Vf()
* \f$ \dfrac{\partial V}{\partial T} \f$ is VTf()
* \f$ a_0 \f$ is `p.gdeg` and \f$ g_* \f$ is `p.gstar`

## Dimensionless vs dimensionful friction parameter

The fluid is coupled to the field through a dissipative term proportional
to field gradient.
\f[
\left[ \partial_\mu T^{\mu\nu}\right]_\text{field}=
-\left[ \partial_\mu T^{\mu\nu}\right]_\text{fluid}=
\eta U^\mu \partial_\mu \phi \partial^\nu \phi
\f]
with \f$ T^{\mu \nu} \f$ the energy momentum tensor and \f$ \eta \f$
the friction parameter.

The default choice for the friction parameter is that \f$ \eta \f$ is a
dimensionful constant set by `C` in the input parameter file.

This can be changed using the `-DDIMENSIONLESS` compiler flag. Then
\f[
\eta = \tilde{\eta} \frac{\phi^2}{T}
\f]
with \f$\tilde{\eta}\f$ a constant set by `C` in the input parameter file.


## Order of evolution:

This is non unique and can be modified. This was settled on as it
appears to have good energy conservation.

-# Scalar field is evolved [evolve_field()]
-# Temperature, pressure and \f$\kappa\f$ calculated [eq_of_state()]
-# Update \f$E\f$ and \f$Z \f$ with field-fluid interaction terms
     [evolve_hydro()]
-# Update \f$ Z \f$ with pressure acceleration [evolve_hydro()]
-# Update covariant 4-velocity spatial terms \f$ U_i \f$ [evolve_hydro()]
-# Update contravariant 3-velocity  \f$ V^i \f$ [evolve_hydro()]
-# Update \f$ E \f$ with PdV work terms [evolve_hydro()]
-# Evolve metric perturbations \f$ u_{ij} \f$ [evolve_uij()]
-# Advection of \f$E\f$ and \f$ Z \f$ [advect_E() & advect_Z()]
-# Find temperature again [find_Ta()]

## Initial conditions with a given power spectrum INITPS

To launch a simulation with a random gaussian velocity field following a given power spectrum
\f[
\frac{\mathrm{d} \langle U^2 \rangle}{\mathrm{d}\ln k}
\f].

Requirements :
- `-DFFT` flag in the Makefile, compile with support for FFTs
- `-DBAG` flag in the Makefile, compiles with the bag model equation of state.
- `initial initps` in the initialization file to enable initps
- `initpsfile initial.txt div` to specify the type of initial conditions. `initial.txt` is the location of the input power spectrum. `div`(resp. `rot`, `all`) stands for longitudinal (resp. vortical, non-projected) initial conditions
- `initpsbins 200` the number of points in the power spectrum

### Mock input spectrum

`init_file.py` is a Python3 file to mock input power spectrum.
It produces power spectra according to the formula :

\f[\frac{{\rm d} \langle v^2 \rangle}{{\rm d} \ln k} =
C \frac{k^p}{(k_{peak}^s + k^s)^{(q-p)/s}} \exp\left(-\frac{k}{k_{max}}\right)\f]

with parameters
- `RMS_VELOCITY`: float
    Root mean square velocity of the fluid
- `INTEGRAL_SCALE`: float
   Integral scale of the fluid, in physical units
- `IR_SLOPE`: float
    Slope of the power spectrum in the infrared, $p$ in the formula
- `UV_SLOPE`: float
    Slope of the power spectrum in the ultraviolet, $q$ in the formula
- `SKEW`: float
    Parameter $s$ in the formula
- `K_MAX`: float
    Wavenumber associated with the ultraviolet cutoff
- `LATTICE_SIZE`: int
    Number of lattice sites in one direction (assumes Lx == Ly == Lz)
- `LATTICE_SPACING`: float
    Size of a given lattice site
- `CUTOFF`: float
    Hard cutoff, just in case

It relies on the libraries `numpy` and `scipy`.

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


* `-DUSE_MPI`: compile with MPI support (necessary for parallelisation).

* `-DSCALAR`: compile without the fluid.

* `-DDIMENSIONLESS`: Switch between dimensionless and dimensionful
damping couplings.

* `-DBAG` : Use bag model for equation of state/potential instead of
EIKR formalism.

* `-DVANLEER` Use Van Leer advection instead of donor cell.
