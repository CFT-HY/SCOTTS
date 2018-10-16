# 3D hydro code

## The model 

[Include links to papers here]

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

* `-DDIMENSIONLESS`: Switch between dimensionless and dimensionful
damping couplings.

* `-DBAG` : Use bag model for equation of state/potential instead of
EIKR formalism.
