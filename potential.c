/* potential.c
 *
 * The potential, and various derivatives thereof.
 */
#include "hydro.h"


// Potential
double Vf(hydro_params p, double T, double this_phi) {
  return  0.5*p.gamma*(T*T - p.T0*p.T0)*this_phi*this_phi
    - p.alpha*T*this_phi*this_phi*this_phi/3.0
    + 0.25*p.lambda*this_phi*this_phi*this_phi*this_phi + V0;
}


// First derivative wrt phi
double Vdf(hydro_params p, double T, double this_phi) {
  return p.gamma*(T*T - p.T0*p.T0)*this_phi
    - p.alpha*T*this_phi*this_phi
    + p.lambda*this_phi*this_phi*this_phi;
}


// First derivative wrt T
double VTf(hydro_params p, double T, double this_phi) {
  return p.gamma*T*this_phi*this_phi 
    - p.alpha*this_phi*this_phi*this_phi/3.0;
}


// Second derivative wrt T
double VTTf(hydro_params p, double T, double this_phi) {
  return p.gamma*this_phi*this_phi;
}


/*  Vpot(...)
 *
 * Calculate potential for all sites. Except that doing it this way
 * does not encourage the compiler to fuse the loops.
 */
void Vpot(hydro_params p,
	  double *T,
	  double *phi, double *Vprecalc) {
  int x;

  for(x=0; x<p.fieldN; x++)  {
    Vprecalc[x] = Vf(p, T[x], phi[x]);
  }
}


/* Vdpot(...)
 *
 * Calculate first deriviative of potential wrt phi for all sites.
 */
void Vdpot(hydro_params p,
	    double *T,
	    double *phi, double *Vprecalc) {
  int x;

  for(x=0; x<p.fieldN; x++)  {
    Vprecalc[x] = Vdf(p, T[x], phi[x]);
  }
}
