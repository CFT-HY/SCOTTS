/* potential.c
 *
 * The potential, and various derivatives thereof.
 */
#include "hydro.h"


// Potential
double Vf(hydro_params params, double T, double this_phi) {
  return  0.5*params.gamma*(T*T - params.T0*params.T0)*this_phi*this_phi
    - params.alpha*T*this_phi*this_phi*this_phi/3.0
    + 0.25*params.lambda*this_phi*this_phi*this_phi*this_phi;
}


// First derivative wrt phi
double Vdf(hydro_params params, double T, double this_phi) {
  return params.gamma*(T*T - params.T0*params.T0)*this_phi
    - params.alpha*T*this_phi*this_phi
    + params.lambda*this_phi*this_phi*this_phi;
}


// First derivative wrt T
double VTf(hydro_params params, double T, double this_phi) {
  return params.gamma*T*this_phi*this_phi 
    - params.alpha*this_phi*this_phi*this_phi/3.0;
}


// Second derivative wrt T
double VTTf(hydro_params params, double T, double this_phi) {
  return params.gamma*this_phi*this_phi;
}


/*  Vpot(...)
 *
 * Calculate potential for all sites. Except that doing it this way
 * does not encourage the compiler to fuse the loops.
 */
void Vpot(hydro_params params,
	  double *T,
	  double *phi, double *Vprecalc) {
  int x;

  for(x=0; x<params.N; x++)  {
    Vprecalc[x] = Vf(params, T[x], phi[x]);
  }
}


/* Vdpot(...)
 *
 * Calculate first deriviative of potential wrt phi for all sites.
 */
void Vdpot(hydro_params params,
	    double *T,
	    double *phi, double *Vprecalc) {
  int x;

  for(x=0; x<params.N; x++)  {
    Vprecalc[x] = Vdf(params, T[x], phi[x]);
  }
}
