/* potential.c
 *
 * The potential, and various derivatives thereof.
 */
#include "hydro.h"


// Potential
double Vf(double alpha, double gamma, double lambda,
	  double T, double T0, double this_phi) {
  return  0.5*gamma*(T*T - T0*T0)*this_phi*this_phi
    - alpha*T*this_phi*this_phi*this_phi/3.0
    + 0.25*lambda*this_phi*this_phi*this_phi*this_phi;
}


// First derivative wrt phi
double Vdf(double alpha, double gamma, double lambda,
	   double T, double T0, double this_phi) {
  return gamma*(T*T - T0*T0)*this_phi
    - alpha*T*this_phi*this_phi
    + lambda*this_phi*this_phi*this_phi;
}


// First derivative wrt T
double VTf(double alpha, double gamma, double lambda,
	   double T, double T0, double this_phi) {
  return gamma*T*this_phi*this_phi 
    - alpha*this_phi*this_phi*this_phi/3.0;
}


// Second derivative wrt T
double VTTf(double alpha, double gamma, double lambda,
	    double T, double T0, double this_phi) {
  return gamma*this_phi*this_phi;
}


/*  Vpot(...)
 *
 * Calculate potential for all sites. Except that doing it this way
 * does not encourage the compiler to fuse the loops.
 */
void Vpot(int N, double alpha, double gamma,
	  double lambda, double *T, double T0,
	  double *phi, double *Vprecalc) {
  int x;

  for(x=0; x<N; x++)  {
    Vprecalc[x] = Vf(alpha, gamma, lambda, T[x], T0, phi[x]);
  }
}


/* Vdpot(...)
 *
 * Calculate first deriviative of potential wrt phi for all sites.
 */
void Vdpot(int N, double alpha, double gamma,
	    double lambda, double *T, double T0,
	    double *phi, double *Vprecalc) {
  int x;

  for(x=0; x<N; x++)  {
    Vprecalc[x] = Vdf(alpha, gamma, lambda, T[x], T0, phi[x]);
  }
}
