/* output.c
 * 
 * Things that we could in principle calculate offline, but
 * that are nice to know.
 */
#include "hydro.h"

/* double wallpos(...)
 *
 * Estimated position of wall.
 * Straight from fortran _BUT_ no interpolation of position.
 */
double wallpos(double *xc, int N, double *phi, double *Q,
	       double dt, double alpha, double gamma, double lambda,
	       double *T, double T0) {

  int x, xmax;

  double Vmax, Vtest;

  Vmax = Vf(alpha,gamma,lambda,T[0],T0,phi[0]);
  xmax = 0;

  // Just search for maxmimum
  for(x=1; x<N; x++) {
    Vtest = Vf(alpha,gamma,lambda,T[x],T0,phi[x]);
    if(Vtest > Vmax) {
      xmax = x;
      Vmax = Vtest;
    }
  }

  return xc[xmax];

}
