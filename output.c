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
double wallpos(double *xc, double *phi,
	       double *T, hydro_params params) {

  int x, xmax;

  double Vmax, Vtest;

  Vmax = Vf(params,T[0],phi[0]);
  xmax = 0;

  // Just search for maxmimum
  for(x=1; x<params.N; x++) {
    Vtest = Vf(params,T[x],phi[x]);
    if(Vtest > Vmax) {
      xmax = x;
      Vmax = Vtest;
    }
  }

  return xc[xmax];

}
