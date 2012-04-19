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
double wallpos(hydro_fields f, hydro_params p) {

  int x, xmax;

  double Vmax, Vtest;

  Vmax = Vf(p, f.T[0], f.phi[0]);
  xmax = 0;

  // Just search for maxmimum
  for(x = 1; x < p.N; x++) {
    Vtest = Vf(p, f.T[x], f.phi[x]);
    if(Vtest > Vmax) {
      xmax = x;
      Vmax = Vtest;
    }
  }

  return f.xc[xmax];

}




double get_gamma_max(hydro_fields f, hydro_params p) {

  int x, xmax;

  double gmax = f.gb[0];
  
  double gtest;

  // Just search for maxmimum
  for(x = 1; x < p.N; x++) {
    gtest = f.gb[x];
    if(gtest > gmax) {
      xmax = x;
      gmax = gtest;

    }
  }

  return gmax;

}
