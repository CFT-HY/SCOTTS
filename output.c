/* output.c
 * 
 * Things that we could in principle calculate offline, but
 * that are nice to know.
 */
#include "hydro.h"



double get_gamma_max(hydro_fields f, hydro_params p) {

  int x, y, z, xmax;

  double gmax = f.W[0][0][0];
  
  double gtest;

  // Just search for maxmimum
  for(x = 1; x <= p.slicex; x++) {
    for(y = 1; y <= p.slicey; y++) {
      for(z = 0; z < p.Lz; z++) {
	gtest = f.W[x][y][z];
	if(gtest > gmax) {
	  gmax = gtest;
	}

      }
    }
  }
  
  return gmax;

}



void dump(double *field, hydro_params p) {
  int x;

  fprintf(stderr,"%g", field[0]);
  for(x=1;x<p.N;x++) {
    fprintf(stderr,", %g", field[x]);
  }
  fprintf(stderr,"\n");
}

