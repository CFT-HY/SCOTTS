/* transport.c
 *
 * Do advection. Two versions, depending on whether the
 * field being advected lives on lattice sites or between them.
 */
#include "hydro.h"

// straight from fortran
void donor_r(double dt, double dx, int N, double *v,
	     double *xe, double *xc, double *field, int **nb) {

  double s = dt/dx;

  int x;

  // Slow: see comments about this in eos.c
  double *F = (double *)malloc(N*sizeof(double));

  for(x=0; x<N; x++) {
    if(v[x] >= 0.0)
      F[x] = v[x]*xe[x]*xe[x]*field[x];
    else
      F[x] = v[x]*xe[x]*xe[x]*field[nb[x][0]];
  }

  // BC's wrap
  F[0] = 0.0;

  for(x=0; x<N; x++)
    field[x] = field[x] - s*(F[x] - F[nb[x][1]])/(xc[x]*xc[x]);

  free(F);

}




// straight from fortran
void donor_Z_r(double dt, double dx, int N, double *v,
	       double *xe, double *xc, double *field, int **nb) {

  double s = dt/dx;

  double vc;
  int x;

  double *F = (double *)malloc(N*sizeof(double));

  for(x=0; x<N; x++) {
    vc = 0.5*(v[nb[x][1]] + v[x]);

    if(vc >= 0.0)
      F[x] = vc*xc[x]*xc[x]*field[nb[x][1]];
    else
      F[x] = vc*xc[x]*xc[x]*field[x];
  }

  for(x=0; x<N; x++)
    field[x] = field[x] - s*(F[nb[x][0]] - F[x])/(xe[x]*xe[x]);

  // BC's wrap
  field[0] = 0.0;

  free(F);

}
