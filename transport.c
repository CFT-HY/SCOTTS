/* transport.c
 *
 * Do advection. Two versions, depending on whether the
 * field being advected lives on lattice sites or between them.
 */
#include "hydro.h"

// straight from fortran

// advect zone quantities (ie state variable, eg E)
void donor_r(double dt, double dx, int N, double *v,
	     double *xe, double *xc, double *field, int **nb) {

  double s = dt/dx;

  int x;

  // Flux field
  double *F = (double *)malloc(N*sizeof(double));
  // (Slow: see comments about this in eos.c)



  // Calculate flux
  // see advection chapter (4) PDF included
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
void donor_z(double dt, double dx, int N, double *v,
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




/*
 * Fancier transport - van Leer.
 */

void transport_r(double dt, double dx, int N, double *v,
	     double *xe, double *xc, double *field, int **nb) {

  int x;

  double *F = (double *)malloc(N*sizeof(double));
  double *delta = (double *)malloc(N*sizeof(double));

  double s = dt/dx;
  
  double r;

  for(x=0; x<N; x++) {
    r = (field[x] - field[nb[x][1]])*(field[nb[x][0]] - field[x]);

    if(r>0)
      delta[x] = 2*r/(field[nb[x][0]] - field[nb[x][1]]);
    else
      delta[x] = 0.0;

  }

  delta[0] = 0.0;

  for(x=0; x<N; x++) {
    if(v[x] > 0)
      F[x] = v[x]*xe[x]*xe[x]*(field[x] + 0.5*(1.0-v[x]*s)*delta[x]);
    else
      F[x] = v[x]*xe[x]*xe[x]*(field[nb[x][0]] - 0.5*(1.0+v[x]*s)*delta[nb[x][0]]);
  }

  F[0] = 0.0;  

  // "Advect D"
  for(x=0; x<N; x++)
    field[x] = field[x] - s*(F[x] - F[nb[x][1]])/(xc[x]*xc[x]);




  free(delta);
  free(F);
}








void transport_z(double dt, double dx, int N, double *v,
		 double *xe, double *xc, double *field, int **nb) {

  int x;

  double *F = (double *)malloc(N*sizeof(double));
  double *delta = (double *)malloc(N*sizeof(double));

  double s = dt/dx;
  
  double r;
  double vc;

  for(x=0; x<N; x++) {
    r = (field[x] - field[nb[x][1]])*(field[nb[x][0]] - field[x]);

    if(r>0)
      delta[x] = 2*r/(field[nb[x][0]] - field[nb[x][1]]);
    else
      delta[x] = 0.0;

  }

  // Probably not
  //  delta[0] = 0.0;

  for(x=0; x<N; x++) {
    vc = 0.5*(v[nb[x][1]] + v[x]);

    if(vc > 0)
      F[x] = vc*xc[x]*xc[x]*(field[nb[x][1]] + 0.5*(1.0-vc*s)*delta[nb[x][1]]);
    else
      F[x] = vc*xc[x]*xc[x]*(field[nb[x][0]] - 0.5*(1.0+vc*s)*delta[x]);
  }

  F[0] = 0.0;  

  // "Advect Z"
  for(x=0; x<N; x++)
    field[x] = field[x] - s*(F[nb[x][0]] - F[x])/(xe[x]*xe[x]);




  free(delta);
  free(F);
}
