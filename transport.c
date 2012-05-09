/* transport.c
 *
 * Do advection. Two versions, depending on whether the
 * field being advected lives on lattice sites or between them.
 */
#include "hydro.h"


void donor_E_dir(hydro_fields f, int **nb, hydro_params p, int dir) {

  int x;

  // Flux field
  double *F = (double *) malloc(p.N*sizeof(double));

  for(x = 0; x < p.N; x++) {
      if(f.V[dir][x] >= 0.0)
	F[x] = f.V[dir][x]*f.E[x];
      else
	F[x] = f.V[dir][x]*f.E[nb[x][2*dir]];
  }

  
  for(x = 0; x < p.N; x++)
    f.E[x] = f.E[x] - p.dt*(F[x] - F[nb[x][2*dir + 1]]);


  free(F);
}


void donor_E(hydro_fields f, int **nb, hydro_params p) {

  int order = lrand48();
  order = order % 3;

  donor_E_dir(f, nb, p, order);
  donor_E_dir(f, nb, p, (order + 1) % 3);
  donor_E_dir(f, nb, p, (order + 2) % 3);

}




void donor_Z_dir(hydro_fields f, int **nb, hydro_params p, int dir) {


  double *F = (double *)malloc(p.N*sizeof(double));

  double vc;
  int x;

  // x-direction
  for(x = 0; x < p.N; x++) {

    vc = 0.5*(f.V[dir][nb[x][2*dir + 1]] + f.V[dir][x]);

    if(vc >= 0.0)
      F[x] = vc*f.Zx[nb[x][2*dir+1]];
    else
      F[x] = vc*f.Zx[x];
  }
  
  for(x = 0; x < p.N; x++) {
    f.Zx[x] = f.Zx[x] - p.dt*(F[nb[x][dir]] - F[x])/(p.dx);

  }  


  // y-direction
  for(x = 0; x < p.N; x++) {

    vc = 0.5*(f.V[dir][nb[x][2*dir + 1]] + f.V[dir][x]);

    if(vc >= 0.0)
      F[x] = vc*f.Zy[nb[x][2*dir+1]];
    else
      F[x] = vc*f.Zy[x];
  }
  
  for(x = 0; x < p.N; x++) {
    f.Zy[x] = f.Zy[x] - p.dt*(F[nb[x][dir]] - F[x])/(p.dx);

  }  


  // z-direction
  for(x = 0; x < p.N; x++) {

    vc = 0.5*(f.V[dir][nb[x][2*dir + 1]] + f.V[dir][x]);

    if(vc >= 0.0)
      F[x] = vc*f.Zz[nb[x][2*dir+1]];
    else
      F[x] = vc*f.Zz[x];
  }
  
  for(x = 0; x < p.N; x++) {
    f.Zz[x] = f.Zz[x] - p.dt*(F[nb[x][dir]] - F[x])/(p.dx);

  }  


  free(F);

}





void donor_Z(hydro_fields f, int **nb, hydro_params p) {
  donor_Z_dir(f, nb, p, 0);
  donor_Z_dir(f, nb, p, 1);
  donor_Z_dir(f, nb, p, 2);
}



/*
 * Fancier transport - van Leer.
 */
/*
void transport_E(hydro_fields f, int **nb, hydro_params p) {

  int x;

  double *F = (double *)malloc(p.N*sizeof(double));
  double *delta = (double *)malloc(p.N*sizeof(double));

  double s = p.dt/p.dx;
  
  double r;

  for(x = 0; x < p.N; x++) {
    r = (f.E[x] - f.E[nb[x][1]])*(f.E[nb[x][0]] - f.E[x]);

    if(r > 0)
      delta[x] = 2*r/(f.E[nb[x][0]] - f.E[nb[x][1]]);
    else
      delta[x] = 0.0;

  }

  delta[0] = 0.0;

  for(x = 0; x < p.N; x++) {
    if(f.v[x] > 0)
      F[x] = f.v[x]*1.0*(f.E[x] 
				     + 0.5*(1.0-f.v[x]*s)*delta[x]);
    else
      F[x] = f.v[x]*1.0*(f.E[nb[x][0]] 
				     - 0.5*(1.0+f.v[x]*s)*delta[nb[x][0]]);
  }

  F[0] = 0.0;  

  // "Advect D"
  for(x = 0; x < p.N; x++)
    f.E[x] = f.E[x] - s*(F[x] - F[nb[x][1]])/(1.0);

  free(delta);
  free(F);
}








void transport_Z(hydro_fields f, int **nb, hydro_params p) {

  int x;

  double *F = (double *)malloc(p.N*sizeof(double));
  double *delta = (double *)malloc(p.N*sizeof(double));

  double s = p.dt/p.dx;
  
  double r;
  double vc;

  for(x=0; x<p.N; x++) {
    r = (f.Z[x] - f.Z[nb[x][1]])*(f.Z[nb[x][0]] - f.Z[x]);

    if(r>0)
      delta[x] = 2*r/(f.Z[nb[x][0]] - f.Z[nb[x][1]]);
    else
      delta[x] = 0.0;

  }

  // Probably not
  //  delta[0] = 0.0;

  for(x = 0; x < p.N; x++) {
    vc = 0.5*(f.v[nb[x][1]] + f.v[x]);

    if(vc > 0)
      F[x] = vc*1.0*(f.Z[nb[x][1]] 
				 + 0.5*(1.0-vc*s)*delta[nb[x][1]]);
    else
      F[x] = vc*1.0*(f.Z[nb[x][0]]
				 - 0.5*(1.0+vc*s)*delta[x]);
  }

  F[0] = 0.0;  

  // "Advect Z"
  for(x = 0; x < p.N; x++)
    f.Z[x] = f.Z[x] - s*(F[nb[x][0]] - F[x])/(1.0);




  free(delta);
  free(F);
}
*/
