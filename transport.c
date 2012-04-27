/* transport.c
 *
 * Do advection. Two versions, depending on whether the
 * field being advected lives on lattice sites or between them.
 */
#include "hydro.h"

// straight from fortran

// advect zone quantities (ie state variable, eg E)
/*
void donor_E(hydro_fields f, int **nb, hydro_params p) {

  double s = p.dt/p.dx;

  int x;

  // Flux field
  double *F = (double *) malloc(p.N*sizeof(double));
  // (Slow: see comments about this in eos.c)



  // Calculate flux
  // see advection chapter (4) PDF included
  // 1.0 in what follows will be upgraded to area of a cube later
  for(x = 0; x < p.N; x++) {
    if(f.v[x] >= 0.0)
      F[x] = f.v[x]*1.0*f.E[x];
    else
      F[x] = f.v[x]*1.0*f.E[nb[x][0]];
  }

  // BC's wrap
  F[0] = 0.0;


  for(x = 0; x < p.N; x++)
    f.E[x] = f.E[x] - s*(F[x] - F[nb[x][1]])/(1.0);

  free(F);

}
*/


void donor_E(hydro_fields f, int **nb, hydro_params p) {

  double s = p.dt/p.dx;

  int x;

  // Flux field
  double *F = (double *) malloc(p.N*sizeof(double));
  // (Slow: see comments about this in eos.c)

  

  // Calculate flux
  // see advection chapter (4) PDF included
  // 1.0 in what follows will be upgraded to area of a cube later

  for(x = 0; x < p.N; x++) {
      if(f.Vx[x] >= 0.0)
	F[x] = f.Vx[x]*1.0*f.E[x];
      else
	F[x] = f.Vx[x]*1.0*f.E[nb[x][0]];
  }

  for(x = 0; x < p.N; x++)
    f.E[x] = f.E[x] - s*(F[x] - F[nb[x][1]])/(p.dx);

  for(x=0; x<p.N; x++) {
    if(f.Vy[x] >= 0.0)
      F[x] = f.Vy[x]*1.0*f.E[x];
    else
      F[x] = f.Vy[x]*1.0*f.E[nb[x][2]];
  }

  for(x = 0; x < p.N; x++)
    f.E[x] = f.E[x] - s*(F[x] - F[nb[x][3]])/(p.dx);

  for(x=0; x<p.N; x++) {
    if(f.Vz[x] >= 0.0)
	F[x] = f.Vz[x]*1.0*f.E[x];
    else
      F[x] = f.Vz[x]*1.0*f.E[nb[x][4]];
  }

  for(x = 0; x < p.N; x++)
    f.E[x] = f.E[x] - s*(F[x] - F[nb[x][5]])/(p.dx);

  free(F);
}



/*
// straight from fortran
void donor_Z(hydro_fields f, int **nb, hydro_params p) {

  double s = p.dt/p.dx;

  double vc;
  int x;

  double *F = (double *)malloc(p.N*sizeof(double));

  for(x=0; x<p.N; x++) {
    vc = 0.5*(f.v[nb[x][1]] + f.v[x]);

    if(vc >= 0.0)
      F[x] = vc*1.0*f.Z[nb[x][1]];
    else
      F[x] = vc*1.0*f.Z[x];
  }

  for(x=0; x<p.N; x++)
    f.Z[x] = f.Z[x] - s*(F[nb[x][0]] - F[x])/(1.0);

  // BC's wrap
  f.Z[0] = 0.0;

  free(F);

}
*/

void donor_Z(hydro_fields f, int **nb, hydro_params p) {

  double s = p.dt/p.dx;


  int x, y, z;

  double *F = (double *)malloc(p.N*sizeof(double));

  double Vb;

  for(x = 0; x < p.N; x++) {
    if(f.Vx[x] >= 0.0)
      F[x] = f.Vx[x]*1.0*f.E[x];
    else
      F[x] = f.Vx[x]*1.0*f.E[nb[x][0]];
  }
  
  for(x = 0; x < p.N; x++)
    f.Zx[x] = f.Zx[x] - s*(F[x] - F[nb[x][1]])/(p.dx);
  
  for(x=0; x<p.N; x++) {
    if(f.Vy[x] >= 0.0)
      F[x] = f.Vy[x]*1.0*f.E[x];
    else
      F[x] = f.Vy[x]*1.0*f.E[nb[x][2]];
  }
	
  for(x = 0; x < p.N; x++)
    f.Zx[x] = f.Zx[x] - s*(F[x] - F[nb[x][3]])/(p.dx);
	
  for(x=0; x<p.N; x++) {
    if(f.Vz[x] >= 0.0)
      F[x] = f.Vz[x]*1.0*f.E[x];
    else
      F[x] = f.Vz[x]*1.0*f.E[nb[x][4]];
  }
  
  for(x = 0; x < p.N; x++)
    f.Zx[x] = f.Zx[x] - s*(F[x] - F[nb[x][5]])/(p.dx);





  // Y

  for(x = 0; x < p.N; x++) {
    if(f.Vx[x] >= 0.0)
      F[x] = f.Vx[x]*1.0*f.E[x];
    else
      F[x] = f.Vx[x]*1.0*f.E[nb[x][0]];
  }
  
  for(x = 0; x < p.N; x++)
    f.Zy[x] = f.Zy[x] - s*(F[x] - F[nb[x][1]])/(p.dx);
  
  for(x=0; x<p.N; x++) {
    if(f.Vy[x] >= 0.0)
      F[x] = f.Vy[x]*1.0*f.E[x];
    else
      F[x] = f.Vy[x]*1.0*f.E[nb[x][2]];
  }
	
  for(x = 0; x < p.N; x++)
    f.Zy[x] = f.Zy[x] - s*(F[x] - F[nb[x][3]])/(p.dx);
	
  for(x=0; x<p.N; x++) {
    if(f.Vz[x] >= 0.0)
      F[x] = f.Vz[x]*1.0*f.E[x];
    else
      F[x] = f.Vz[x]*1.0*f.E[nb[x][4]];
  }
  
  for(x = 0; x < p.N; x++)
    f.Zy[x] = f.Zy[x] - s*(F[x] - F[nb[x][5]])/(p.dx);


  // Z

  for(x = 0; x < p.N; x++) {
    if(f.Vx[x] >= 0.0)
      F[x] = f.Vx[x]*1.0*f.E[x];
    else
      F[x] = f.Vx[x]*1.0*f.E[nb[x][0]];
  }
  
  for(x = 0; x < p.N; x++)
    f.Zz[x] = f.Zz[x] - s*(F[x] - F[nb[x][1]])/(p.dx);
  
  for(x=0; x<p.N; x++) {
    if(f.Vy[x] >= 0.0)
      F[x] = f.Vy[x]*1.0*f.E[x];
    else
      F[x] = f.Vy[x]*1.0*f.E[nb[x][2]];
  }
	
  for(x = 0; x < p.N; x++)
    f.Zz[x] = f.Zz[x] - s*(F[x] - F[nb[x][3]])/(p.dx);
	
  for(x=0; x<p.N; x++) {
    if(f.Vz[x] >= 0.0)
      F[x] = f.Vz[x]*1.0*f.E[x];
    else
      F[x] = f.Vz[x]*1.0*f.E[nb[x][4]];
  }
  
  for(x = 0; x < p.N; x++)
    f.Zz[x] = f.Zz[x] - s*(F[x] - F[nb[x][5]])/(p.dx);
  


	/*
	Vb = 0.5*(f.Vx[iix(x,y,z)] + f.Vx[iix(x+1,y,z)]);



	if(Vb >= 0.0)
	  F[x] = Vb*1.0*f.Zx[iix(x-1,y,z)];
	else
	  F[x] = Vb*1.0*f.Zx[iix(x,y,z)];
	*/


  free(F);

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
