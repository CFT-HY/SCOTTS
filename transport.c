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
    f.E[x] = f.E[x] - p.dt*(F[x] - F[nb[x][2*dir + 1]])/(p.dx);


  free(F);
}

double minof3(double a, double b, double c) {
  if (a < b) {
    if (b < c) {
      return a;
    } else {
      if (c < a) {
	return c;
      } else {
	return a;
      }
    }
  } else if (b < c) {
    return b;
  }

  return c;  
}

double maxof3(double a, double b, double c) {
  if (a > b) {
    if (b > c) {
      return a;
    } else {
      if (c > a) {
	return c;
      } else {
	return a;
      }
    }
  } else if (b > c) {
    return b;
  }

  return c;  
}

double minof2(double a, double b) {
  if (a < b)
    return a;
  return b;
}

void transport_E_dir(hydro_fields f, int **nb, hydro_params p, int dir) {


  int x;

  double *delta = (double *)malloc(p.N*sizeof(double));

  double Emin, Emax, dEmin, dEmax, dE;
  double dEsgn;
 
  for(x = 0; x < p.N; x++) {


    Emin = minof3(f.E[x], f.E[nb[x][2*dir + 1]], f.E[nb[x][2*dir]]);
    Emax = maxof3(f.E[x], f.E[nb[x][2*dir + 1]], f.E[nb[x][2*dir]]);

    dE = f.E[x] - f.E[nb[x][2*dir + 1]];
    dEmin = 2.0*minof2(Emax - f.E[x],f.E[x] - Emin);
    dEmax = maxof3(Emax - f.E[x], f.E[x] - Emin, fabs(dE));

    if(dE < 0)
      dEsgn = -1.0;
    else
      dEsgn = 1.0;

    delta[x] = dEsgn*(minof2(2*dEmin, dEmax))/(p.dx);


  }


  // Eq 2.116
  for(x = 0; x < p.N; x++) {
    if(f.V[dir][x] > 0)
      f.deltaM[dir][x] = (f.E[nb[x][2*dir + 1]] + 0.5*delta[nb[x][2*dir + 1]]*(p.dx - f.V[dir][x]*p.dt))*f.V[dir][x];
    else
      f.deltaM[dir][x] = (f.E[x] - 0.5*delta[x]*(p.dx + f.V[dir][x]*p.dt))*f.V[dir][x];

  }

  // Eq 2.115
  for(x = 0; x < p.N; x++) {
    f.E[x] = f.E[x] - p.dt*(f.deltaM[dir][nb[x][2*dir]] - f.deltaM[dir][x])/p.dx;
  }



  free(delta);
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



void transport_Z_dir(hydro_fields f, int **nb, hydro_params p, int dir) {


  int x;

  double *F = (double *)malloc(p.N*sizeof(double));
  double *delta = (double *)malloc(p.N*sizeof(double));

  double *Vbody = (double *)malloc(p.N*sizeof(double));
  double *Ubody = (double *)malloc(p.N*sizeof(double));

  double **deltaMI = (double **)malloc(3*sizeof(double *));
  deltaMI[0] = (double *)malloc(p.N*sizeof(double));
  deltaMI[1] = (double *)malloc(p.N*sizeof(double));
  deltaMI[2] = (double *)malloc(p.N*sizeof(double));

  double **deltaMIb = (double **)malloc(3*sizeof(double *));
  deltaMIb[0] = (double *)malloc(p.N*sizeof(double));
  deltaMIb[1] = (double *)malloc(p.N*sizeof(double));
  deltaMIb[2] = (double *)malloc(p.N*sizeof(double));



  double Umin, Umax, dUmin, dUmax, dU;
  double dUsgn;

  // Eq 2.119 -- flux
  for(x=0; x<p.N; x++) {
    deltaMI[0][x] = f.deltaM[0][x]*(f.kappa[x] + f.kappa[nb[x][1]])*0.5;
    deltaMI[1][x] = f.deltaM[1][x]*(f.kappa[x] + f.kappa[nb[x][3]])*0.5;
    deltaMI[2][x] = f.deltaM[2][x]*(f.kappa[x] + f.kappa[nb[x][5]])*0.5;
  }


  // Average flux over node-centred cube (these are the f_i's)
  for(x=0; x<p.N; x++) {
    deltaMIb[0][x] = 0.25*(deltaMI[0][x] + deltaMI[0][nb[x][3]] + deltaMI[0][nb[x][5]] + deltaMI[0][nb[nb[x][3]][5]])*p.dt;
    deltaMIb[1][x] = 0.25*(deltaMI[1][x] + deltaMI[1][nb[x][1]] + deltaMI[1][nb[x][5]] + deltaMI[1][nb[nb[x][1]][5]])*p.dt;
    deltaMIb[2][x] = 0.25*(deltaMI[2][x] + deltaMI[2][nb[x][1]] + deltaMI[2][nb[x][3]] + deltaMI[2][nb[nb[x][1]][3]])*p.dt;
  }



  // Eq 2.120 -- already obtained? also 2.121, 2.122
  for(x = 0; x < p.N; x++) {
    Vbody[x] = 0.125*(f.V[dir][x] + f.V[dir][nb[x][0]] + f.V[dir][nb[x][2]] + f.V[dir][nb[x][4]] + f.V[dir][nb[nb[x][0]][2]] + f.V[dir][nb[nb[x][0]][4]] 
		     + f.V[dir][nb[nb[x][2]][4]] + f.V[dir][nb[nb[nb[x][0]][2]][4]])*p.dt;


    Ubody[x] = 0.125*(f.U[dir][x] + f.U[dir][nb[x][0]] + f.U[dir][nb[x][2]] + f.U[dir][nb[x][4]] + f.U[dir][nb[nb[x][0]][2]] + f.U[dir][nb[nb[x][0]][4]] 
		     + f.U[dir][nb[nb[x][2]][4]] + f.U[dir][nb[nb[nb[x][0]][2]][4]])*p.dt;
    

  }

  // 2.123-2.127
  for(x = 0; x < p.N; x++) {


    Umin = minof3(f.U[dir][x], f.U[dir][nb[x][2*dir + 1]], f.U[dir][nb[x][2*dir]]);
    Umax = maxof3(f.U[dir][x], f.U[dir][nb[x][2*dir + 1]], f.U[dir][nb[x][2*dir]]);


    dUmin = 2.0*minof2(Umax - f.U[dir][x],f.U[dir][x] - Umin);
    dU = Ubody[x] - Ubody[nb[x][2*dir + 1]];
    
    if(dU < 0)
      dUsgn = -1.0;
    else
      dUsgn = 1.0;

    delta[x] = dUsgn*(minof2(dUmin, fabs(dU)))/(p.dx);


  }


  // Eq 2.116
  for(x = 0; x < p.N; x++) {
    if(Vbody[x] > 0)
      F[x] = (f.U[dir][x] + 0.5*delta[x]*(p.dx - Vbody[x]*p.dt));
    else
      F[x] = (f.U[dir][x] - 0.5*delta[nb[x][2*dir]]*(p.dx - Vbody[x]*p.dt));

  }


  // Eq 2.115
  for(x = 0; x < p.N; x++) {
    f.E[x] = f.E[x] - p.dt*(deltaMIb[0][x]*F[x] - deltaMIb[0][nb[x][0]]*F[nb[x][0]])/p.dx;
    f.E[x] = f.E[x] - p.dt*(deltaMIb[1][x]*F[x] - deltaMIb[1][nb[x][2]]*F[nb[x][2]])/p.dx;
    f.E[x] = f.E[x] - p.dt*(deltaMIb[2][x]*F[x] - deltaMIb[2][nb[x][4]]*F[nb[x][4]])/p.dx;
  }



  free(delta);
  free(F);
  free(Vbody);
  free(Ubody);
  free(deltaMI[0]);
  free(deltaMI[1]);
  free(deltaMI[2]);
  free(deltaMI);
  free(deltaMIb[0]);
  free(deltaMIb[1]);
  free(deltaMIb[2]);
  free(deltaMIb);


}



void advect_E(hydro_fields f, int **nb, hydro_params p) {

  int order; 
  order = lrand48() % 3;
  //  order = 0;

    
  transport_E_dir(f, nb, p, order);
  transport_E_dir(f, nb, p, (order + 1) % 3);
  transport_E_dir(f, nb, p, (order + 2) % 3);
    

  /*
  donor_E_dir(f, nb, p, order);
  donor_E_dir(f, nb, p, (order + 1) % 3);
  donor_E_dir(f, nb, p, (order + 2) % 3);
  */

}




void advect_Z(hydro_fields f, int **nb, hydro_params p) {
  /* 
  donor_Z_dir(f, nb, p, 0);
  donor_Z_dir(f, nb, p, 1);
  donor_Z_dir(f, nb, p, 2);
  */

  int order; 
  order = lrand48() % 3;

  transport_Z_dir(f, nb, p, order);
  transport_Z_dir(f, nb, p, (order + 1) % 3);
  transport_Z_dir(f, nb, p, (order + 2) % 3);

}



/*
 * Fancier transport - van Leer.
 */
/*







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
