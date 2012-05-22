/* transport.c
 *
 * Do advection. Two versions, depending on whether the
 * field being advected lives on lattice sites or between them.
 */
#include "hydro.h"



/* 1D port */
// advect zone quantities (ie state variable, eg E) 
void donor_E(hydro_fields f, int **nb, hydro_params p) {

  int x;

  // Flux field
  double *F = (double *) malloc(p.N*sizeof(double));
  // (Slow: see comments about this in eos.c)

  // Calculate flux
  // see advection chapter (4) PDF included
  // 1.0 in what follows will be upgraded to area of a cube later

  for(x = 0; x < p.N; x++) {
    if(f.Vx[x] >= 0.0)
      F[x] = f.Vx[x]*f.E[nb[x][1]];
    else
      F[x] = f.Vx[x]*f.E[x];
  }


  
  for(x = 0; x < p.N; x++)
    f.E[x] = f.E[x] + p.dt*(F[x] - F[nb[x][0]])/p.dx;
  

  free(F);

}


double phi(double x) {
  //return maxof3(0.0,minof2(1.0,2.0*x),minof2(2.0,x));
  //  return (x+fabs(x))/(1+fabs(x));
  return 1.0;
  //  return 0.0;
  //  return 0.5*(1.0 + x);
  //  return x;
  //  return 0.0;
}


void transport_E_dir(hydro_fields f, int **nb, hydro_params p, int dir) {


  int x;

  //  double *delta = (double *)malloc(p.N*sizeof(double));
  double *r = (double *)malloc(p.N*sizeof(double));
  double *F = (double *)malloc(p.N*sizeof(double));
  
  double Emin, Emax, dEmin, dEmax, dE;
  double dEsgn;
  double a, b;
  double theta;
 
  for(x = 0; x < p.N; x++) {
    /*
    Emin = minof3(f.E[x], f.E[nb[x][2*dir + 1]], f.E[nb[x][2*dir]]);
    Emax = maxof3(f.E[x], f.E[nb[x][2*dir + 1]], f.E[nb[x][2*dir]]);
    
    dE = (0.5*(f.E[x] + f.E[nb[x][0]]) - 0.5*(f.E[nb[x][2*dir + 1]] + f.E[x]));
    dEmin = 2.0*minof2(Emax - f.E[x],f.E[x] - Emin);
    dEmax = maxof3(Emax - f.E[x], f.E[x] - Emin, fabs(dE));
    
    if(dE < 0)
      dEsgn = -1.0;
    else
      dEsgn = 1.0;
    delta[x] = dEsgn*(minof2(2.0*dEmin, dEmax))/(p.dx);
    */
    
    if(fabs(f.E[nb[x][2*dir+1]] - f.E[nb[nb[x][2*dir+1]][2*dir+1]]) < 0.000001) {
      r[x] = 0.0;
    } else if(fabs(f.E[x] - f.E[nb[x][2*dir+1]]) < 0.000001) {
      r[x] = 0.0;
    } else if(f.V[dir][x] > 0) {
      r[x] = (f.E[nb[x][2*dir+1]] - f.E[nb[nb[x][2*dir+1]][2*dir+1]])
	/(f.E[x] - f.E[nb[x][2*dir+1]]);
    } else {
      r[x] = (f.E[nb[x][2*dir]] - f.E[x])
	/(f.E[x] - f.E[nb[x][2*dir+1]]);
    }
  }

  // Eq 2.116
  for(x = 0; x < p.N; x++) {
    if(f.V[dir][x] > 0)
      theta = 1.0;
    else
      theta = -1.0;

    /*
    if(f.V[dir][x] > 0)
      f.deltaM[dir][x] = (f.E[nb[x][2*dir + 1]] + 0.5*delta[nb[x][2*dir + 1]]
			  *(p.dx - 0.0*f.V[dir][x]*p.dt))*f.V[dir][x]*p.dt;
    theta = 1.0;
    else
      f.deltaM[dir][x] = (f.E[x] - 0.5*delta[x]
			  *(p.dx + 0.0*f.V[dir][x]*p.dt))*f.V[dir][x]*p.dt;
    */

    f.deltaM[dir][x] = 0.5*f.Vx[x]*((1.0+theta)*f.E[nb[x][2*dir+1]] + (1.0-theta)*f.E[x])
      + 0.5*fabs(f.Vx[x])*(1.0-fabs(f.Vx[x]*p.dt/p.dx))*phi(r[x])*(f.E[x] - f.E[nb[x][2*dir+1]]);

    if(isnan(0.5*fabs(f.Vx[x])*(1.0-fabs(f.Vx[x]*p.dt/p.dx))*phi(r[x])*(f.E[x] - f.E[nb[x][2*dir+1]]))) {

      fprintf(stderr,"went nan! x=%d\n",x);
      fprintf(stderr,"%lf %lf %lf<-%lf %lf\n", 0.5*fabs(f.Vx[x]),(1.0-fabs(f.Vx[x]*p.dt/p.dx)),phi(r[x]),r[x],(f.E[x] - f.E[nb[x][2*dir+1]]));
    }

  }

  // Eq 2.115
  for(x = 0; x < p.N; x++) {
    f.E[x] = f.E[x] + p.dt*(f.deltaM[dir][x] - f.deltaM[dir][nb[x][2*dir]])/p.dx;
    if(isnan(f.E[x])) {
      fprintf(stderr,"x=%d E went nan: deltaM = %lf, %lf\n", x, f.deltaM[dir][x], f.deltaM[dir][nb[x][2*dir]]);
    }
  }

  // free(delta);
  free(r);
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
    
  transport_E_dir(f, nb, p, order);
  transport_E_dir(f, nb, p, (order + 1) % 3);
  transport_E_dir(f, nb, p, (order + 2) % 3);    

}




void advect_Z(hydro_fields f, int **nb, hydro_params p) {

  int order; 
  order = lrand48() % 3;

  transport_Z_dir(f, nb, p, order);
  transport_Z_dir(f, nb, p, (order + 1) % 3);
  transport_Z_dir(f, nb, p, (order + 2) % 3);

}

