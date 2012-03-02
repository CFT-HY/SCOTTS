/* evolve.c
 *
 * Everything to do with evolving fields and fluids with the
 * equations of motion.
 */
#include "hydro.h"


// straight from fortran
void evolve_backstep(double dt, double dx, int N, double *phi,
		     double *pifull, double *xe, double *xc,
		     double alpha, double gamma, double lambda,
		     double *T, double T0, double *pi, int **nb) {

  int x;

  for(x=0; x<N; x++)
    pi[x] = pifull[x] - 0.5*dt/(xc[x]*xc[x]*dx*dx)
      * (xe[x]*xe[x]*(phi[nb[x][0]] - phi[x]) 
	 - xe[nb[x][1]]*xe[nb[x][1]]*(phi[x] - phi[nb[x][1]]))
      + 0.125*(dt/dx/xc[x])*(dt/dx/xc[x])
      * (xe[x]*xe[x]*(pifull[nb[x][0]] - pifull[x]) 
	 - xe[nb[x][1]]*xe[nb[x][1]]*(pifull[x] - pifull[nb[x][1]]))
      + 0.5*dt*Vdf(alpha, gamma, lambda, T[x], T0, phi[x]-0.25*dt*pifull[x]);

}


// straight from fortran, except set v = 0 explicitly rather than in wrap
void evolve_field(double dt, double dx, double C, int N,
		  double *gb, double *v, double *xe, double *xc,
		  double *pi, double *phi, double alpha, double gamma,
		  double lambda, double *T, double T0,
		  double *phiold, double *pifull, int **nb) {
  

  // Slow and not strictly necessary, see comments in eos.c
  double *Vdold = (double *)malloc(N*sizeof(double));

  int x;

  double piold, s;


  Vdpot(N, alpha, gamma, lambda, T, T0, phi, Vdold);

  // Move conjugate momentum
  for(x=0; x<N; x++) {
    piold = pi[x];

    s = -0.5*dt*C*gb[x];

    // Nearly impenetrable!
    pi[x] = ( (1+s)*pi[x] 
	      + dt*( (xe[x]*xe[x]*(phi[nb[x][0]] - phi[x])
		      - xe[nb[x][1]]*xe[nb[x][1]]
		      *(phi[x] - phi[nb[x][1]]))
		     /(xc[x]*xc[x]*dx*dx) 
		     - Vdold[x] - 0.5*C*gb[x]
		     *(v[nb[x][1]]*(phi[x]-phi[nb[x][1]])
		       + v[x]*(phi[nb[x][0]] - phi[x]))/dx
		     ))/(1-s);
    
    pifull[x] = piold + 1.5*(pi[x] - piold);
 
  }

  // Move field
  for(x=0; x<N; x++) {
    phiold[x] = phi[x];
    phi[x] = phi[x] + dt*pi[x];
  }


  free(Vdold);
}



/* evolve_hydro(...)
 *
 * Straight from fortran. Original comment:
 *
 *   --  Does hydro except transport.  E  and  Z evolved.  v  and gb
 *       obtained.  Uses artificial viscosity  Q,  adjusted with  Cav.
 *   --  The order chosen is not unique, but is what I used in
 *       dissertation.  Should experiment with different orders,
 *       especially with position of field-fluid interaction,
 *       now placed first (otherwise order is from CW).
 *
 */
void evolve_hydro(double dt, double dx, double *kappa, double C, double Cav,
		  int N, double *phiold, double *phi, double *pi, double *p,
		  double *xe, double *xc,
		  double alpha, double gamma, double lambda,
		  double T0, double a,
		  double *E, double *Z, double *v, double *gb,
		  double *T, double *Q, int **nb) {

  int x;
  
  double *Vdmid = (double *)malloc(N*sizeof(double));
  double *gbold = (double *)malloc(N*sizeof(double));
  
  double *phiav = (double *)malloc(N*sizeof(double));
  double *dxphi = (double *)malloc(N*sizeof(double));
  
  double *gbv = (double *)malloc(N*sizeof(double));
  
  double gpi, dv, gv, inner, s;


  for(x=0; x<N; x++) {
    dxphi[x] = 0.5*(phiold[nb[x][0]] + phi[nb[x][0]] - phiold[x] - phi[x])/dx;
    phiav[x] = 0.5*(phiold[x] + phi[x]);
  }

  // Reflective BC's wrap
  dxphi[0] = 0.0;

  // Precompute potential
  Vdpot(N, alpha, gamma, lambda, T, T0, phiav, Vdmid);

  
  // Field-fluid interaction and artificial viscosity for 'E'
  for(x=0; x<N; x++) {

    Z[x] = Z[x] - dt*0.5*(C*(gb[x] + gb[nb[x][0]])
			  *(0.5*(pi[x] + pi[nb[x][0]]) + v[x]*dxphi[x] )
			  +(Vdmid[x] + Vdmid[nb[x][0]])) * dxphi[x];


    gpi = gb[x] * (pi[x] + 0.5*(v[nb[x][1]]*dxphi[nb[x][1]] + v[x]*dxphi[x]));

    E[x] = E[x] + dt*(C*gpi*gpi + gpi*Vdmid[x]);

    dv = v[x] - v[nb[x][1]];

    // D term deleted in the following four lines
    if (dv < 0)
      Q[x] = Cav*E[x]*dv*dv;
    else
      Q[x] = 0.0;

    E[x] = E[x] - dt*Q[x]*gb[x]*dv/dx;

  }



  // Pressure acceleration and artificial viscosity for 'Z'
  // The pressure used is advected old pressure
  for(x=0; x<N; x++) {

    Z[x] = Z[x] - dt*(p[nb[x][0]]- p[x] + Q[nb[x][0]] - Q[x])/dx;


    // D term deleted
    gv = 2.0*Z[x] / (kappa[x]*(E[x] + E[nb[x][0]]) );

    v[x] = gv/sqrt(1.0 + gv*gv);
  }

  // Boundary conditions require this, if we don't do "wrap"
  Z[0] = 0.0;
  v[0] = 0.0;


  // Obtain boost factor gb and pressure work on boost
  for(x=0; x<N; x++) {
    gbold[x] = gb[x];

    // D term deleted
    inner = (0.5*(Z[nb[x][1]] + Z[x])) / (kappa[x]*E[x]);
    gb[x] = sqrt(1.0 + inner*inner);

    s = (kappa[x] - 1)*(gb[x] - gbold[x])/(gb[x] + gbold[x]);

    E[x] = E[x]*(1-s)/(1+s);

  }


  // Pressure work on divergence
  // Original comment: "Like CW, use time-averaged  gb,  but new  v."
  for(x=0; x<N; x++)
    gbv[x] = v[x]*xe[x]*xe[x]
      *(gbold[x] + gb[x] + gb[nb[x][0]] + gb[nb[x][0]])/4.0;


  // Reflective BC's wrap
  gbv[0] = 0.0;

  for(x=0; x<N; x++) {
    s = (kappa[x] - 1.0) * dt/ dx
      * (gbv[x] - gbv[nb[x][1]]) / (gbold[x] + gb[x])/(xc[x]*xc[x]);

    E[x] = E[x]*(1-s)/(1+s);
  }


  // Clean up memory. Surely we don't need all these temporary arrays?
  free(Vdmid);
  free(gbold);
  free(phiav);
  free(dxphi);
  free(gbv);
}
