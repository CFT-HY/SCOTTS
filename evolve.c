/* evolve.c
 *
 * Everything to do with evolving fields and fluids with the
 * equations of motion.
 */
#include "hydro.h"


// straight from fortran
void evolve_backstep(hydro_fields f, int **nb, hydro_params p) {

  int x;

  for(x = 0; x < p.N; x++)
    f.pi[x] = f.pifull[x] - 0.5*(p.dt/(p.dx*p.dx))
      * ((f.phi[nb[x][0]] - f.phi[x])
	 - (f.phi[x] - f.phi[nb[x][1]]))
      + 0.125*(p.dt/(p.dx))*(p.dt/(p.dx))
      * ((f.pifull[nb[x][0]] - f.pifull[x]) 
	 - (f.pifull[x] - f.pifull[nb[x][1]]))
      + 0.5*p.dt*Vdf(p, f.T[x], f.phi[x] - 0.25*p.dt*f.pifull[x]);

}


// straight from fortran, except set v = 0 explicitly rather than in wrap
void evolve_field(hydro_fields f, int **nb, hydro_params p) {
  
  // Slow and not strictly necessary, see comments in eos.c
  double *Vdold = (double *)malloc(p.N*sizeof(double));

  int x;

  double piold, s;


  Vdpot(p, f.T, f.phi, Vdold);

  // Move conjugate momentum (leapfrog)
  for(x = 0; x < p.N; x++) {
    piold = f.pi[x];

    // first-order viscosity term
    s = -0.5*p.dt*p.C*f.gb[x];

    // Evolve momentum, badly
    /* pi[x] = ( (1+s)*pi[x]
	      + dt*( (xe[x]*xe[x]*(phi[nb[x][0]] - phi[x])
		      - xe[nb[x][1]]*xe[nb[x][1]]
		      *(phi[x] - phi[nb[x][1]]))
		     /(xc[x]*xc[x]*dx*dx) 
		     - Vdold[x] - 0.5*C*gb[x]
		     *(v[nb[x][1]]*(phi[x]-phi[nb[x][1]])
		       + v[x]*(phi[nb[x][0]] - phi[x]))/dx
		     ))/(1-s);
    */    


    // first order term 
    f.pi[x] = (1+s)*f.pi[x]/(1-s);

    // gradient term
    f.pi[x] = f.pi[x] 
      - 0.5*p.C*f.gb[x]*(f.v[nb[x][1]]*(f.phi[x]-f.phi[nb[x][1]])
			 + f.v[x]*(f.phi[nb[x][0]] - f.phi[x]))/p.dx;

    // splitting like this *seems* to conserve energy better than 
    // the original scheme above

    // vanilla scalar field gradient and potential terms
    // but notice funky centred derivative
    f.pi[x] = f.pi[x] 
      + p.dt*(((f.phi[nb[x][0]] - f.phi[x]) 
	       - (f.phi[x] - f.phi[nb[x][1]]))/(p.dx*p.dx)
	      - Vdold[x]);
    

    // pifull is (1.5*pi - 0.5*piold), not sure why
    f.pifull[x] = piold + 1.5*(f.pi[x] - piold);
 
  }

  // Move field (leapfrog)
  for(x = 0; x < p.N; x++) {
    f.phiold[x] = f.phi[x];
    f.phi[x] = f.phi[x] + p.dt*f.pi[x];
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
void evolve_hydro(hydro_fields f, int **nb, hydro_params p) {

  int x;
  
  double *Vdmid = (double *)malloc(p.N*sizeof(double));
  double *gbold = (double *)malloc(p.N*sizeof(double));
  
  double *phiav = (double *)malloc(p.N*sizeof(double));
  double *dxphi = (double *)malloc(p.N*sizeof(double));
  
  double *gbv = (double *)malloc(p.N*sizeof(double));

  double *Q = (double *)malloc(p.N*sizeof(double));
  
  double gpi, dv, gv, inner, s;


  // phi lives inside zones, Z lives on zone boundaries
  // phi is half a timestep ahead of Z (being zonal)
  // therefore the phi that Z sees is phiav, and the gradient of phi is dxphi
  for(x = 0; x < p.N; x++) {
    dxphi[x] = 0.5*(f.phiold[nb[x][0]] + f.phi[nb[x][0]] 
		    - f.phiold[x] - f.phi[x])/p.dx;
    phiav[x] = 0.5*(f.phiold[x] + f.phi[x]);
  }

  // Reflective BC's wrap
  dxphi[0] = 0.0;

  // Precompute potential on TIMESLICE
  Vdpot(p, f.T, phiav, Vdmid);

  
  // Field-fluid interaction and artificial viscosity for 'E'
  for(x = 0; x < p.N; x++) {

   
    // evolve Z (eq 10 of 9512202)
    f.Z[x] = f.Z[x] - p.dt*0.5*(p.C*(f.gb[x] + f.gb[nb[x][0]]) \
				// eta times boundary centred gamma
				*(0.5*(f.pi[x] + f.pi[nb[x][0]]) 
				   + f.v[x]*dxphi[x] )			\
				// boundary centred pi, v*grad(phi)
				+(Vdmid[x] + Vdmid[nb[x][0]])) * dxphi[x];
                                //  SPAT'Y recentred potential times grad(phi)



    // zone centred gamma times (already zone centred) pi
    // plus v*grad(phi) which is all zone centred
    gpi = f.gb[x]*(f.pi[x] + 0.5*(f.v[nb[x][1]]*dxphi[nb[x][1]] 
				  + f.v[x]*dxphi[x]));
    // in above, could change nb[x][1] to nb[x][0] from fortran
    // for both v and dxphi

    // evolve E (eq 9 of 9512202)
    f.E[x] = f.E[x] + p.dt*(p.C*gpi*gpi + gpi*Vdmid[x]);
    // first term is coupling to scalar field, second is potential

    dv = f.v[x] - f.v[nb[x][1]];

    // D term deleted in the following four lines
    if (dv < 0)
      Q[x] = p.Cav*f.E[x]*dv*dv;
    else
      Q[x] = 0.0;

    // Artificial viscosity term
    f.E[x] = f.E[x] - p.dt*Q[x]*f.gb[x]*dv/p.dx;

  }



  // Pressure acceleration (and artificial viscosity for 'Z')
  // W&M sec 2.4.12
  for(x = 0; x < p.N; x++) {

    // Q term related to artificial viscosity; otherwise W&M eq (2.103)
    f.Z[x] = f.Z[x] - p.dt*(f.p[nb[x][0]] - f.p[x] + Q[nb[x][0]] - Q[x])/p.dx;

    // update velocity v; denominator is W&M eq (2.85)
    // but note grid is Eulerian and D=0
    // then gv is the four-velocity (2.84)
    // (what we call kappa they call sigma, ish?)
    gv = 2.0*f.Z[x] / (f.kappa[x]*(f.E[x] + f.E[nb[x][0]]) );

    // three velocity
    f.v[x] = gv/sqrt(1.0 + gv*gv);
  }


  // Boundary conditions require this, if we don't do "wrap"
  f.Z[0] = 0.0;
  f.v[0] = 0.0;


  // Original comment: "Obtain boost factor gb and pressure work on boost"
  for(x = 0; x < p.N; x++) {
    gbold[x] = f.gb[x];

    // Calculate new zone-centred boost factor
    // this is just a repeat of what we did above...
    inner = (0.5*(f.Z[nb[x][1]] + f.Z[x])) / (f.kappa[x]*f.E[x]);

    // inner*inner should equal U^2 + U[nb[x][1]]^2
    // to give the zone centred gamma factor
    // instead this uses the arithmetic mean (cf W&M eq 2.88)
    f.gb[x] = sqrt(1.0 + inner*inner);

    // This is W&M eq (2.89), poor man's way of doing the power
    s = (f.kappa[x] - 1)*(f.gb[x] - gbold[x])/(f.gb[x] + gbold[x]);
    // sort of E*exp(-2.0*s)
    f.E[x] = f.E[x]*(1-s)/(1+s);

  }


  // Pressure work on divergence
  // Original comment: "Like CW, use time-averaged  gb,  but new  v."
  for(x = 0; x < p.N; x++) {

    // Velocity times area (boundary coord squared)
    //times poor man's volume gamma
    gbv[x] = f.v[x]
      *(f.gb[x] + f.gb[nb[x][0]])/2.0;
    // this is W*V, we calculate
  
      // previously (why gbold is here, no idea):
      //      *(gbold[x] + gb[x] + gb[nb[x][0]] + gb[nb[x][0]])/4.0;

  }


  // Reflective BC's wrap
  gbv[0] = 0.0;

  // (grad w)/dx
  for(x = 0; x < p.N; x++) {

    // W&M (2.92) and (2.93) combined: divergence
    s = (gbv[x] - gbv[nb[x][1]])/p.dx;
    // divide by zone volume times zone centred boost (previously gb[x] + gbold[x], not sure why)
    s = s/1.0;

    // by this stage s should have contributions from inside (2.91) and (2.93)
    // here we do the kappa multiplication and divide by 2 in preparation for the poor man's exponential
    // we also divide by the common zone centred gamma factor
    s =  s*(f.kappa[x] - 1.0) * p.dt/(2.0*f.gb[x]);

    // smells like W&M eq (2.93)
    // similar to E*exp(-2.0*s);
    f.E[x] = f.E[x]*(1-s)/(1+s);
  }


  // Clean up memory. Surely we don't need all these temporary arrays?
  free(Vdmid);
  free(gbold);
  free(phiav);
  free(dxphi);
  free(gbv);
  free(Q);
}
