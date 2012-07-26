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
    f.pi[x] = f.pifull[x] - 0.5*p.dt
      *(f.phi[nb[x][0]] + f.phi[nb[x][2]] + f.phi[nb[x][4]] 
	- 6.0*f.phi[x] 
	+ f.phi[nb[x][1]] + f.phi[nb[x][3]] + f.phi[nb[x][5]])/(p.dx*p.dx)
      + 0.125*p.dt*p.dt
      *(f.pifull[nb[x][0]] - 2.0*f.pifull[x] + f.pifull[nb[x][1]])/(p.dx*p.dx)
      + 0.5*p.dt*Vdf(p, f.T[x], f.phi[x] - 0.25*p.dt*f.pifull[x]);


  halo_field(f.pi, p);
}


// straight from fortran, except set v = 0 explicitly rather than in wrap
void evolve_field(hydro_fields f, int **nb, hydro_params p) {
  
  int x;

  double piold, s;

  // Move conjugate momentum (leapfrog)
  for(x = 0; x < p.N; x++) {
    piold = f.pi[x];

    // first-order viscosity term
    s = -0.5*p.dt*p.C*f.W[x];
    f.pi[x] = (1+s)*f.pi[x]/(1-s);

    // gradient term
    f.pi[x] = f.pi[x] - (
			 + 0.5*p.dt*p.C*f.W[x]*(f.V[0][x]*(f.phi[x] - f.phi[nb[x][1]])
						+ f.V[0][nb[x][0]]*(f.phi[nb[x][0]] - f.phi[x]))/p.dx
			 + 0.5*p.dt*p.C*f.W[x]*(f.V[1][x]*(f.phi[x] - f.phi[nb[x][3]])
			 			+ f.V[1][nb[x][2]]*(f.phi[nb[x][2]] - f.phi[x]))/p.dx
			 + 0.5*p.dt*p.C*f.W[x]*(f.V[2][x]*(f.phi[x] - f.phi[nb[x][5]])
						+ f.V[2][nb[x][4]]*(f.phi[nb[x][4]] - f.phi[x]))/p.dx
			 )/(1-s);

    // scalar field gradient and potential terms
    f.pi[x] = f.pi[x]
      + p.dt
      *((f.phi[nb[x][0]] + f.phi[nb[x][2]] + f.phi[nb[x][4]]
	 - 6.0*f.phi[x] 
	 + f.phi[nb[x][1]] + f.phi[nb[x][3]] + f.phi[nb[x][5]]
	 )/(p.dx*p.dx)
	- Vdf(p, f.T[x], f.phi[x]));
    

    // pifull is (1.5*pi - 0.5*piold), not sure why
    f.pifull[x] = piold + 1.5*(f.pi[x] - piold);
 
  }

  // Move field (leapfrog)
  for(x = 0; x < p.N; x++) {
    f.phiold[x] = f.phi[x];
    f.phi[x] = f.phi[x] + p.dt*f.pi[x];
  }

  halo_field(f.pi, p);
  halo_field(f.pifull, p);

  halo_field(f.phi, p);
  halo_field(f.phiold, p);

}




/* evolve_hydro(...)
 *
 * Straight from fortran. Original comment:
 *
 *   --  Does hydro except transport.  E  and  Z evolved.  v  and W
 *       obtained.  Uses artificial viscosity  Q,  adjusted with  Cav.
 *   --  The order chosen is not unique, but is what I used in
 *       dissertation.  Should experiment with different orders,
 *       especially with position of field-fluid interaction,
 *       now placed first (otherwise order is from CW).
 *
 */
void evolve_hydro(hydro_fields f, int **nb, hydro_params p) {

  int x;
  
  double *Vdmid = (double *)malloc(p.fieldN*sizeof(double));
  double *Wold = (double *)malloc(p.fieldN*sizeof(double));
  
  double *phiav = (double *)malloc(p.fieldN*sizeof(double));
  double **dxphi = (double **) malloc(3*sizeof(double *));
  dxphi[0] = (double *)malloc(p.fieldN*sizeof(double));
  dxphi[1] = (double *)malloc(p.fieldN*sizeof(double));
  dxphi[2] = (double *)malloc(p.fieldN*sizeof(double));
  

  double *Wfacex = (double *)malloc(p.fieldN*sizeof(double));
  double *Wfacey = (double *)malloc(p.fieldN*sizeof(double));
  double *Wfacez = (double *)malloc(p.fieldN*sizeof(double));

  
  double gpi, dv, gv, inner, s;

  double p_bar_x_plus, p_bar_x_minus, p_bar_y_plus, p_bar_y_minus, p_bar_z_plus, p_bar_z_minus;

  double utildex, utildey, utildez;
  double qx, qy, qz, gradv;
  double divv;

  double sigmabar;

  double ubarx, ubary, ubarz, W;


  for(x = 0; x < p.N; x++) {
    phiav[x] = 0.5*(f.phiold[x] + f.phi[x]);

    dxphi[0][x] =  0.5*(f.phiold[nb[x][0]] + f.phi[nb[x][0]]
			- f.phiold[x] - f.phi[x])/p.dx;

    dxphi[1][x] =  0.5*(f.phiold[nb[x][2]] + f.phi[nb[x][2]]
			- f.phiold[x] - f.phi[x])/p.dx;

    dxphi[2][x] =  0.5*(f.phiold[nb[x][4]] + f.phi[nb[x][4]]
			- f.phiold[x] - f.phi[x])/p.dx;

  }

  halo_field(dxphi[0], p);
  halo_field(dxphi[1], p);
  halo_field(dxphi[2], p);

  Vdpot(p, f.T, phiav, Vdmid);


  halo_field(Vdmid, p);

  for(x = 0; x < p.N; x++) {

    double vdnb = 0.125*(Vdmid[nb[x][1]] + Vdmid[x] + Vdmid[nb[x][3]] + Vdmid[nb[x][DIR_13]]
			 + Vdmid[nb[x][5]] + Vdmid[nb[x][DIR_35]] + Vdmid[nb[x][DIR_15]] + Vdmid[nb[x][DIR_135]]);


    double pinb = 0.125*(f.pi[nb[x][1]] + f.pi[x] + f.pi[nb[x][3]] + f.pi[nb[x][DIR_13]]
			 + f.pi[nb[x][5]] + f.pi[nb[x][DIR_35]] + f.pi[nb[x][DIR_15]] + f.pi[nb[x][DIR_135]]);


    double wnb = 0.125*(f.W[nb[x][1]] + f.W[x] + f.W[nb[x][3]] + f.W[nb[x][DIR_13]]
			 + f.W[nb[x][5]] + f.W[nb[x][DIR_35]] + f.W[nb[x][DIR_15]] + f.W[nb[x][DIR_135]]);

   
    double dxphinb0 = 0.25*(dxphi[0][nb[x][1]] + dxphi[0][nb[x][DIR_13]]
			   + dxphi[0][nb[x][DIR_15]] + dxphi[0][nb[x][DIR_135]]);


    double dxphinb1 = 0.25*(dxphi[1][nb[x][3]] + dxphi[1][nb[x][DIR_13]]
			   + dxphi[1][nb[x][DIR_35]] + dxphi[1][nb[x][DIR_135]]);

    double dxphinb2 = 0.25*(dxphi[2][nb[x][5]] + dxphi[2][nb[x][DIR_15]]
			   + dxphi[2][nb[x][DIR_35]] + dxphi[2][nb[x][DIR_135]]);

   
    
    f.Z[0][x] = f.Z[0][x] - p.dt*p.C*wnb*pinb*dxphinb0;
    f.Z[1][x] = f.Z[1][x] - p.dt*p.C*wnb*pinb*dxphinb1;
    f.Z[2][x] = f.Z[2][x] - p.dt*p.C*wnb*pinb*dxphinb2;
    
    f.Z[1][x] = f.Z[1][x] - p.dt*p.C*wnb*f.V[1][x]*dxphinb1*dxphinb1;
    f.Z[0][x] = f.Z[0][x] - p.dt*p.C*wnb*f.V[0][x]*dxphinb0*dxphinb0;    
    f.Z[2][x] = f.Z[2][x] - p.dt*p.C*wnb*f.V[2][x]*dxphinb2*dxphinb2;

    f.Z[0][x] = f.Z[0][x] - p.dt*vdnb*dxphinb0;     
    f.Z[1][x] = f.Z[1][x] - p.dt*vdnb*dxphinb1;
    f.Z[2][x] = f.Z[2][x] - p.dt*vdnb*dxphinb2;
				       
    gpi = f.W[x]*(f.pi[x] + 0.5*(f.V[0][x]*dxphi[0][nb[x][1]]
				 + f.V[0][nb[x][0]]*dxphi[0][x]));
    f.E[x] = f.E[x] + p.dt*(p.C*gpi*gpi + gpi*Vdmid[x]);


    gpi = f.W[x]*(f.pi[x] + 0.5*(f.V[1][x]*dxphi[1][nb[x][3]]
				 + f.V[1][nb[x][2]]*dxphi[1][x]));
    f.E[x] = f.E[x] + p.dt*(p.C*gpi*gpi + gpi*Vdmid[x]);

    gpi = f.W[x]*(f.pi[x] + 0.5*(f.V[2][x]*dxphi[2][nb[x][5]]
				 + f.V[2][nb[x][4]]*dxphi[2][x]));
    f.E[x] = f.E[x] + p.dt*(p.C*gpi*gpi + gpi*Vdmid[x]);
    
  }


  // halo_field(f.Z[0], p);
  // halo_field(f.Z[1], p);
  // halo_field(f.Z[2], p);
  //  halo_field(f.E, p);



  // Pressure acceleration
  // W&M sec 2.4.12, 3.5.1, DONE
  for(x = 0; x < p.N; x++) {

    // equation (3.68)
    // repeated calculation of these unnecessary, so split it out later
    p_bar_x_plus = (f.p[x]
		    + f.p[nb[x][3]] 
		    + f.p[nb[x][5]]
		    + f.p[nb[x][DIR_35]]
		    )/4.0;
    
    p_bar_x_minus = (f.p[nb[x][1]] 
		     + f.p[nb[x][DIR_13]]
		     + f.p[nb[x][DIR_15]]
		     + f.p[nb[x][DIR_135]]
		     )/4.0;
    
    p_bar_y_plus = (f.p[x]
		    + f.p[nb[x][1]] 
		    + f.p[nb[x][5]]
		    + f.p[nb[x][DIR_15]]
		    )/4.0;
    
    
    p_bar_y_minus = (f.p[nb[x][3]]
		     + f.p[nb[x][DIR_13]] 
		     + f.p[nb[x][DIR_35]]
		     + f.p[nb[x][DIR_135]]
		     )/4.0;
    
    
    p_bar_z_plus = (f.p[x]
		    + f.p[nb[x][1]] 
		    + f.p[nb[x][3]]
		    + f.p[nb[x][DIR_13]]
		    )/4.0;
    
    
    p_bar_z_minus = (f.p[nb[x][5]]
		     + f.p[nb[x][DIR_15]] 
		     + f.p[nb[x][DIR_35]]
		     + f.p[nb[x][DIR_135]]
		     )/4.0;
    
    
	
    
    f.Z[0][x] = f.Z[0][x] - (p_bar_x_plus - p_bar_x_minus)*p.dt/p.dx;
    
    f.Z[1][x] = f.Z[1][x] - (p_bar_y_plus - p_bar_y_minus)*p.dt/p.dx;
    
    f.Z[2][x] = f.Z[2][x] - (p_bar_z_plus - p_bar_z_minus)*p.dt/p.dx;

  }


  //  halo_field(f.Z[0], p);
  //  halo_field(f.Z[1], p);
  //  halo_field(f.Z[2], p);



  // update velocity v; denominator is W&M eq (2.85)
  // but note grid is Eulerian and D=0
  // then gv is the four-velocity (2.84)
  // (what we call kappa they call sigma, ish?)
  //
  // Section 3.4.5, equations 3.5.7, 3.5.8
  for(x = 0; x < p.N; x++) {

    /*
    if(x == 0)
      fprintf(stderr,"old Ux = %lf\n", f.U[0][1]);
    */

    sigmabar = (f.kappa[x]*f.E[x] 
		+ f.kappa[nb[x][1]]*f.E[nb[x][1]]
		+ f.kappa[nb[x][DIR_13]]*f.E[nb[x][DIR_13]]
		+ f.kappa[nb[x][DIR_135]]*f.E[nb[x][DIR_135]]
		+ f.kappa[nb[x][DIR_35]]*f.E[nb[x][DIR_35]]
		+ f.kappa[nb[x][DIR_15]]*f.E[nb[x][DIR_15]]
		+ f.kappa[nb[x][3]]*f.E[nb[x][3]]
		+ f.kappa[nb[x][5]]*f.E[nb[x][5]]
		)/8.0;

   
    f.U[0][x] = f.Z[0][x]/sigmabar;
    f.U[1][x] = f.Z[1][x]/sigmabar;
    f.U[2][x] = f.Z[2][x]/sigmabar;

    //    fprintf(stderr,"%lf %lf\n", f.Z[0][x], sigmabar);
    /*
    if(x == 0)
      fprintf(stderr,"Ux = %lf\n", f.U[0][1]);
    */
  }

  halo_field(f.U[0], p);
  halo_field(f.U[1], p);
  halo_field(f.U[2], p);


  // section 3.4.5 continued
  // face-centred x-cpt of 3-velocity
  for(x = 0; x < p.N; x++) {

    ubarx = (f.U[0][x]
	     + f.U[0][nb[x][2]] 
	     + f.U[0][nb[x][4]]
	     + f.U[0][nb[x][DIR_24]]
	     )/4.0;
    
    ubary = (f.U[1][x]
	     + f.U[1][nb[x][2]]
	     + f.U[1][nb[x][4]]
	     + f.U[1][nb[x][DIR_24]]
	     )/4.0;
    
    ubarz = (f.U[2][x]
	     + f.U[2][nb[x][2]]
	     + f.U[2][nb[x][4]]
	     + f.U[2][nb[x][DIR_24]]
	     )/4.0;
    
    Wfacex[x] = sqrt(1.0 + ubarx*ubarx + ubary*ubary + ubarz*ubarz);
    
    f.V[0][x] = ubarx/Wfacex[x];
  }




  // y-cpt
  for(x = 0; x < p.N; x++) {

	ubarx = (f.U[0][x]
		 + f.U[0][nb[x][0]] 
		 + f.U[0][nb[x][4]]
		 + f.U[0][nb[x][DIR_04]]
		 )/4.0;

	ubary = (f.U[1][x]
		 + f.U[1][nb[x][0]] 
		 + f.U[1][nb[x][4]]
		 + f.U[1][nb[x][DIR_04]]
		 )/4.0;

	ubarz = (f.U[2][x]
		 + f.U[2][nb[x][0]] 
		 + f.U[2][nb[x][4]]
		 + f.U[2][nb[x][DIR_04]]
		 )/4.0;

	Wfacey[x] = sqrt(1.0 + ubarx*ubarx + ubary*ubary + ubarz*ubarz);

	f.V[1][x] = ubary/Wfacey[x];

  }
  // z-cpt
  for(x = 0; x < p.N; x++) {

    ubarx = (f.U[0][x]
	     + f.U[0][nb[x][2]] 
	     + f.U[0][nb[x][0]]
	     + f.U[0][nb[x][DIR_02]]
	     )/4.0;
    
    ubary = (f.U[1][x]
	     + f.U[1][nb[x][2]] 
	     + f.U[1][nb[x][0]]
	     + f.U[1][nb[x][DIR_02]]
	     )/4.0;
    
    ubarz = (f.U[2][x]
	     + f.U[2][nb[x][2]] 
	     + f.U[2][nb[x][0]]
	     + f.U[2][nb[x][DIR_02]]
	     )/4.0;
    
    Wfacez[x] = sqrt(1.0 
		     + ubarx*ubarx
		     + ubary*ubary
		     + ubarz*ubarz);
    
    f.V[2][x] = ubarz/Wfacez[x];
  }


  halo_field(f.V[0], p);
  halo_field(f.V[1], p);
  halo_field(f.V[2], p);



  // End section 3.4.5

  //  fprintf(stderr,"Vx[1] = %lf\n", f.V[0][1]);
  //  getchar();

  // Section 3.4.6
  for(x = 0; x < p.N; x++) {

    utildex = (f.U[0][x] 
	       + f.U[0][nb[x][0]]
	       + f.U[0][nb[x][2]]
	       + f.U[0][nb[x][4]]
	       + f.U[0][nb[x][DIR_04]]
	       + f.U[0][nb[x][DIR_24]]
	       + f.U[0][nb[x][DIR_02]]
	       + f.U[0][nb[x][DIR_024]]
	       )/8.0;
    
    utildey = (f.U[1][x] 
	       + f.U[1][nb[x][0]]
	       + f.U[1][nb[x][2]]
	       + f.U[1][nb[x][4]]
	       + f.U[1][nb[x][DIR_04]]
	       + f.U[1][nb[x][DIR_24]]
	       + f.U[1][nb[x][DIR_02]]
	       + f.U[1][nb[x][DIR_024]]
	       )/8.0;
    
    utildez = (f.U[2][x] 
	       + f.U[2][nb[x][0]]
	       + f.U[2][nb[x][2]]
	       + f.U[2][nb[x][4]]
	       + f.U[2][nb[x][DIR_04]]
	       + f.U[2][nb[x][DIR_24]]
	       + f.U[2][nb[x][DIR_02]]
	       + f.U[2][nb[x][DIR_024]]
	       )/8.0;


    //    if(x < 3 || x > p.N-3)
    //      fprintf(stderr,"utildex %d = %lf\n", x, utildex);
     
    Wold[x] = f.W[x];
    f.W[x] = sqrt(1.0 
		  + utildex*utildex 
		  + utildey*utildey 
		  + utildez*utildez);

    
    //	s = (f.kappa[iix(x,y,z,p)] - 1)*(f.W[iix(x,y,z,p)] - Wold[iix(x,y,z,p)])
    //	  /(f.W[iix(x,y,z,p)] + Wold[iix(x,y,z,p)]);
    //	f.E[iix(x,y,z,p)] = f.E[iix(x,y,z,p)]*(1-s)/(1+s);
    // sort of E*exp(-2.0*s)
    
    
    // Top of p90 ch 3
    f.E[x] = f.E[x]*pow(Wold[x]/f.W[x],f.kappa[x]-1.0);





  }



  //  halo_field(f.W, p);

  halo_field(Wfacex, p);
  halo_field(Wfacey, p);
  halo_field(Wfacez, p);

  // Section 3.4.4 -- pressure terms
  for(x = 0; x < p.N; x++) {
	qx = (f.V[0][nb[x][0]] - f.V[0][x])/p.dx;
	qy = (f.V[1][nb[x][2]] - f.V[1][x])/p.dx;
	qz = (f.V[2][nb[x][4]] - f.V[2][x])/p.dx;

	divv = qx + qy + qz;

		
	// s = 0.5*(f.kappa[x] - 1.0)*divv*p.dt/2.0;
	// f.E[x] = f.E[x]*(1-s)/(1+s);
	

	// Is it 1.0 or 0.5? -- not the major problem I think
	// also note that dx=1.0 still get problems
	// so it is not some factor of dx somewhere
	// nor is it this particular expression:
	f.E[x] = f.E[x]*exp(-1.0*(f.kappa[x] - 1.0)*divv*p.dt);

  }




  //  halo_field(f.E, p);



  // Last bit of 3.4.6
  for(x = 0; x < p.N; x++) {
    qx = (f.V[0][x] + f.V[0][nb[x][0]])
      *(Wfacex[nb[x][0]] - Wfacex[x])*p.dt/(2.0*p.dx);
    
    qy = (f.V[1][x] + f.V[1][nb[x][2]])
      *(Wfacey[nb[x][2]] - Wfacey[x])*p.dt/(2.0*p.dx);
    
    qz = (f.V[2][x] + f.V[2][nb[x][4]])
      *(Wfacez[nb[x][4]] - Wfacez[x])*p.dt/(2.0*p.dx);
    
    gradv = (qx+qy+qz)/f.W[x];
    
    //	s = (f.kappa[x] - 1.0)*gradv/2.0;
    //	f.E[iix(x,y,z,p)] = f.E[iix(x,y,z,p)]*(1-s)/(1+s);
    f.E[x] = f.E[x]*exp(-1.0*(f.kappa[x]-1.0)*gradv);
       
  }


  halo_field(f.E, p);



  // Clean up memory. Surely we don't need all these temporary arrays?


  free(Vdmid);
  free(Wold);

  free(phiav);

  free(dxphi[0]);
  free(dxphi[1]);
  free(dxphi[2]);
  free(dxphi);

  free(Wfacex);
  free(Wfacey); 
  free(Wfacez);


}

void artificial_viscosity(hydro_fields f, int **nb, hydro_params p) {
  int x;

  double Ux, Uxp, Uy, Uyp, Uz, Uzp;

  double *zoneDiv = (double *)malloc(p.N*sizeof(double));  

  double Q = 1.0;
  double L = p.dt;

  double Umin, Umax, Udix, Udiv;
  double *c = (double *)malloc(p.N*sizeof(double));  
  double *qx = (double *)malloc(p.N*sizeof(double));  
  double *q = (double *)malloc(p.N*sizeof(double));  
  double *tx = (double *)malloc(p.N*sizeof(double));  
  double g, F;
  double delta, deltam, deltav;
  double a,b,C;
  double qt;
  double A = 0.25;


  double cs = 0.025;


  // Zone centred divergence (W&M p93, bottom)
  for(x=0;x<p.N;x++) {

    Ux = f.U[0][x] + f.U[0][nb[x][2]] + f.U[0][nb[x][4]] + f.U[0][nb[x][DIR_24]];
    Uxp = f.U[0][nb[x][0]] + f.U[0][nb[x][DIR_02]] + f.U[0][nb[x][DIR_04]]
      + f.U[0][nb[x][DIR_024]];

    Uy = f.U[1][x] + f.U[1][nb[x][0]] + f.U[1][nb[x][4]] + f.U[1][nb[x][DIR_04]];
    Uyp = f.U[1][nb[x][2]] + f.U[1][nb[x][DIR_02]] + f.U[1][nb[x][DIR_24]]
      + f.U[1][nb[x][DIR_024]];

    Uz = f.U[2][x] + f.U[2][nb[x][0]] + f.U[2][nb[x][2]] + f.U[2][nb[x][DIR_02]];
    Uzp = f.U[2][nb[x][4]] + f.U[2][nb[x][DIR_04]] + f.U[2][nb[x][DIR_24]]
      + f.U[2][nb[x][DIR_024]];

    zoneDiv[x] =  0.25*(Uxp - Ux + Uyp - Uy + Uzp - Uz)/p.dx;

  }




  // x-direction art visc
  for(x=0; x<p.N; x++) {
    Umin = minof3(f.U[0][nb[x][1]],f.U[0][x],f.U[0][nb[x][0]]);
    Umax = maxof3(f.U[0][nb[x][1]],f.U[0][x],f.U[0][nb[x][0]]);
    Udix = minof2(Umax-f.U[0][x],f.U[0][x]-Umin);

    if((f.U[0][nb[x][0]]-f.U[0][x])*(f.U[0][x]-f.U[0][nb[x][1]]) > 0) {
      Udiv = Udix;
    } else {
      Udiv = 0;
    }

    c[x] = 0.5*(f.U[0][x] - f.U[0][nb[x][0]] + f.U[0][nb[x][1]] - f.U[0][x])/p.dx;
  }

  for(x=0; x<p.N; x++) {
    g = sqrt(1.0 + f.U[0][x]*f.U[0][x] + f.U[1][x]*f.U[1][x] + f.U[2][x]*f.U[2][x]);
    F = f.kappa[x]*f.E[x] + f.kappa[nb[x][5]]*f.E[nb[x][5]]
      + f.kappa[nb[x][3]]*f.E[nb[x][3]] + f.kappa[nb[x][DIR_35]]*f.E[nb[x][DIR_35]];
    delta = f.U[0][nb[x][0]] - f.U[0][x] - 0.5*p.dx*(c[x] + c[nb[x][0]])/g;



    if(delta*(f.U[0][nb[x][0]]-f.U[0][x]) > 0.0)
      deltam = delta;
    else
      deltam = 0.0;

    deltav = minof2(deltam,0.0);

    qx[x] = F*(Q*deltav - L*cs);
    tx[x] = F*deltav*(f.V[0][x] + f.V[0][nb[x][0]]);
  }


  for(x=0;x<p.N;x++) {
    if(zoneDiv[x] > 0)
      q[x] = 0.0;
    else
      q[x] = 0.25*(qx[x] + qx[nb[x][4]] + qx[nb[x][2]] + qx[nb[x][DIR_24]]);

    a = tx[x] + tx[nb[x][2]] + tx[nb[x][4]] + tx[nb[x][DIR_24]];
    b = 0.0;
    C = 0.0;

    qt = q[x] - (a + b + C)*A/32.0;

    if(zoneDiv[x] < 0)
      qt = 0.0;

    f.E[x] = f.E[x] + (q[x] - qt)*zoneDiv[x];

  }


  free(zoneDiv);
  free(c);
  free(qx);
  free(q);
  free(tx);

}
