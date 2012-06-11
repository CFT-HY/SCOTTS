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
    f.pi[x] = f.pi[x] 
      - 0.5*p.C*f.W[x]*(f.Vx[nb[x][1]]*(f.phi[x]-f.phi[nb[x][1]])
			+ f.Vx[x]*(f.phi[nb[x][0]] - f.phi[x]))/p.dx;    

    // scalar field gradient and potential terms
    f.pi[x] = f.pi[x]
      + p.dt
      *((f.phi[nb[x][0]] + f.phi[nb[x][2]] + f.phi[nb[x][4]]
	 - 6.0*f.phi[x] 
	 + f.phi[nb[x][1]] + f.phi[nb[x][3]] + f.phi[nb[x][5]])/(p.dx*p.dx)
	- Vdf(p, f.T[x], f.phi[x]));
    

    // pifull is (1.5*pi - 0.5*piold), not sure why
    f.pifull[x] = piold + 1.5*(f.pi[x] - piold);
 
  }

  // Move field (leapfrog)
  for(x = 0; x < p.N; x++) {
    f.phiold[x] = f.phi[x];
    f.phi[x] = f.phi[x] + p.dt*f.pi[x];
  }

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
  
  double *Vdmid = (double *)malloc(p.N*sizeof(double));
  double *Wold = (double *)malloc(p.N*sizeof(double));
  
  double *phiav = (double *)malloc(p.N*sizeof(double));
  double *dxphi = (double *)malloc(p.N*sizeof(double));
  
  double *Wv = (double *)malloc(p.N*sizeof(double));


  double *Wfacex = (double *)malloc(p.N*sizeof(double));
  double *Wfacey = (double *)malloc(p.N*sizeof(double));
  double *Wfacez = (double *)malloc(p.N*sizeof(double));

  double *Q = (double *)malloc(p.N*sizeof(double));
  
  double gpi, dv, gv, inner, s;

  double p_bar_x_plus, p_bar_x_minus, p_bar_y_plus, p_bar_y_minus, p_bar_z_plus, p_bar_z_minus;

  double utildex, utildey, utildez;
  double qx, qy, qz, gradv;
  double divv;

  double sigmabar;

  double ubarx, ubary, ubarz, W;

  Vdpot(p, f.T, f.phi, Vdmid);

  for(x = 0; x < p.N; x++) {
    f.Zx[x] = f.Zx[x] - p.dt*0.5*(p.C*(f.W[x] + f.W[nb[x][0]])
				  *(0.5*(f.pi[x]+f.pi[nb[x][0]])
				    + f.Vx[x]*(f.phi[nb[x][0]] - f.phi[x])/p.dx)
				  + (Vdmid[x] + Vdmid[nb[x][0]]))
				       *(f.phi[nb[x][0]] - f.phi[x])/p.dx;
    
    /*
      Z(j) = Z(j) - dt*0.5* ( C*(gb(j)+gb(j+1)) *
      .      ( 0.5*(pi(j)+pi(j+1)) + v(j)*dxfi(j) )
      .      +(Vdmid(j)+Vdmid(j+1)) ) * dxfi(j)
    */


    //    gpi = gb(j) * ( pi(j) + 0.5*(v(j-1)*dxfi(j-1)+v(j)*dxfi(j)) )
    //      E(j) = E(j) + dt * (C*gpi**2 + Vdmid(j)*gpi)

    gpi = f.W[x]*(f.pi[x] + 0.5*(f.Vx[nb[x][1]]*(f.phi[x]-f.phi[nb[x][1]])/p.dx
				 + f.Vx[x]*(f.phi[nb[x][0]] - f.phi[x])/p.dx));
    f.E[x] = f.E[x] + p.dt*(p.C*gpi*gpi + gpi*Vdmid[x]);

  }

  // Pressure acceleration
  // W&M sec 2.4.12, 3.5.1, DONE
  for(x = 0; x < p.N; x++) {

    // equation (3.68)
    // repeated calculation of these unnecessary, so split it out later
    p_bar_x_plus = (f.p[x]
		    + f.p[nb[x][3]] 
		    + f.p[nb[x][5]]
		    + f.p[nb[nb[x][3]][5]]
		    )/4.0;
    
    p_bar_x_minus = (f.p[nb[x][1]] 
		     + f.p[nb[nb[x][1]][3]]
		     + f.p[nb[nb[x][1]][5]]
		     + f.p[nb[nb[nb[x][1]][3]][5]]
		     )/4.0;
    
    p_bar_y_plus = (f.p[x]
		    + f.p[nb[x][1]] 
		    + f.p[nb[x][5]]
		    + f.p[nb[nb[x][1]][5]]
		    )/4.0;
    
    
    p_bar_y_minus = (f.p[nb[x][3]]
		     + f.p[nb[nb[x][1]][3]] 
		     + f.p[nb[nb[x][3]][5]]
		     + f.p[nb[nb[nb[x][1]][3]][5]]
		     )/4.0;
    
    
    p_bar_z_plus = (f.p[x]
		    + f.p[nb[x][1]] 
		    + f.p[nb[x][3]]
		    + f.p[nb[nb[x][1]][3]]
		    )/4.0;
    
    
    p_bar_z_minus = (f.p[nb[x][5]]
		     + f.p[nb[nb[x][1]][5]] 
		     + f.p[nb[nb[x][3]][5]]
		     + f.p[nb[nb[nb[x][1]][3]][5]]
		     )/4.0;
    
    
	
    
    f.Zx[x] = f.Zx[x] - (p_bar_x_plus - p_bar_x_minus)*p.dt/p.dx;
    // No difference:
    //    f.Zx[x] = f.Zx[x] - (f.p[x] - f.p[nb[x][1]])*p.dt/p.dx;
    
    f.Zy[x] = f.Zy[x] - (p_bar_y_plus - p_bar_y_minus)*p.dt/p.dx;
    
    f.Zz[x] = f.Zz[x] - (p_bar_z_plus - p_bar_z_minus)*p.dt/p.dx;

  }

  //  fprintf(stderr,"p[0] = %lf, p[1] = %lf\n", f.p[0], f.p[1]);
  //  fprintf(stderr,"Z[0] = %lf, Z[1] = %lf\n", f.Zx[0], f.Zx[1]);

  // update velocity v; denominator is W&M eq (2.85)
  // but note grid is Eulerian and D=0
  // then gv is the four-velocity (2.84)
  // (what we call kappa they call sigma, ish?)
  //
  // Section 3.4.5, equations 3.5.7, 3.5.8
  for(x = 0; x < p.N; x++) {

    f.Ux[x]  = 2.0*f.Zx[x] / (f.kappa[x]*(f.E[x] + f.E[nb[x][0]]) );

    /*
    if(x == 0)
      fprintf(stderr,"old Ux = %lf\n", f.Ux[1]);
    */

    sigmabar = (f.kappa[x]*f.E[x] 
		+ f.kappa[nb[x][1]]*f.E[nb[x][1]]
		+ f.kappa[nb[nb[x][1]][3]]*f.E[nb[nb[x][1]][3]]
		+ f.kappa[nb[nb[nb[x][1]][3]][5]]*f.E[nb[nb[nb[x][1]][3]][5]]
		+ f.kappa[nb[nb[x][3]][5]]*f.E[nb[nb[x][3]][5]]
		+ f.kappa[nb[nb[x][1]][5]]*f.E[nb[nb[x][1]][5]]
		+ f.kappa[nb[x][3]]*f.E[nb[x][3]]
		+ f.kappa[nb[x][5]]*f.E[nb[x][5]]
		)/8.0;

   
    f.Ux[x] = f.Zx[x]/sigmabar;
    f.Uy[x] = f.Zy[x]/sigmabar;
    f.Uz[x] = f.Zz[x]/sigmabar;

    //    fprintf(stderr,"%lf %lf\n", f.Zx[x], sigmabar);
    /*
    if(x == 0)
      fprintf(stderr,"Ux = %lf\n", f.Ux[1]);
    */
  }


  // section 3.4.5 continued
  // face-centred x-cpt of 3-velocity
  for(x = 0; x < p.N; x++) {

    ubarx = (f.Ux[x]
	     + f.Ux[nb[x][2]] 
	     + f.Ux[nb[x][4]]
	     + f.Ux[nb[nb[x][2]][4]]
	     )/4.0;
    
    ubary = (f.Uy[x]
	     + f.Uy[nb[x][2]]
	     + f.Uy[nb[x][4]]
	     + f.Uy[nb[nb[x][2]][4]]
	     )/4.0;
    
    ubarz = (f.Uz[x]
	     + f.Uz[nb[x][2]]
	     + f.Uz[nb[x][4]]
	     + f.Uz[nb[nb[x][2]][4]]
	     )/4.0;
    
    Wfacex[x] = sqrt(1.0 + ubarx*ubarx + ubary*ubary + ubarz*ubarz);
    
    f.Vx[x] = ubarx/Wfacex[x];
  }
  // y-cpt
  for(x = 0; x < p.N; x++) {

	ubarx = (f.Ux[x]
		 + f.Ux[nb[x][0]] 
		 + f.Ux[nb[x][4]]
		 + f.Ux[nb[nb[x][0]][4]]
		 )/4.0;

	ubary = (f.Uy[x]
		 + f.Uy[nb[x][0]] 
		 + f.Uy[nb[x][4]]
		 + f.Uy[nb[nb[x][0]][4]]
		 )/4.0;

	ubarz = (f.Uz[x]
		 + f.Uz[nb[x][0]] 
		 + f.Uz[nb[x][4]]
		 + f.Uz[nb[nb[x][0]][4]]
		 )/4.0;

	Wfacey[x] = sqrt(1.0 + ubarx*ubarx + ubary*ubary + ubarz*ubarz);

	f.Vy[x] = ubary/Wfacey[x];

  }
  // z-cpt
  for(x = 0; x < p.N; x++) {

    ubarx = (f.Ux[x]
	     + f.Ux[nb[x][2]] 
	     + f.Ux[nb[x][0]]
	     + f.Ux[nb[nb[x][0]][2]]
	     )/4.0;
    
    ubary = (f.Uy[x]
	     + f.Uy[nb[x][2]] 
	     + f.Uy[nb[x][0]]
	     + f.Uy[nb[nb[x][0]][2]]
	     )/4.0;
    
    ubarz = (f.Uz[x]
	     + f.Uz[nb[x][2]] 
	     + f.Uz[nb[x][0]]
	     + f.Uz[nb[nb[x][0]][2]]
	     )/4.0;
    
    Wfacez[x] = sqrt(1.0 
		     + ubarx*ubarx
		     + ubary*ubary
		     + ubarz*ubarz);
    
    f.Vz[x] = ubarz/Wfacez[x];
  }

  // End section 3.4.5

  //  fprintf(stderr,"Vx[1] = %lf\n", f.Vx[1]);
  //  getchar();


  // Section 3.4.6
  for(x = 0; x < p.N; x++) {

    utildex = (f.Ux[x] 
	       + f.Ux[nb[x][0]]
	       + f.Ux[nb[x][2]]
	       + f.Ux[nb[x][4]]
	       + f.Ux[nb[nb[x][0]][4]]
	       + f.Ux[nb[nb[x][2]][4]]
	       + f.Ux[nb[nb[x][0]][2]]
	       + f.Ux[nb[nb[nb[x][0]][2]][4]]
	       )/8.0;
    
    utildey = (f.Uy[x] 
	       + f.Uy[nb[x][0]]
	       + f.Uy[nb[x][2]]
	       + f.Uy[nb[x][4]]
	       + f.Uy[nb[nb[x][0]][4]]
	       + f.Uy[nb[nb[x][2]][4]]
	       + f.Uy[nb[nb[x][0]][2]]
	       + f.Uy[nb[nb[nb[x][0]][2]][4]]
	       )/8.0;
    
    utildez = (f.Uz[x] 
	       + f.Uz[nb[x][0]]
	       + f.Uz[nb[x][2]]
	       + f.Uz[nb[x][4]]
	       + f.Uz[nb[nb[x][0]][4]]
	       + f.Uz[nb[nb[x][2]][4]]
	       + f.Uz[nb[nb[x][0]][2]]
	       + f.Uz[nb[nb[nb[x][0]][2]][4]]
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


  /*
  fprintf(stderr,"N-2 powarg = %lf/%lf\n", Wold[p.N-2],f.W[p.N-2]);
  fprintf(stderr,"N-1 powarg = %lf/%lf\n", Wold[p.N-1],f.W[p.N-1]);
  fprintf(stderr,"0 powarg = %lf/%lf\n", Wold[0],f.W[0]);
  fprintf(stderr,"1 powarg = %lf/%lf\n", Wold[1],f.W[1]);
  fprintf(stderr,"2 powarg = %lf/%lf\n", Wold[2],f.W[2]);

  fprintf(stderr,"powe! E[0] = %lf\n", f.E[0]);
  */

  // Section 3.4.4 -- pressure terms
  for(x = 0; x < p.N; x++) {
	qx = (f.Vx[nb[x][0]] - f.Vx[x])/p.dx;
	qy = (f.Vy[nb[x][2]] - f.Vy[x])/p.dx;
	qz = (f.Vz[nb[x][4]] - f.Vz[x])/p.dx;

	divv = qx + qy + qz;

		
	// s = 0.5*(f.kappa[x] - 1.0)*divv*p.dt/2.0;
	// f.E[x] = f.E[x]*(1-s)/(1+s);
	

	// Is it 1.0 or 0.5? -- not the major problem I think
	// also note that dx=1.0 still get problems
	// so it is not some factor of dx somewhere
	// nor is it this particular expression:
	f.E[x] = f.E[x]*exp(-1.0*(f.kappa[x] - 1.0)*divv*p.dt);

  }







  // Last bit of 3.4.6
  for(x = 0; x < p.N; x++) {
    qx = (f.Vx[x] + f.Vx[nb[x][0]])
      *(Wfacex[nb[x][0]] - Wfacex[x])*p.dt/(2.0*p.dx);
    
    qy = (f.Vy[x] + f.Vy[nb[x][2]])
      *(Wfacey[nb[x][2]] - Wfacey[x])*p.dt/(2.0*p.dx);
    
    qz = (f.Vz[x] + f.Vz[nb[x][4]])
      *(Wfacez[nb[x][4]] - Wfacez[x])*p.dt/(2.0*p.dx);
    
    gradv = (qx+qy+qz)/f.W[x];
    
    //	s = (f.kappa[x] - 1.0)*gradv/2.0;
    //	f.E[iix(x,y,z,p)] = f.E[iix(x,y,z,p)]*(1-s)/(1+s);
    f.E[x] = f.E[x]*exp(-1.0*(f.kappa[x]-1.0)*gradv);
       
  }



  // Clean up memory. Surely we don't need all these temporary arrays?


  free(Vdmid);
  free(Wold);
  free(Wfacex);
  free(Wfacey); 
  free(Wfacez);
  free(Wv);
  free(phiav);
  free(dxphi);
  free(Q);
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

    Ux = f.Ux[x] + f.Ux[nb[x][2]] + f.Ux[nb[x][4]] + f.Ux[nb[nb[x][2]][4]];
    Uxp = f.Ux[nb[x][0]] + f.Ux[nb[nb[x][2]][0]] + f.Ux[nb[nb[x][4]][0]]
      + f.Ux[nb[nb[nb[x][2]][4]][0]];

    Uy = f.Uy[x] + f.Uy[nb[x][0]] + f.Uy[nb[x][4]] + f.Uy[nb[nb[x][0]][4]];
    Uyp = f.Uy[nb[x][2]] + f.Uy[nb[nb[x][0]][2]] + f.Uy[nb[nb[x][4]][2]]
      + f.Uy[nb[nb[nb[x][0]][4]][2]];

    Uz = f.Uz[x] + f.Uz[nb[x][0]] + f.Uz[nb[x][2]] + f.Uz[nb[nb[x][0]][2]];
    Uzp = f.Uz[nb[x][4]] + f.Uz[nb[nb[x][0]][4]] + f.Uz[nb[nb[x][2]][4]]
      + f.Uz[nb[nb[nb[x][0]][2]][4]];

    zoneDiv[x] =  0.25*(Uxp - Ux + Uyp - Uy + Uzp - Uz)/p.dx;

  }




  // x-direction art visc
  for(x=0; x<p.N; x++) {
    Umin = minof3(f.Ux[nb[x][1]],f.Ux[x],f.Ux[nb[x][0]]);
    Umax = maxof3(f.Ux[nb[x][1]],f.Ux[x],f.Ux[nb[x][0]]);
    Udix = minof2(Umax-f.Ux[x],f.Ux[x]-Umin);

    if((f.Ux[nb[x][0]]-f.Ux[x])*(f.Ux[x]-f.Ux[nb[x][1]]) > 0) {
      Udiv = Udix;
    } else {
      Udiv = 0;
    }

    c[x] = 0.5*(f.Ux[x] - f.Ux[nb[x][0]] + f.Ux[nb[x][1]] - f.Ux[x])/p.dx;
  }

  for(x=0; x<p.N; x++) {
    g = sqrt(1.0 + f.Ux[x]*f.Ux[x] + f.Uy[x]*f.Uy[x] + f.Uz[x]*f.Uz[x]);
    F = f.kappa[x]*f.E[x] + f.kappa[nb[x][5]]*f.E[nb[x][5]]
      + f.kappa[nb[x][3]]*f.E[nb[x][3]] + f.kappa[nb[nb[x][5]][3]]*f.E[nb[nb[x][5]][3]];
    delta = f.Ux[nb[x][0]] - f.Ux[x] - 0.5*p.dx*(c[x] + c[nb[x][0]])/g;



    if(delta*(f.Ux[nb[x][0]]-f.Ux[x]) > 0.0)
      deltam = delta;
    else
      deltam = 0.0;

    deltav = minof2(deltam,0.0);

    qx[x] = F*(Q*deltav - L*cs);
    tx[x] = F*deltav*(f.Vx[x] + f.Vx[nb[x][0]]);
  }


  for(x=0;x<p.N;x++) {
    if(zoneDiv[x] > 0)
      q[x] = 0.0;
    else
      q[x] = 0.25*(qx[x] + qx[nb[x][4]] + qx[nb[x][2]] + qx[nb[nb[x][2]][4]]);

    a = tx[x] + tx[nb[x][2]] + tx[nb[x][4]] + tx[nb[nb[x][2]][4]];
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
