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
#warning evolve_field p.C term not done!
    /*
    f.pi[x] = f.pi[x] 
      - 0.5*p.C*f.W[x]*(f.v[nb[x][1]]*(f.phi[x]-f.phi[nb[x][1]])
			 + f.v[x]*(f.phi[nb[x][0]] - f.phi[x]))/p.dx;
    */

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

  int x, y, z;
  
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



  // Section 3.4.4 -- pressure terms
  for(x = 0; x < p.N; x++) {
	qx = (f.Vx[nb[x][0]] - f.Vx[x])/p.dx;
	qy = (f.Vy[nb[x][2]] - f.Vy[x])/p.dx;
	qz = (f.Vz[nb[x][4]] - f.Vz[x])/p.dx;

	divv = qx + qy + qz;

	/*	
	s = 0.5*(f.kappa[x] - 1.0)*divv*p.dt/2.0;
	f.E[x] = f.E[x]*(1-s)/(1+s);
	*/

	f.E[x] = f.E[x]*exp(-1.0*(f.kappa[x] - 1.0)*divv*p.dt);

  }


  // update velocity v; denominator is W&M eq (2.85)
  // but note grid is Eulerian and D=0
  // then gv is the four-velocity (2.84)
  // (what we call kappa they call sigma, ish?)
  //
  // Section 3.4.5, equations 3.5.7, 3.5.8
  for(x = 0; x < p.N; x++) {

    sigmabar = (f.kappa[x]*f.E[x] 
		+ f.kappa[nb[x][1]]*f.E[nb[x][1]]
		+ f.kappa[nb[nb[x][1]][3]]*f.E[nb[nb[x][1]][3]]
		+ f.kappa[nb[nb[nb[x][1]][3]][5]]*f.E[nb[nb[nb[x][1]][3]][5]]
		+ f.kappa[nb[nb[x][3]][5]]*f.E[nb[nb[x][3]][5]]
		+ f.kappa[nb[nb[x][1]][5]]*f.E[nb[nb[x][1]][5]]
		+ f.kappa[nb[x][3]]*f.E[nb[x][3]]
		+ f.kappa[nb[x][5]]*f.E[nb[x][5]])/8.0;
    
    f.Ux[x] = f.Zx[x]/sigmabar;
    f.Uy[x] = f.Zy[x]/sigmabar;
    f.Uz[x] = f.Zz[x]/sigmabar;
  }


  // section 3.4.5 continued
  // face-centred x-cpt of 3-velocity
  for(x = 0; x < p.N; x++) {

    ubarx = (f.Ux[x] + f.Ux[nb[x][2]] 
	     + f.Ux[nb[x][4]] + f.Ux[nb[nb[x][2]][4]])/4.0;
    
    ubary = (f.Uy[x] + f.Uy[nb[x][2]]
	     + f.Uy[nb[x][4]] + f.Uy[nb[nb[x][2]][4]])/4.0;
    
    ubarz = (f.Uz[x] + f.Uz[nb[x][2]]
	     + f.Uz[nb[x][4]] + f.Uz[nb[nb[x][2]][4]])/4.0;
    
    Wfacex[x] = sqrt(1.0 + ubarx*ubarx + ubary*ubary + ubarz*ubarz);
    
    f.Vx[x] = ubarx/Wfacex[x];
  }
  // y-cpt
  for(x = 0; x < p.N; x++) {

	ubarx = (f.Ux[x] + f.Ux[nb[x][0]] 
		 + f.Ux[nb[x][4]] + f.Ux[nb[nb[x][0]][4]])/4.0;

	ubary = (f.Uy[x] + f.Uy[nb[x][0]] 
		 + f.Uy[nb[x][4]] + f.Uy[nb[nb[x][0]][4]])/4.0;

	ubarz = (f.Uz[x] + f.Uz[nb[x][0]] 
		 + f.Uz[nb[x][4]] + f.Uz[nb[nb[x][0]][4]])/4.0;

	Wfacey[x] = sqrt(1.0 + ubarx*ubarx + ubary*ubary + ubarz*ubarz);

	f.Vy[x] = ubary/Wfacey[x];

  }
  // z-cpt
  for(x = 0; x < p.N; x++) {

    ubarx = (f.Ux[x] + f.Ux[nb[x][2]] 
	     + f.Ux[nb[x][0]] + f.Ux[nb[nb[x][0]][2]])/4.0;
    
    ubary = (f.Uy[x] + f.Uy[nb[x][2]] 
	     + f.Uy[nb[x][0]] + f.Uy[nb[nb[x][0]][2]])/4.0;
    
    ubarz = (f.Uz[x] + f.Uz[nb[x][2]] 
	     + f.Uz[nb[x][0]] + f.Uz[nb[nb[x][0]][2]])/4.0;
    
    Wfacez[x] = sqrt(1.0 + ubarx*ubarx + ubary*ubary + ubarz*ubarz);
    
    f.Vz[x] = ubarz/Wfacez[x];
  }

  // End section 3.4.5




  // Section 3.4.6
  for(x = 0; x < p.L; x++) {
    for(y = 0; y < p.L; y++) {
      for(z = 0; z < p.L; z++) {

	utildex = (f.Ux[iix(x,y,z,p)] 
		   + f.Ux[iix(x+1,y,z,p)] + f.Ux[iix(x,y+1,z,p)] + f.Ux[iix(x,y,z+1,p)]
		   + f.Ux[iix(x+1,y,z+1,p)] + f.Ux[iix(x,y+1,z+1,p)] + f.Ux[iix(x+1,y+1,z,p)]
		   + f.Ux[iix(x+1,y+1,z+1,p)])/8.0;

	utildey = (f.Uy[iix(x,y,z,p)] 
		   + f.Uy[iix(x+1,y,z,p)] + f.Uy[iix(x,y+1,z,p)] + f.Uy[iix(x,y,z+1,p)]
		   + f.Uy[iix(x+1,y,z+1,p)] + f.Uy[iix(x,y+1,z+1,p)] + f.Uy[iix(x+1,y+1,z,p)]
		   + f.Uy[iix(x+1,y+1,z+1,p)])/8.0;

	utildez = (f.Uz[iix(x,y,z,p)] 
		   + f.Uz[iix(x+1,y,z,p)] + f.Uz[iix(x,y+1,z,p)] + f.Uz[iix(x,y,z+1,p)]
		   + f.Uz[iix(x+1,y,z+1,p)] + f.Uz[iix(x,y+1,z+1,p)] + f.Uz[iix(x+1,y+1,z,p)]
		   + f.Uz[iix(x+1,y+1,z+1,p)])/8.0;

	Wold[iix(x,y,z,p)] = f.W[iix(x,y,z,p)];
	f.W[iix(x,y,z,p)] = sqrt(1.0 + utildex*utildex + utildey*utildey + utildez*utildez);

	//	s = (f.kappa[iix(x,y,z,p)] - 1)*(f.W[iix(x,y,z,p)] - Wold[iix(x,y,z,p)])
	//	  /(f.W[iix(x,y,z,p)] + Wold[iix(x,y,z,p)]);
	//	f.E[iix(x,y,z,p)] = f.E[iix(x,y,z,p)]*(1-s)/(1+s);
	// sort of E*exp(-2.0*s)


	// Top of p90 ch 3
	f.E[iix(x,y,z,p)] = f.E[iix(x,y,z,p)]*pow(Wold[iix(x,y,z,p)]/f.W[iix(x,y,z,p)],f.kappa[iix(x,y,z,p)]-1.0);






      }
    }
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




  // Pressure acceleration (and artificial viscosity for 'Z')
  // W&M sec 2.4.12, 3.5.1, DONE
  for(x = 0; x < p.L; x++) {
    for(y = 0; y < p.L; y++) {
      for(z = 0; z < p.L; z++) {

	// equation (3.68)
	// repeated calculation of these unnecessary, so split it out later
	p_bar_x_plus = (f.p[iix(x,y,z,p)] + f.p[iix(x,y-1,z,p)] 
			+ f.p[iix(x,y,z-1,p)] + f.p[iix(x,y-1,z-1,p)])/4.0;

	p_bar_x_minus = (f.p[iix(x-1,y,z,p)] + f.p[iix(x-1,y-1,z,p)] 
			 + f.p[iix(x-1,y,z-1,p)] + f.p[iix(x-1,y-1,z-1,p)])/4.0;

	p_bar_y_plus = (f.p[iix(x,y,z,p)] + f.p[iix(x-1,y,z,p)] 
			+ f.p[iix(x,y,z-1,p)] + f.p[iix(x-1,y,z-1,p)])/4.0;


	p_bar_y_minus = (f.p[iix(x,y-1,z,p)] + f.p[iix(x-1,y-1,z,p)] 
			 + f.p[iix(x,y-1,z-1,p)] + f.p[iix(x-1,y-1,z-1,p)])/4.0;


	p_bar_z_plus = (f.p[iix(x,y,z,p)] + f.p[iix(x-1,y,z,p)] 
			+ f.p[iix(x,y-1,z,p)] + f.p[iix(x-1,y-1,z,p)])/4.0;


	p_bar_z_minus = (f.p[iix(x,y,z-1,p)] + f.p[iix(x-1,y,z-1,p)] 
			 + f.p[iix(x,y-1,z-1,p)] + f.p[iix(x-1,y-1,z-1,p)])/4.0;


	

	f.Zx[iix(x,y,z,p)] = f.Zx[iix(x,y,z,p)] - (p_bar_x_plus - p_bar_x_minus)*p.dt/p.dx;

	f.Zy[iix(x,y,z,p)] = f.Zy[iix(x,y,z,p)] - (p_bar_y_plus - p_bar_y_minus)*p.dt/p.dx;

	f.Zz[iix(x,y,z,p)] = f.Zz[iix(x,y,z,p)] - (p_bar_z_plus - p_bar_z_minus)*p.dt/p.dx;
      }
    }
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
