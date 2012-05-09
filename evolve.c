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



  // Section 3.4.4

  for(x = 0; x < p.L; x++) {
    for(y = 0; y < p.L; y++) {
      for(z = 0; z < p.L; z++) {
	qx = (f.Vx[iix(x+1,y,z,p)] - f.Vx[iix(x,y,z,p)])/p.dx;
	qy = (f.Vy[iix(x,y+1,z,p)] - f.Vy[iix(x,y,z,p)])/p.dx;
	qz = (f.Vz[iix(x,y,z+1,p)] - f.Vz[iix(x,y,z,p)])/p.dx;

	divv = 0.5*(qx+qy+qz);
	

	s = (f.kappa[x] - 1.0)*divv*p.dt/2.0;

	f.E[iix(x,y,z,p)] = f.E[iix(x,y,z,p)]*(1-s)/(1+s);



      }
    }
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


	
	//	if(f.E[iix(x,y,z,p)] > 10)


	f.Zx[iix(x,y,z,p)] = (p_bar_x_plus - p_bar_x_minus)*p.dt/p.dx;

	f.Zy[iix(x,y,z,p)] = (p_bar_y_plus - p_bar_y_minus)*p.dt/p.dx;

	f.Zz[iix(x,y,z,p)] = (p_bar_z_plus - p_bar_z_minus)*p.dt/p.dx;
      }
    }
  }





  // update velocity v; denominator is W&M eq (2.85)
  // but note grid is Eulerian and D=0
  // then gv is the four-velocity (2.84)
  // (what we call kappa they call sigma, ish?)
  //
  // Section 3.4.5, equations 3.5.7, 3.5.8
  for(x = 0; x < p.L; x++) {
    for(y = 0; y < p.L; y++) {
      for(z = 0; z < p.L; z++) {  

	sigmabar = (f.kappa[iix(x,y,z,p)]*f.E[iix(x,y,z,p)] 
		    + f.kappa[iix(x-1,y,z,p)]*f.E[iix(x-1,y,z,p)]
		    + f.kappa[iix(x-1,y-1,z,p)]*f.E[iix(x-1,y-1,z,p)]
		    + f.kappa[iix(x-1,y-1,z-1,p)]*f.E[iix(x-1,y-1,z-1,p)]
		    + f.kappa[iix(x,y-1,z-1,p)]*f.E[iix(x,y-1,z-1,p)]
		    + f.kappa[iix(x-1,y,z-1,p)]*f.E[iix(x-1,y,z-1,p)]
		    + f.kappa[iix(x,y-1,z,p)]*f.E[iix(x,y-1,z,p)]
		    + f.kappa[iix(x,y,z-1,p)]*f.E[iix(x,y,z-1,p)])/8.0;

	f.Ux[iix(x,y,z,p)] = f.Zx[iix(x,y,z,p)]/sigmabar;
	f.Uy[iix(x,y,z,p)] = f.Zy[iix(x,y,z,p)]/sigmabar;
	f.Uz[iix(x,y,z,p)] = f.Zz[iix(x,y,z,p)]/sigmabar;
      }
    }
  }




  // face-centred x-cpt of 3-velocity
  for(x = 0; x < p.L; x++) {
    for(y = 0; y < p.L; y++) {
      for(z = 0; z < p.L; z++) {  

	ubarx = (f.Ux[iix(x,y,z,p)] + f.Ux[iix(x,y+1,z,p)] 
		 + f.Ux[iix(x,y,z+1,p)] + f.Ux[iix(x,y+1,z+1,p)])/4.0;

	ubary = (f.Uy[iix(x,y,z,p)] + f.Uy[iix(x,y+1,z,p)] 
		 + f.Uy[iix(x,y,z+1,p)] + f.Uy[iix(x,y+1,z+1,p)])/4.0;

	ubarz = (f.Uz[iix(x,y,z,p)] + f.Uz[iix(x,y+1,z,p)] 
		 + f.Uz[iix(x,y,z+1,p)] + f.Uz[iix(x,y+1,z+1,p)])/4.0;

	Wfacex[iix(x,y,z,p)] = sqrt(1.0 + ubarx*ubarx + ubary*ubary + ubarz*ubarz);

	f.Vx[iix(x,y,z,p)] = ubarx/Wfacex[iix(x,y,z,p)];
      }
    }
  }
  // y-cpt
  for(x = 0; x < p.L; x++) {
    for(y = 0; y < p.L; y++) {
      for(z = 0; z < p.L; z++) {  

	ubarx = (f.Ux[iix(x,y,z,p)] + f.Ux[iix(x+1,y,z,p)] 
		 + f.Ux[iix(x,y,z+1,p)] + f.Ux[iix(x+1,y,z+1,p)])/4.0;

	ubary = (f.Uy[iix(x,y,z,p)] + f.Uy[iix(x+1,y,z,p)] 
		 + f.Uy[iix(x,y,z+1,p)] + f.Uy[iix(x+1,y,z+1,p)])/4.0;

	ubarz = (f.Uz[iix(x,y,z,p)] + f.Uz[iix(x+1,y,z,p)] 
		 + f.Uz[iix(x,y,z+1,p)] + f.Uz[iix(x+1,y,z+1,p)])/4.0;

	Wfacey[iix(x,y,z,p)] = sqrt(1.0 + ubarx*ubarx + ubary*ubary + ubarz*ubarz);

	f.Vy[iix(x,y,z,p)] = ubary/Wfacey[iix(x,y,z,p)];
      }
    }
  }
  // z-cpt
  for(x = 0; x < p.L; x++) {
    for(y = 0; y < p.L; y++) {
      for(z = 0; z < p.L; z++) {  

	ubarx = (f.Ux[iix(x,y,z,p)] + f.Ux[iix(x,y+1,z,p)] 
		 + f.Ux[iix(x+1,y,z,p)] + f.Ux[iix(x+1,y+1,z,p)])/4.0;

	ubary = (f.Uy[iix(x,y,z,p)] + f.Uy[iix(x,y+1,z,p)] 
		 + f.Uy[iix(x+1,y,z,p)] + f.Uy[iix(x+1,y+1,z,p)])/4.0;

	ubarz = (f.Uz[iix(x,y,z,p)] + f.Uz[iix(x,y+1,z,p)] 
		 + f.Uz[iix(x+1,y,z,p)] + f.Uz[iix(x+1,y+1,z,p)])/4.0;

	Wfacez[iix(x,y,z,p)] = sqrt(1.0 + ubarx*ubarx + ubary*ubary + ubarz*ubarz);

	f.Vz[iix(x,y,z,p)] = ubarz/Wfacez[iix(x,y,z,p)];
      }
    }
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

	s = (f.kappa[iix(x,y,z,p)] - 1)*(f.W[iix(x,y,z,p)] - Wold[iix(x,y,z,p)])
	  /(f.W[iix(x,y,z,p)] + Wold[iix(x,y,z,p)]);
	// sort of E*exp(-2.0*s)



	f.E[iix(x,y,z,p)] = f.E[iix(x,y,z,p)]*(1-s)/(1+s);



      }
    }
  }




  for(x = 0; x < p.L; x++) {
    for(y = 0; y < p.L; y++) {
      for(z = 0; z < p.L; z++) {
	qx = (f.Vx[iix(x,y,z,p)] + f.Vx[iix(x+1,y,z,p)])
	  *(Wfacex[iix(x+1,y,z,p)] - Wfacex[iix(x,y,z,p)])*p.dt/(2.0*p.dx);

	qy = (f.Vy[iix(x,y,z,p)] + f.Vy[iix(x,y+1,z,p)])
	  *(Wfacey[iix(x,y+1,z,p)] - Wfacey[iix(x,y,z,p)])*p.dt/(2.0*p.dx);

	qz = (f.Vz[iix(x,y,z,p)] + f.Vz[iix(x,y,z+1,p)])
	  *(Wfacez[iix(x,y,z+1,p)] - Wfacez[iix(x,y,z,p)])*p.dt/(2.0*p.dx);

	gradv = (qx+qy+qz)/f.W[iix(x,y,z,p)];

	s = (f.kappa[x] - 1.0)*gradv/2.0;


	f.E[iix(x,y,z,p)] = f.E[iix(x,y,z,p)]*(1-s)/(1+s);



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
