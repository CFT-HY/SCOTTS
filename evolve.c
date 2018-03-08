/* evolve.c
 *
 * Everything to do with evolving fields and fluids with the
 * equations of motion.
 */
#include "hydro.h"


/* void evolve_backstep(hydro_fields f, hydro_params p)
 *
 * Evolve the scalar field momentum back a halfstep.
 * Taken straight from the 1+1 spherical fortran code.
 */
void evolve_backstep(hydro_fields f, hydro_params p) {

  int x, y, z;

  for(x = 1; x <= p.slicex; x++) {
    for(y = 1; y <= p.slicey; y++) {
      for(z = 0; z < p.Lz; z++) {
	f.pi[x][y][z] = f.pifull[x][y][z] - 0.5*p.dt
	  *(f.phi[x+1][y][z] + f.phi[x][y+1][z] 
	    + f.phi[x][y][((z+1)%p.Lz)] 
	    - 6.0*f.phi[x][y][z] 
	    + f.phi[x-1][y][z] + f.phi[x][y-1][z]
	    + f.phi[x][y][((z-1+p.Lz)%p.Lz)])/(p.dx*p.dx)
	  + 0.125*p.dt*p.dt
	  *(f.pifull[x+1][y][z] - 2.0*f.pifull[x][y][z]
	    + f.pifull[x-1][y][z])/(p.dx*p.dx)
#ifndef SCALAR
	  + 0.5*p.dt*Vdf(p, f.T[x][y][z], f.phi[x][y][z] 
			 - 0.25*p.dt*f.pifull[x][y][z]);
#else
	  + 0.5*p.dt*Vdf(p, p.Tconst, f.phi[x][y][z] 
			 - 0.25*p.dt*f.pifull[x][y][z]);
#endif
      }
    }
  }

  halo_field(f.pi, p);
}


/* void evolve_field(hydro_fields f, hydro_params p)
 *
 * Evolve the scalar field forward one timestep.
 * Taken more or less straight from the 1+1 spherical fortran code.
 */
void evolve_field(hydro_fields f, hydro_params p) {
  
  int x, y, z;

  float piold, s;

  // Move conjugate momentum (leapfrog)
  for(x = 1; x <= p.slicex; x++) {
    for(y = 1; y <= p.slicey; y++) {
      for(z = 0; z < p.Lz; z++) {
	
	
	piold = f.pi[x][y][z];

#ifndef SCALAR
	// first-order viscosity term
	s = -0.5*p.dt*(p.C*f.phi[x][y][z]*f.phi[x][y][z]/f.T[x][y][z])*f.W[x][y][z];

	f.pi[x][y][z] = (1+s)*f.pi[x][y][z]/(1-s);
	
	// gradient term
	f.pi[x][y][z] = f.pi[x][y][z] +s*(
					
					  (f.V[0][x][y][z]*(f.phi[x][y][z]
							    - f.phi[x-1][y][z])
					   + f.V[0][x+1][y][z]*(f.phi[x+1][y][z]
								- f.phi[x][y][z]))
					 
					  +(f.V[1][x][y][z]*(f.phi[x][y][z]
							     - f.phi[x][y-1][z])
					    + f.V[1][x][y+1][z]*(f.phi[x][y+1][z]
								 - f.phi[x][y][z]))
					  
					  +(f.V[2][x][y][z]*(f.phi[x][y][z]
							     - f.phi[x][y][((z-1+p.Lz)%p.Lz)])
					    + f.V[2][x][y][((z+1)%p.Lz)]*(f.phi[x][y][((z+1)%p.Lz)]
									  - f.phi[x][y][z]))	
					  )/(p.dx*(1-s));
#endif // !SCALAR
	
	// scalar field gradient and potential terms
	f.pi[x][y][z] = f.pi[x][y][z]
	  + p.dt
	  *((f.phi[x+1][y][z] + f.phi[x][y+1][z]
	     + f.phi[x][y][((z+1)%p.Lz)]
	     - 6.0*f.phi[x][y][z] 
	     + f.phi[x-1][y][z] + f.phi[x][y-1][z]
	     + f.phi[x][y][((z-1+p.Lz)%p.Lz)]
	     )/(p.dx*p.dx)
#ifndef SCALAR
	    - Vdf(p, f.T[x][y][z], f.phi[x][y][z]))/(1-s);
#else	
	    - Vdf(p, p.Tconst, f.phi[x][y][z]));
#endif // !SCALAR	

#ifndef SCALAR
      // pifull is (1.5*pi - 0.5*piold), not sure why
      f.pifull[x][y][z] = piold + 1.5*(f.pi[x][y][z] - piold);
#else
      // more reasonable?
      f.pifull[x][y][z] = 0.5*(piold + f.pi[x][y][z]);
#endif // !SCALAR	

	
      }
    }
  }

  
  // Move field (leapfrog)
  for(x = 1; x <= p.slicex; x++) {
    for(y = 1; y <= p.slicey; y++) {
      for(z = 0; z < p.Lz; z++) {
	f.phiold[x][y][z] = f.phi[x][y][z];
	f.phi[x][y][z] = f.phi[x][y][z] + p.dt*f.pi[x][y][z];
      }
    }
  }
  
  halo_field(f.pi, p);
  halo_field(f.pifull, p);
  
  halo_field(f.phi, p);
  halo_field(f.phiold, p);
  
}




/* void evolve_hydro(hydro_fields f, hydro_params p)
 *
 * Basically straight from fortran. Original comment:
 *
 *   --  Does hydro except transport.  E  and  Z evolved.  v  and W
 *       obtained.  Uses artificial viscosity  Q,  adjusted with  Cav.
 *   --  The order chosen is not unique, but is what I used in
 *       dissertation.  Should experiment with different orders,
 *       especially with position of field-fluid interaction,
 *       now placed first (otherwise order is from CW).
 *
 */
void evolve_hydro(hydro_fields f, hydro_params p) {

#ifndef SCALAR

  int x, y, z;
  
  float ***Vdmid = make_field(p);
  float ***Wold = make_field(p);
  float ***phiav = make_field(p);

  float ****dxphi = make_vector(p);

  float ***Wfacex =  make_field(p);
  float ***Wfacey =  make_field(p);
  float ***Wfacez =  make_field(p);

  
  float gpi, dv, gv, s;

  float p_bar_x_plus, p_bar_x_minus;
  float p_bar_y_plus, p_bar_y_minus;
  float p_bar_z_plus, p_bar_z_minus;

  float utildex, utildey, utildez;
  float qx, qy, qz, gradv;
  float divv;

  float sigmabar;

  float ubarx, ubary, ubarz, W;

  float vdnb, pinb, phinb, wnb, Tnb, dxphinb0, dxphinb1, dxphinb2;



  for(x = 1; x <= p.slicex; x++) {
    for(y = 1; y <= p.slicey; y++) {
      for(z = 0; z < p.Lz; z++) {


	
	phiav[x][y][z] = 0.5*(f.phiold[x][y][z] + f.phi[x][y][z]);
	
	dxphi[0][x][y][z] =  0.5*(f.phiold[x+1][y][z]
				  + f.phi[x+1][y][z]
				  - f.phiold[x][y][z] - f.phi[x][y][z])/p.dx;
	
	dxphi[1][x][y][z] =  0.5*(f.phiold[x][y+1][z]
				  + f.phi[x][y+1][z]
				  - f.phiold[x][y][z] - f.phi[x][y][z])/p.dx;
	
	dxphi[2][x][y][z] =  0.5*(f.phiold[x][y][((z+1)%p.Lz)]
				  + f.phi[x][y][((z+1)%p.Lz)]
				  - f.phiold[x][y][z] - f.phi[x][y][z])/p.dx;
	
      }
    }
  }

  halo_field(dxphi[0], p);
  halo_field(dxphi[1], p);
  halo_field(dxphi[2], p);

  Vdpot(p, f.T[0][0], phiav[0][0], Vdmid[0][0]);


  halo_field(Vdmid, p);



  for(x = 1; x <= p.slicex; x++) {
    for(y = 1; y <= p.slicey; y++) {
      for(z = 0; z < p.Lz; z++) {

	vdnb = 0.125*(Vdmid[x-1][y][z]
		      + Vdmid[x][y][z]
		      + Vdmid[x][y-1][z]
		      + Vdmid[x-1][y-1][z]
		      + Vdmid[x][y][((z-1+p.Lz)%p.Lz)]
		      + Vdmid[x][y-1][((z+p.Lz-1)%p.Lz)]
		      + Vdmid[x-1][y][((z+p.Lz-1)%p.Lz)]
		      + Vdmid[x-1][y-1][((z+p.Lz-1)%p.Lz)]);
	
	
	pinb = 0.125*(f.pi[x-1][y][z]
		      + f.pi[x][y][z]
		      + f.pi[x][y-1][z]
		      + f.pi[x-1][y-1][z]
		      + f.pi[x][y][((z-1+p.Lz)%p.Lz)]
		      + f.pi[x][y-1][((z+p.Lz-1)%p.Lz)]
		      + f.pi[x-1][y][((z+p.Lz-1)%p.Lz)]
		      + f.pi[x-1][y-1][((z+p.Lz-1)%p.Lz)]);


	Tnb = 0.125*(f.T[x-1][y][z]
		      + f.T[x][y][z]
		      + f.T[x][y-1][z]
		      + f.T[x-1][y-1][z]
		      + f.T[x][y][((z-1+p.Lz)%p.Lz)]
		      + f.T[x][y-1][((z+p.Lz-1)%p.Lz)]
		      + f.T[x-1][y][((z+p.Lz-1)%p.Lz)]
		      + f.T[x-1][y-1][((z+p.Lz-1)%p.Lz)]);

	phinb = 0.125*(f.phi[x-1][y][z]
		      + f.phi[x][y][z]
		      + f.phi[x][y-1][z]
		      + f.phi[x-1][y-1][z]
		      + f.phi[x][y][((z-1+p.Lz)%p.Lz)]
		      + f.phi[x][y-1][((z+p.Lz-1)%p.Lz)]
		      + f.phi[x-1][y][((z+p.Lz-1)%p.Lz)]
		      + f.phi[x-1][y-1][((z+p.Lz-1)%p.Lz)]);

	
	
	wnb = 0.125*(f.W[x-1][y][z]
		     + f.W[x][y][z]
		     + f.W[x][y-1][z]
		     + f.W[x-1][y-1][z]
		     + f.W[x][y][((z-1+p.Lz)%p.Lz)]
		     + f.W[x][y-1][((z+p.Lz-1)%p.Lz)]

		     + f.W[x-1][y][((z+p.Lz-1)%p.Lz)]
		     + f.W[x-1][y-1][((z+p.Lz-1)%p.Lz)]);
	
	
	dxphinb0 = 0.25*(dxphi[0][x-1][y][z]
			 + dxphi[0][x-1][y-1][z]
			 + dxphi[0][x-1][y][((z+p.Lz-1)%p.Lz)]
			 + dxphi[0][x-1][y-1][((z+p.Lz-1)%p.Lz)]);
	
	dxphinb1 = 0.25*(dxphi[1][x][y-1][z]
			 + dxphi[1][x-1][y-1][z]
			 + dxphi[1][x][y-1][((z+p.Lz-1)%p.Lz)]
			 + dxphi[1][x-1][y-1][((z+p.Lz-1)%p.Lz)]);
	
	dxphinb2 = 0.25*(dxphi[2][x][y][((z-1+p.Lz)%p.Lz)]
			 + dxphi[2][x-1][y][((z+p.Lz-1)%p.Lz)]
			 + dxphi[2][x][y-1][((z+p.Lz-1)%p.Lz)]
			 + dxphi[2][x-1][y-1][((z+p.Lz-1)%p.Lz)]);


	
	// this is the one that causes trouble       
	f.Z[0][x][y][z] = f.Z[0][x][y][z] - p.dt*(p.C*phinb*phinb/Tnb)*wnb*pinb*dxphinb0;
	f.Z[1][x][y][z] = f.Z[1][x][y][z] - p.dt*(p.C*phinb*phinb/Tnb)*wnb*pinb*dxphinb1;
	f.Z[2][x][y][z] = f.Z[2][x][y][z] - p.dt*(p.C*phinb*phinb/Tnb)*wnb*pinb*dxphinb2;
    
	f.Z[1][x][y][z] = f.Z[1][x][y][z] 
	  - p.dt*(p.C*phinb*phinb/Tnb)*wnb*(f.V[1][x][y][z]*dxphinb1
					    + f.V[0][x][y][z]*dxphinb0
					    + f.V[2][x][y][z]*dxphinb2)*dxphinb1;
	
	f.Z[0][x][y][z] = f.Z[0][x][y][z] 
	  - p.dt*(p.C*phinb*phinb/Tnb)*wnb*(f.V[1][x][y][z]*dxphinb1
					    + f.V[0][x][y][z]*dxphinb0
					    + f.V[2][x][y][z]*dxphinb2)*dxphinb0;    
    
	f.Z[2][x][y][z] = f.Z[2][x][y][z]
	  - p.dt*(p.C*phinb*phinb/Tnb)*wnb*(f.V[1][x][y][z]*dxphinb1
					    + f.V[0][x][y][z]*dxphinb0
					    + f.V[2][x][y][z]*dxphinb2)*dxphinb2;
	
	f.Z[0][x][y][z] = f.Z[0][x][y][z] - p.dt*vdnb*dxphinb0;     
	f.Z[1][x][y][z] = f.Z[1][x][y][z] - p.dt*vdnb*dxphinb1;
	f.Z[2][x][y][z] = f.Z[2][x][y][z] - p.dt*vdnb*dxphinb2;
				       
	gpi = f.W[x][y][z]*(f.pi[x][y][z]
			    + 0.5*(f.V[0][x][y][z]
				   *dxphi[0][x-1][y][z]
				   + f.V[0][x+1][y][z]
				   *dxphi[0][x][y][z])
			    + 0.5*(f.V[1][x][y][z]
				   *dxphi[1][x][y-1][z]
				   + f.V[1][x][y+1][z]
				   *dxphi[1][x][y][z])
			    + 0.5*(f.V[2][x][y][z]
				   *dxphi[2][x][y][((z-1+p.Lz)%p.Lz)]
				   + f.V[2][x][y][((z+1)%p.Lz)]
				   *dxphi[2][x][y][z]));
      
	f.E[x][y][z] = f.E[x][y][z] + p.dt*((p.C*phiav[x][y][z]*phiav[x][y][z]/f.T[x][y][z])*gpi*gpi + gpi*Vdmid[x][y][z]);

    
      }
    }
  }


  // halo_field(f.Z[0], p);
  // halo_field(f.Z[1], p);
  // halo_field(f.Z[2], p);
  //  halo_field(f.E, p);



  // Pressure acceleration
  // W&M sec 2.4.12, 3.5.1, DONE
  for(x = 1; x <= p.slicex; x++) {
    for(y = 1; y <= p.slicey; y++) {
      for(z = 0; z < p.Lz; z++) {

	// equation (3.68)
	// repeated calculation of these unnecessary, so split it out later
	p_bar_x_plus = (f.p[x][y][z]
			+ f.p[x][y-1][z] 
			+ f.p[x][y][((z-1+p.Lz)%p.Lz)]
			+ f.p[x][y-1][((z+p.Lz-1)%p.Lz)]
			)/4.0;
	
	p_bar_x_minus = (f.p[x-1][y][z] 
			 + f.p[x-1][y-1][z]
			 + f.p[x-1][y][((z+p.Lz-1)%p.Lz)]
			 + f.p[x-1][y-1][((z+p.Lz-1)%p.Lz)]
			 )/4.0;
	
	p_bar_y_plus = (f.p[x][y][z]
			+ f.p[x-1][y][z] 
			+ f.p[x][y][((z-1+p.Lz)%p.Lz)]
			+ f.p[x-1][y][((z+p.Lz-1)%p.Lz)]
			)/4.0;
	
	
	p_bar_y_minus = (f.p[x][y-1][z]
			 + f.p[x-1][y-1][z] 
			 + f.p[x][y-1][((z+p.Lz-1)%p.Lz)]
			 + f.p[x-1][y-1][((z+p.Lz-1)%p.Lz)]
			 )/4.0;
	
	
	p_bar_z_plus = (f.p[x][y][z]
			+ f.p[x-1][y][z] 
			+ f.p[x][y-1][z]
			+ f.p[x-1][y-1][z]
		    )/4.0;
	
	
	p_bar_z_minus = (f.p[x][y][((z-1+p.Lz)%p.Lz)]
			 + f.p[x-1][y][((z+p.Lz-1)%p.Lz)] 
			 + f.p[x][y-1][((z+p.Lz-1)%p.Lz)]
			 + f.p[x-1][y-1][((z+p.Lz-1)%p.Lz)]
			 )/4.0;
	
	
	
	
	f.Z[0][x][y][z] = f.Z[0][x][y][z] - (p_bar_x_plus
					     - p_bar_x_minus)*p.dt/p.dx;
	
	f.Z[1][x][y][z] = f.Z[1][x][y][z] - (p_bar_y_plus
					     - p_bar_y_minus)*p.dt/p.dx;
	
	f.Z[2][x][y][z] = f.Z[2][x][y][z] - (p_bar_z_plus
					     - p_bar_z_minus)*p.dt/p.dx;

      }
    }
  }

  //  halo_field(f.Z[0], p);
  //  halo_field(f.Z[1], p);
  //  halo_field(f.Z[2], p);



  // update velocity v; denominator is W&M eq (2.85)
  // but note grid is Eulerian and D=0
  // then gv is the four-velocity (2.84)
  // (what we call kappa they call sigma, ish?)
  //
  // Section 3.4.5, equations 3.57, 3.58

  for(x = 1; x <= p.slicex; x++) {
    for(y = 1; y <= p.slicey; y++) {
      for(z = 0; z < p.Lz; z++) {

    /*
    if(x == 0)
      fprintf(stderr,"old Ux = %lf\n", f.U[0][1]);
    */

	sigmabar = (f.kappa[x][y][z]
		    *f.E[x][y][z] 
		    + f.kappa[x-1][y][z]
		    *f.E[x-1][y][z]
		    + f.kappa[x-1][y-1][z]
		    *f.E[x-1][y-1][z]
		    + f.kappa[x-1][y-1][((z+p.Lz-1)%p.Lz)]
		    *f.E[x-1][y-1][((z+p.Lz-1)%p.Lz)]
		    + f.kappa[x][y-1][((z+p.Lz-1)%p.Lz)]
		    *f.E[x][y-1][((z+p.Lz-1)%p.Lz)]
		    + f.kappa[x-1][y][((z+p.Lz-1)%p.Lz)]
		    *f.E[x-1][y][((z+p.Lz-1)%p.Lz)]
		    + f.kappa[x][y-1][z]
		    *f.E[x][y-1][z]
		    + f.kappa[x][y][((z-1+p.Lz)%p.Lz)]
		    *f.E[x][y][((z-1+p.Lz)%p.Lz)]
		)/8.0;

   
    f.U[0][x][y][z] = f.Z[0][x][y][z]/sigmabar;
    f.U[1][x][y][z] = f.Z[1][x][y][z]/sigmabar;
    f.U[2][x][y][z] = f.Z[2][x][y][z]/sigmabar;
    
    //    fprintf(stderr,"%lf %lf\n", f.Z[0][x], sigmabar);
    /*
      if(x == 0)
      fprintf(stderr,"Ux = %lf\n", f.U[0][1]);
    */
      }
    }
  }

  halo_field(f.U[0], p);
  halo_field(f.U[1], p);
  halo_field(f.U[2], p);


  // section 3.4.5 continued
  // face-centred x-cpt of 3-velocity
  for(x = 1; x <= p.slicex; x++) {
    for(y = 1; y <= p.slicey; y++) {
      for(z = 0; z < p.Lz; z++) {

	ubarx = (f.U[0][x][y][z]
		 + f.U[0][x][y+1][z] 
		 + f.U[0][x][y][((z+1)%p.Lz)]
		 + f.U[0][x][y+1][((z+1)%p.Lz)]
		 )/4.0;
	
	ubary = (f.U[1][x][y][z]
		 + f.U[1][x][y+1][z]
		 + f.U[1][x][y][((z+1)%p.Lz)]
		 + f.U[1][x][y+1][((z+1)%p.Lz)]
		 )/4.0;
	
	ubarz = (f.U[2][x][y][z]
		 + f.U[2][x][y+1][z]
		 + f.U[2][x][y][((z+1)%p.Lz)]
		 + f.U[2][x][y+1][((z+1)%p.Lz)]
		 )/4.0;
	
	Wfacex[x][y][z] = sqrt(1.0 + ubarx*ubarx + ubary*ubary + ubarz*ubarz);
	
	f.V[0][x][y][z] = ubarx/Wfacex[x][y][z];
      }
    }
  }
  



  // y-cpt
  for(x = 1; x <= p.slicex; x++) {
    for(y = 1; y <= p.slicey; y++) {
      for(z = 0; z < p.Lz; z++) {
	
	ubarx = (f.U[0][x][y][z]
		 + f.U[0][x+1][y][z] 
		 + f.U[0][x][y][((z+1)%p.Lz)]
		 + f.U[0][x+1][y][((z+1)%p.Lz)]
		 )/4.0;

	ubary = (f.U[1][x][y][z]
		 + f.U[1][x+1][y][z] 
		 + f.U[1][x][y][((z+1)%p.Lz)]
		 + f.U[1][x+1][y][((z+1)%p.Lz)]
		 )/4.0;

	ubarz = (f.U[2][x][y][z]
		 + f.U[2][x+1][y][z] 
		 + f.U[2][x][y][((z+1)%p.Lz)]
		 + f.U[2][x+1][y][((z+1)%p.Lz)]
		 )/4.0;

	Wfacey[x][y][z] = sqrt(1.0 + ubarx*ubarx + ubary*ubary + ubarz*ubarz);

	f.V[1][x][y][z] = ubary/Wfacey[x][y][z];
	
      }
    }
  }

  // z-cpt
  for(x = 1; x <= p.slicex; x++) {
    for(y = 1; y <= p.slicey; y++) {
      for(z = 0; z < p.Lz; z++) {

	ubarx = (f.U[0][x][y][z]
		 + f.U[0][x][y+1][z] 
		 + f.U[0][x+1][y][z]
		 + f.U[0][x+1][y+1][z]
		 )/4.0;
	
	ubary = (f.U[1][x][y][z]
		 + f.U[1][x][y+1][z] 
		 + f.U[1][x+1][y][z]
		 + f.U[1][x+1][y+1][z]
		 )/4.0;
	
	ubarz = (f.U[2][x][y][z]
		 + f.U[2][x][y+1][z] 
		 + f.U[2][x+1][y][z]
		 + f.U[2][x+1][y+1][z]
		 )/4.0;
	
	Wfacez[x][y][z] = sqrt(1.0 
			       + ubarx*ubarx
			       + ubary*ubary
			       + ubarz*ubarz);
    
	f.V[2][x][y][z] = ubarz/Wfacez[x][y][z];
      }
    }
  }
  

  halo_field(f.V[0], p);
  halo_field(f.V[1], p);
  halo_field(f.V[2], p);



  // End section 3.4.5

  //  fprintf(stderr,"Vx[1] = %lf\n", f.V[0][1]);
  //  getchar();

  // Section 3.4.6
  for(x = 1; x <= p.slicex; x++) {
    for(y = 1; y <= p.slicey; y++) {
      for(z = 0; z < p.Lz; z++) {
	
	
	utildex = (f.U[0][x][y][z] 
		   + f.U[0][x+1][y][z]
		   + f.U[0][x][y+1][z]
		   + f.U[0][x][y][((z+1)%p.Lz)]
		   + f.U[0][x+1][y][((z+1)%p.Lz)]
		   + f.U[0][x][y+1][((z+1)%p.Lz)]
		   + f.U[0][x+1][y+1][z]
		   + f.U[0][x+1][y+1][((z+1)%p.Lz)]
		   )/8.0;
	
	utildey = (f.U[1][x][y][z] 
		   + f.U[1][x+1][y][z]
		   + f.U[1][x][y+1][z]
		   + f.U[1][x][y][((z+1)%p.Lz)]
		   + f.U[1][x+1][y][((z+1)%p.Lz)]
		   + f.U[1][x][y+1][((z+1)%p.Lz)]
		   + f.U[1][x+1][y+1][z]
		   + f.U[1][x+1][y+1][((z+1)%p.Lz)]
		   )/8.0;
    
	utildez = (f.U[2][x][y][z] 
		   + f.U[2][x+1][y][z]
		   + f.U[2][x][y+1][z]
		   + f.U[2][x][y][((z+1)%p.Lz)]
		   + f.U[2][x+1][y][((z+1)%p.Lz)]
		   + f.U[2][x][y+1][((z+1)%p.Lz)]
		   + f.U[2][x+1][y+1][z]
		   + f.U[2][x+1][y+1][((z+1)%p.Lz)]
		   )/8.0;

	
	//    if(x < 3 || x > p.N-3)
	//      fprintf(stderr,"utildex %d = %lf\n", x, utildex);
     
	Wold[x][y][z] = f.W[x][y][z];
	f.W[x][y][z] = sqrt(1.0 
			    + utildex*utildex 
			    + utildey*utildey 
			    + utildez*utildez);
	
	// sort of E*exp(-2.0*s)	
	//	s = (f.kappa[iix(x,y,z,p)] - 1)*(f.W[iix(x,y,z,p)]
	//					 - Wold[iix(x,y,z,p)])
	//	  /(f.W[iix(x,y,z,p)] + Wold[iix(x,y,z,p)]);
	//	f.E[iix(x,y,z,p)] = f.E[iix(x,y,z,p)]*(1-s)/(1+s);
    
    
	// Top of p90 ch 3
	f.E[x][y][z] = f.E[x][y][z]*pow(Wold[x][y][z]/f.W[x][y][z],
					f.kappa[x][y][z]-1.0);


      }
    }
  }



  //  halo_field(f.W, p);

  halo_field(Wfacex, p);
  halo_field(Wfacey, p);
  halo_field(Wfacez, p);

  // Section 3.4.4 -- pressure terms
  for(x = 1; x <= p.slicex; x++) {
    for(y = 1; y <= p.slicey; y++) {
      for(z = 0; z < p.Lz; z++) {

	qx = (f.V[0][x+1][y][z] - f.V[0][x][y][z])/p.dx;
	qy = (f.V[1][x][y+1][z] - f.V[1][x][y][z])/p.dx;
	qz = (f.V[2][x][y][((z+1)%p.Lz)] - f.V[2][x][y][z])/p.dx;
	
	divv = qx + qy + qz;

		
	// s = 0.5*(f.kappa[x] - 1.0)*divv*p.dt/2.0;
	// f.E[x] = f.E[x]*(1-s)/(1+s);
	

	// Is it 1.0 or 0.5? -- not the major problem I think
	// also note that dx=1.0 still get problems
	// so it is not some factor of dx somewhere
	// nor is it this particular expression:
	f.E[x][y][z] = f.E[x][y][z]
	  *exp(-1.0*(f.kappa[x][y][z] - 1.0)*divv*p.dt);

      }
    }
  }



  //  halo_field(f.E, p);



  // Last bit of 3.4.6
  for(x = 1; x <= p.slicex; x++) {
    for(y = 1; y <= p.slicey; y++) {
      for(z = 0; z < p.Lz; z++) {

	qx = (f.V[0][x][y][z] + f.V[0][x+1][y][z])
	  *(Wfacex[x+1][y][z] - Wfacex[x][y][z])*p.dt/(2.0*p.dx);
	
	qy = (f.V[1][x][y][z] + f.V[1][x][y+1][z])
	  *(Wfacey[x][y+1][z] - Wfacey[x][y][z])*p.dt/(2.0*p.dx);
	
	qz = (f.V[2][x][y][z] + f.V[2][x][y][((z+1)%p.Lz)])
	  *(Wfacez[x][y][((z+1)%p.Lz)] - Wfacez[x][y][z])*p.dt/(2.0*p.dx);
	
	gradv = (qx+qy+qz)/f.W[x][y][z];
	
	//	s = (f.kappa[x] - 1.0)*gradv/2.0;
	//	f.E[iix(x,y,z,p)] = f.E[iix(x,y,z,p)]*(1-s)/(1+s);
	f.E[x][y][z] = f.E[x][y][z]*exp(-1.0*(f.kappa[x][y][z]-1.0)*gradv);
      }
    }      
  }



  halo_field(f.E, p);
  
  
  
  // Clean up memory. Surely we don't need all these temporary arrays?

  
  free_field(p, Vdmid);
  free_field(p, Wold);

  free_field(p, phiav);

  free_vector(p, dxphi);

  free_field(p, Wfacex);
  free_field(p, Wfacey);
  free_field(p, Wfacez);

#endif // SCALAR
}




/* void evolve_uij(hydro_fields f, hydro_params p)
 *
 * Evolve the metric perturbations and conjugate momenta
 * using the field+fluid system as a source. The cutoff parameter (and
 * conditional compilation) was used to decouple the metric perturbations
 * after a certain time and is no longer in use. Eventually it should
 * be removed.
 */
#ifdef CUTOFF
void evolve_uij(hydro_fields f, hydro_params p, float cutoff) {
#else
void evolve_uij(hydro_fields f, hydro_params p) {

  const float cutoff = 1.0;
#endif

  int x, y, z, i;

  float ****Tij = make_tensor(p);


  // c = hbar = mpl = G = 1
  const float G = 1.0;


  stress_energy(f, p, Tij);

  for(i = 0; i < TENSOR_CPTS; i++) {
    for(x = 1; x <= p.slicex; x++) {
      for(y = 1; y <= p.slicey; y++) {
	for(z = 0; z < p.Lz; z++) {

	  f.uij[i][x][y][z] = f.uij[i][x][y][z] + p.dt*f.udotij[i][x][y][z];
	}
      }
    }
  

    halo_field(f.uij[i], p);



    for(x = 1; x <= p.slicex; x++) {
      for(y = 1; y <= p.slicey; y++) {
	for(z = 0; z < p.Lz; z++) {
	  f.udotij[i][x][y][z] = f.udotij[i][x][y][z] + p.dt*
	    (f.uij[i][x+1][y][z]
	     + f.uij[i][x][y+1][z]
	     + f.uij[i][x][y][((z+1)%p.Lz)]
	     - 6.0*f.uij[i][x][y][z]
	     + f.uij[i][x-1][y][z]
	     + f.uij[i][x][y-1][z]
	     + f.uij[i][x][y][((z-1+p.Lz)%p.Lz)])
	    /(p.dx*p.dx);
	  
	  
	  f.udotij[i][x][y][z] = f.udotij[i][x][y][z]
	    + p.dt*16.0*M_PI*G*Tij[i][x][y][z]*cutoff;

	}
      }
    }
  }

  free_tensor(p, Tij);
}

