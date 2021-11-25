/** @file evolve.c
 *
 * Everything to do with evolving fields and fluids with the
 * equations of motion.
 *
 * Contributors:
 * - 2010-2017 David Weir
 * - 2018-     Daniel Cutting
 */
#include "hydro.h"





/** Evolve the scalar field forward one timestep.
 * 
 * Taken more or less straight from the 1+1 spherical fortran code.
 * Uses leapfrog and Crank-Nicolson method to perform the update.
 */
void evolve_field(hydro_fields f, hydro_params p) {
  
  int x, y, z;

  float pi_old, s;
  
  // Move conjugate momentum (leapfrog)
  for(x = 1; x <= p.slicex; x++) {
    for(y = 1; y <= p.slicey; y++) {
      for(z = 0; z < p.Lz; z++) {
	
	
	pi_old = f.pi[x][y][z];


	/* First-order viscosity term:
	 * Factor of 0.5 comes from Crank-Nicolson method and
	 * leapfrog: pi lives on boundary between timesteps but need
	 * pi on timestep so need to average.  Get a term with
	 * (1+s)*pi from last timestep on RHS and have (1-s)*pi on
	 * current step on LHS, then divide through by (1-s).
	 */
#ifdef DIMENSIONLESS
#ifndef SCALAR
	s = -0.5*p.dt*(p.C*f.phi[x][y][z]*f.phi[x][y][z]
		       /f.T[x][y][z])*f.W[x][y][z];
#else
	s = -0.5*p.dt*(p.C*f.phi[x][y][z]*f.phi[x][y][z]
		       /p.Tconst);
#endif //!SCALAR
#else
#ifndef SCALAR
	s = -0.5*p.dt*p.C*f.W[x][y][z];
#else
	s = -0.5*p.dt*p.C;
#endif //!SCALAR
#endif // DIMENSIONLESS
	f.pi[x][y][z] = (1+s)*f.pi[x][y][z]/(1-s);
	
#ifndef SCALAR
	/* Gradient term:
	 * Gradients and velocities live on cell boundaries, but we
	 * want the values in the body of the cell so we average the
	 * two boundaries.
	 */

	
	f.pi[x][y][z] = f.pi[x][y][z] \
	  + 2*s*(
		 0.5*(f.V[0][x][y][z]*(f.phi[x][y][z]
				       - f.phi[x-1][y][z])
		      + f.V[0][x+1][y][z]*(f.phi[x+1][y][z]
					   - f.phi[x][y][z]))/p.dx
		 
		 + 0.5*(f.V[1][x][y][z]*(f.phi[x][y][z]
					 - f.phi[x][y-1][z])
			+ f.V[1][x][y+1][z]*(f.phi[x][y+1][z]
					     - f.phi[x][y][z]))/p.dx
		 
		 + 0.5*(f.V[2][x][y][z]*(f.phi[x][y][z]
					 - f.phi[x][y][((z-1+p.Lz)%p.Lz)])
			+ f.V[2][x][y][((z+1)%p.Lz)]*(f.phi[x][y][((z+1)%p.Lz)]
						      - f.phi[x][y][z]))/p.dx
		 )/(1-s);

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
	- Vdf(p, p.Tconst, f.phi[x][y][z]))/(1-s);
#endif // !SCALAR	


	   
	/* In the leapfrog/Verlet method, pi lives half a timestep behind
	 * phi. Evolving by an additional half-interval between the new and
	 * old pi's gives an approximation to the value of pi on the same
	 * timestep as phi. This is useful for outputting energies, etc., but
	 * is not used in the evolution code.
	 */
	f.pi_future[x][y][z] = f.pi[x][y][z] + 0.5*(f.pi[x][y][z] - pi_old);
      }
    }
  }

  // Move field (leapfrog)
  for(x = 1; x <= p.slicex; x++) {
    for(y = 1; y <= p.slicey; y++) {
      for(z = 0; z < p.Lz; z++) {
	// (We need to keep phi_old for averaging things on the half timestep)
	f.phi_old[x][y][z] = f.phi[x][y][z];
	
	f.phi[x][y][z] = f.phi[x][y][z] + p.dt*f.pi[x][y][z];

      }
    }
  }

  halo_field(f.pi, p);
  halo_field(f.pi_future, p);
  
  halo_field(f.phi, p);
  halo_field(f.phi_old, p);
  
}




/** Partially evolve the hydro forward one step following chapters 2 and 3 in
 * Wilson and Mathews. In this function we:
 * -# Update \f$ E\f$ and \f$ Z \f$ with Field-fluid interaction terms 
 *    (potential and dissipative term)
 *
 * Grid velocity is 0, lapse vector 0 and flat space.
 *
 * Basically straight from fortran. Original comment on full update:
 *
 *   --  Does hydro except transport.  E  and  Z evolved.  v  and W
 *       obtained.  Uses artificial viscosity  Q,  adjusted with  Cav.
 *
 *   --  The order chosen is not unique, but is what I used in
 *       dissertation.  Should experiment with different orders,
 *       especially with position of field-fluid interaction,
 *       now placed first (otherwise order is from CW).
 *
 */
void evolve_hydro_fieldfluid(hydro_fields f, hydro_params p) {

#ifndef SCALAR

  int x, y, z;
  
  float ***Vdmid = make_field(p); // dV/dphi between timesteps in zone.
  float ***phiav = make_field(p); // phi between timesteps in zone.

  //d(phi)/dx_i between timesteps on appropriate spatial boundary.
  float ****dxphi = make_vector(p); 

  float gpi, s;

#ifdef DIMENSIONLESS
  float phinb, Tnb;
#endif //DIMENSIONLESS
  float vdnb, pinb, wnb, dxphinb0, dxphinb1, dxphinb2;
  float velxnb, velynb, velznb;

  
  // Find quantities needed on half timesteps but defined on whole steps
  for(x = 1; x <= p.slicex; x++) {
    for(y = 1; y <= p.slicey; y++) {
      for(z = 0; z < p.Lz; z++) {


	
	phiav[x][y][z] = 0.5*(f.phi_old[x][y][z] + f.phi[x][y][z]);
	
	dxphi[0][x][y][z] =  0.5*(f.phi_old[x+1][y][z]
				  + f.phi[x+1][y][z]
				  - f.phi_old[x][y][z] - f.phi[x][y][z])/p.dx;
	
	dxphi[1][x][y][z] =  0.5*(f.phi_old[x][y+1][z]
				  + f.phi[x][y+1][z]
				  - f.phi_old[x][y][z] - f.phi[x][y][z])/p.dx;
	
	dxphi[2][x][y][z] =  0.5*(f.phi_old[x][y][((z+1)%p.Lz)]
				  + f.phi[x][y][((z+1)%p.Lz)]
				  - f.phi_old[x][y][z] - f.phi[x][y][z])/p.dx;
	
      }
    }
  }

  halo_field(dxphi[0], p);
  halo_field(dxphi[1], p);
  halo_field(dxphi[2], p);

  Vdpot(p, f.T, phiav, Vdmid);


  halo_field(Vdmid, p);

  // Field-fluid interaction through friction term and potential term:
  
  for(x = 1; x <= p.slicex; x++) {
    for(y = 1; y <= p.slicey; y++) {
      for(z = 0; z < p.Lz; z++) {

	// Averaging to get variables at the corners of the zone. 
	
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

#ifdef DIMENSIONLESS
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
#endif //DIMENSIONLESS
	
	
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

	// Construct node based velocity components.
	velxnb = 0.25*(f.V[0][x][y][z] + f.V[0][x][y-1][z]
		       + f.V[0][x][y][(z-1+p.Lz)%p.Lz] + f.V[0][x][y-1][(z-1+p.Lz)%p.Lz]);
	
	velynb = 0.25*(f.V[1][x][y][z] + f.V[1][x-1][y][z]
		       + f.V[1][x][y][(z-1+p.Lz)%p.Lz] + f.V[1][x-1][y][(z-1+p.Lz)%p.Lz]);

	velznb = 0.25*(f.V[2][x][y][z] + f.V[2][x-1][y][z]
		       + f.V[2][x][y-1][z] + f.V[2][x-1][y-1][z]);

	
#ifdef DIMENSIONLESS	
	// this is the one that causes trouble       
	f.Z[0][x][y][z] = f.Z[0][x][y][z] - p.dt*(p.C*phinb*phinb/Tnb)*wnb*pinb*dxphinb0;
	f.Z[1][x][y][z] = f.Z[1][x][y][z] - p.dt*(p.C*phinb*phinb/Tnb)*wnb*pinb*dxphinb1;
	f.Z[2][x][y][z] = f.Z[2][x][y][z] - p.dt*(p.C*phinb*phinb/Tnb)*wnb*pinb*dxphinb2;
    
	f.Z[1][x][y][z] = f.Z[1][x][y][z] 
	  - p.dt*(p.C*phinb*phinb/Tnb)*wnb*(velynb*dxphinb1
					    + velxnb*dxphinb0
					    + velznb*dxphinb2)*dxphinb1;
	
	f.Z[0][x][y][z] = f.Z[0][x][y][z] 
	  - p.dt*(p.C*phinb*phinb/Tnb)*wnb*(velynb*dxphinb1
					    + velxnb*dxphinb0
					    + velznb*dxphinb2)*dxphinb0;    
    
	f.Z[2][x][y][z] = f.Z[2][x][y][z]
	  - p.dt*(p.C*phinb*phinb/Tnb)*wnb*(velynb*dxphinb1
					    + velxnb*dxphinb0
					    + velznb*dxphinb2)*dxphinb2;
	
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
#else
	f.Z[0][x][y][z] = f.Z[0][x][y][z] - p.dt*p.C*wnb*pinb*dxphinb0;
	f.Z[1][x][y][z] = f.Z[1][x][y][z] - p.dt*p.C*wnb*pinb*dxphinb1;
	f.Z[2][x][y][z] = f.Z[2][x][y][z] - p.dt*p.C*wnb*pinb*dxphinb2;
	
	f.Z[1][x][y][z] = f.Z[1][x][y][z] 
	  - p.dt*p.C*wnb*(velynb*dxphinb1
			  + velxnb*dxphinb0
			  + velznb*dxphinb2)*dxphinb1;
	
	f.Z[0][x][y][z] = f.Z[0][x][y][z] 
	  - p.dt*p.C*wnb*(velynb*dxphinb1
			  + velxnb*dxphinb0
			  + velznb*dxphinb2)*dxphinb0;    
	
	f.Z[2][x][y][z] = f.Z[2][x][y][z]
	  - p.dt*p.C*wnb*(velynb*dxphinb1
			  + velxnb*dxphinb0
			  + velznb*dxphinb2)*dxphinb2;
	
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
      
	f.E[x][y][z] = f.E[x][y][z] + p.dt*(p.C*gpi*gpi + gpi*Vdmid[x][y][z]);
#endif //DIMENSIONLESS
      }
    }
  }

  // Technically unnecessary as Z only accessed out of place in transport_Z_dir
  // and donor_Z_dir.
  //halo_field(f.Z[0], p);
  //halo_field(f.Z[1], p);
  //halo_field(f.Z[2], p);
  halo_field(f.E, p);  
  
  // Clean up memory. Surely we don't need all these temporary arrays?

  
  free_field(p, Vdmid);

  free_field(p, phiav);

  free_vector(p, dxphi);


#endif // SCALAR
}



/** Partially evolve the hydro forward one step following chapters 2 and 3 in
 * Wilson and Mathews. In this function we:
 * -# Update \f$ Z \f$ with pressure acceleration
 *
 * Grid velocity is 0, lapse vector 0 and flat space.
 *
 * Basically straight from fortran. Original comment on full update:
 *
 *   --  Does hydro except transport.  E  and  Z evolved.  v  and W
 *       obtained.  Uses artificial viscosity  Q,  adjusted with  Cav.
 *
 *   --  The order chosen is not unique, but is what I used in
 *       dissertation.  Should experiment with different orders,
 *       especially with position of field-fluid interaction,
 *       now placed first (otherwise order is from CW).
 */
void evolve_hydro_pressureacceleration(hydro_fields f, hydro_params p) {

#ifndef SCALAR

  int x, y, z;
  
  float p_bar_x_plus, p_bar_x_minus;
  float p_bar_y_plus, p_bar_y_minus;
  float p_bar_z_plus, p_bar_z_minus;




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
			+ f.p[x][y][((z+p.Lz-1)%p.Lz)]
			+ f.p[x-1][y][((z+p.Lz-1)%p.Lz)]
			)/4.0;
	
	
	p_bar_y_minus = (  f.p[x][y-1][z]
			 + f.p[x-1][y-1][z] 
			 + f.p[x][y-1][((z+p.Lz-1)%p.Lz)]
			 + f.p[x-1][y-1][((z+p.Lz-1)%p.Lz)]
			 )/4.0;
	
	
	p_bar_z_plus = (f.p[x][y][z]
			+ f.p[x-1][y][z] 
			+ f.p[x][y-1][z]
			+ f.p[x-1][y-1][z]
		    )/4.0;
	
	
	p_bar_z_minus = (f.p[x][y][((z+p.Lz-1)%p.Lz)]
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

  // Technically unnecessary as Z only accessed out of place in transport_Z_dir
  // and donor_Z_dir.
  //halo_field(f.Z[0], p);
  //halo_field(f.Z[1], p);
  //halo_field(f.Z[2], p); 

#endif // SCALAR
  
}

/** Partially evolve the hydro forward one step following chapters 2 and 3 in
 * Wilson and Mathews. In this function we:
 * -# Update covariant 4 velocity spatial terms \f$ U_i \f$.
 * -# Update contravariant 3 velocity \f$ V^i \f$.
 * -# Update Lorentz factor \f$ W \f$.
 *
 * Grid velocity is 0, lapse vector 0 and flat space.
 *
 * Basically straight from fortran. Original comment on full update:
 *
 *   --  Does hydro except transport.  E  and  Z evolved.  v  and W
 *       obtained.  Uses artificial viscosity  Q,  adjusted with  Cav.
 *
 *   --  The order chosen is not unique, but is what I used in
 *       dissertation.  Should experiment with different orders,
 *       especially with position of field-fluid interaction,
 *       now placed first (otherwise order is from CW).
 */
void evolve_hydro_velocities(hydro_fields f, hydro_params p) {

#ifndef SCALAR

  int x, y, z;
  

  float utildex, utildey, utildez;
  
  float sigmabar;

  float ubarx, ubary, ubarz;

  // update velocity v; denominator is W&M eq (2.85)
  // but note grid is Eulerian, see readme for conversion 
  // between gamma and kappa and D.
  // Section 3.4.5, equations 3.57, 3.58
  // Updates spatial comps of covariant 4 velocity

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
  // Update contravariant 3-velocity (lives on faces)
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
	
	f.Wfacex[x][y][z] = sqrt(1.0 + ubarx*ubarx + ubary*ubary + ubarz*ubarz);
	
	f.V[0][x][y][z] = ubarx/f.Wfacex[x][y][z];
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

	f.Wfacey[x][y][z] = sqrt(1.0 + ubarx*ubarx + ubary*ubary + ubarz*ubarz);

	f.V[1][x][y][z] = ubary/f.Wfacey[x][y][z];
	
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
	
	f.Wfacez[x][y][z] = sqrt(1.0 
			       + ubarx*ubarx
			       + ubary*ubary
			       + ubarz*ubarz);
    
	f.V[2][x][y][z] = ubarz/f.Wfacez[x][y][z];
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
  // Update Lorentz factor f.W, then update E with first pressure work term.
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
     
	f.Wold[x][y][z] = f.W[x][y][z];
	f.W[x][y][z] = sqrt(1.0 
			    + utildex*utildex 
			    + utildey*utildey 
			    + utildez*utildez);
      }
    }
  }


  halo_field(f.Wold, p);
  halo_field(f.W, p);

  halo_field(f.Wfacex, p);
  halo_field(f.Wfacey, p);
  halo_field(f.Wfacez, p);

    

#endif // SCALAR
}

/** Partially evolve the hydro forward one step following chapters 2 and 3 in
 * Wilson and Mathews. In this function we:
 * -# Update \f$ E \f$ with PdV work terms.
 *
 * Grid velocity is 0, lapse vector 0 and flat space.
 *
 * Basically straight from fortran. Original comment on full update:
 *
 *   --  Does hydro except transport.  E  and  Z evolved.  v  and W
 *       obtained.  Uses artificial viscosity  Q,  adjusted with  Cav.
 *
 *   --  The order chosen is not unique, but is what I used in
 *       dissertation.  Should experiment with different orders,
 *       especially with position of field-fluid interaction,
 *       now placed first (otherwise order is from CW).
 */
void evolve_hydro_pressurework(hydro_fields f, hydro_params p) {

#ifndef SCALAR

  int x, y, z;

  float qx, qy, qz, gradv;
  float divv;
  float Qsq;
  float Qsum;
  float epsilon = 1e-20;
  

  // Section 3.4.6
  // Update Lorentz factor f.W, then update E with first pressure work term.
  for(x = 1; x <= p.slicex; x++) {
    for(y = 1; y <= p.slicey; y++) {
      for(z = 0; z < p.Lz; z++) {
		
	        
	// Top of p90 ch 3
	f.E[x][y][z] = f.E[x][y][z]*pow(f.Wold[x][y][z]/f.W[x][y][z],
					f.kappa[x][y][z]-1.0);


      }
    }
  }

  // Artificial viscosity update to E which is proportional to dW/dt.
  if ((p.kq > 0) || (p.kl > 0)){
    for(x = 1; x <= p.slicex; x++) {
      for(y = 1; y <= p.slicey; y++) {
	for(z = 0; z < p.Lz; z++) {
	  Qsq = (f.Q[0][x][y][z]*f.Q[0][x][y][z] + f.Q[1][x][y][z]*f.Q[1][x][y][z] +
		 f.Q[2][x][y][z]*f.Q[2][x][y][z]);
	  Qsum = (f.Q[0][x][y][z] + f.Q[1][x][y][z] + f.Q[2][x][y][z]+ epsilon);
	  f.E[x][y][z] = f.E[x][y][z] - (f.W[x][y][z] - f.Wold[x][y][z])*Qsq/Qsum;
	}
      }
    }
  }

	


  // Section 3.4.4 -- Update E with another pressure work term
  for(x = 1; x <= p.slicex; x++) {
    for(y = 1; y <= p.slicey; y++) {
      for(z = 0; z < p.Lz; z++) {

	qx = (f.V[0][x+1][y][z] - f.V[0][x][y][z])/p.dx;
	qy = (f.V[1][x][y+1][z] - f.V[1][x][y][z])/p.dx;
	qz = (f.V[2][x][y][((z+1)%p.Lz)] - f.V[2][x][y][z])/p.dx;
	
	divv = qx + qy + qz;
			
	f.E[x][y][z] = f.E[x][y][z]
	  *exp(-1.0*(f.kappa[x][y][z] - 1.0)*divv*p.dt);

      }
    }
  }



  //  halo_field(f.E, p);



  // Last bit of 3.4.6 -- last pressure work term
  for(x = 1; x <= p.slicex; x++) {
    for(y = 1; y <= p.slicey; y++) {
      for(z = 0; z < p.Lz; z++) {

	qx = (f.V[0][x][y][z] + f.V[0][x+1][y][z])
	  *(f.Wfacex[x+1][y][z] - f.Wfacex[x][y][z])*p.dt/(2.0*p.dx);
	
	qy = (f.V[1][x][y][z] + f.V[1][x][y+1][z])
	  *(f.Wfacey[x][y+1][z] - f.Wfacey[x][y][z])*p.dt/(2.0*p.dx);
	
	qz = (f.V[2][x][y][z] + f.V[2][x][y][((z+1)%p.Lz)])
	  *(f.Wfacez[x][y][((z+1)%p.Lz)] - f.Wfacez[x][y][z])*p.dt/(2.0*p.dx);
	
	gradv = (qx+qy+qz)/f.W[x][y][z];
	
	f.E[x][y][z] = f.E[x][y][z]*exp(-1.0*(f.kappa[x][y][z]-1.0)*gradv);
      }
    }      
  }

  

  halo_field(f.E, p);
  
#endif // SCALAR
}

/** Partially evolve the hydro forward one step. In this function we update the
 * artificial viscosity "scalar". We follow the method described in Anninos and
 * Fragile 2002 for the artificial viscosity. We update the energy and momentum
 * using the artificial viscosity terms, apart from the term for E which
 * includes a dW/dt piece. This is done in evolve_hydropressurework.
 *
 * Set \Delta l = p.dx from Anninos and Fragile.
 *
 * Grid velocity is 0, lapse vector 0 and flat space.
 */
void evolve_hydro_artviscosity(hydro_fields f, hydro_params p)
{

#ifndef SCALAR
  int x, y, z;
  int dir, dx, dy, dz;
  float*** dV = make_field(p);
  float phi;
  // Speed of sound.
  // This is only a constant in BAG, needs to be implemented as a
  // function for general eos.
  float c_s = 1 / sqrt(3);

  float Qnbplus;
  float Qnbminus;

  for (dir = 0; dir < 3; dir++) {
    dx = 0;
    dy = 0;
    dz = 0;
    if (dir == 0) {
      dx = 1;
    }
    if (dir == 1) {
      dy = 1;
    }
    if (dir == 2) {
      dz = 1;
    }

    for (x = 1; x <= p.slicex; x++) {
      for (y = 1; y <= p.slicey; y++) {
        for (z = 0; z < p.Lz; z++) {
          dV[x][y][z]
              = (f.V[dir][x + dx][y + dy][(z + dz) % p.Lz] - f.V[dir][x][y][z]);
        }
      }
    }

    halo_field(dV, p);

    for (x = 1; x <= p.slicex; x++) {
      for (y = 1; y <= p.slicey; y++) {
        for (z = 0; z < p.Lz; z++) {
          if (dV[x][y][z] > 0) {
            f.Q[dir][x][y][z] = 0;
          } else {
            if ((dV[x + dx][y + dy][(z + dz) % p.Lz] * dV[x][y][z] < 0)
                || (dV[x - dx][y - dy][(z - dz + p.Lz) % p.Lz] * dV[x][y][z] < 0)
                || (f.V[dir][x + dx][y + dy][(z + dz) % p.Lz] * f.V[dir][x][y][z]
                    < 0)) {
              phi = 0;
            } else {
              phi = 0.5;
            }
            // Possibly phi in else statement should equal 0.5 or 1.0? Or should
            // we use something like minmod?

            // In cosmos++ they have W||Q|| term as
            // well as E*kappa piece.
            f.Q[dir][x][y][z]
                = (f.E[x][y][z] * f.kappa[x][y][z] * dV[x][y][z] * (1 - phi)
                   * (p.kq * dV[x][y][z] * (1 - phi) - p.kl * c_s));
          }
        }
      }
    }
  }
  halo_field(f.Q[0], p);
  halo_field(f.Q[1], p);
  halo_field(f.Q[2], p);

  // Now update E and Z with viscosity terms that don't involve dW/dt.
  for (x = 1; x <= p.slicex; x++) {
    for (y = 1; y <= p.slicey; y++) {
      for (z = 0; z < p.Lz; z++) {

        // First E.
        // x derivative
        f.E[x][y][z] = f.E[x][y][z]
                       - f.Q[0][x][y][z]
                             * (f.V[0][x + 1][y][z] * f.Wfacex[x + 1][y][z]
                                - f.V[0][x][y][z] * f.Wfacex[x][y][z])
                             * p.dt / p.dx;
        // y derivative
        f.E[x][y][z] = f.E[x][y][z]
                       - f.Q[1][x][y][z]
                             * (f.V[1][x][y + 1][z] * f.Wfacey[x][y + 1][z]
                                - f.V[1][x][y][z] * f.Wfacey[x][y][z])
                             * p.dt / p.dx;
        // z derivative
        f.E[x][y][z] = f.E[x][y][z]
                       - f.Q[2][x][y][z]
                             * (f.V[2][x][y][(z + 1) % p.Lz]
                                    * f.Wfacez[x][y][(z + 1) % p.Lz]
                                - f.V[2][x][y][z] * f.Wfacez[x][y][z])
                             * p.dt / p.dx;

        // Now Z.
        // Zx (dQx/dx term)
        Qnbplus = (f.Q[0][x][y][z] + f.Q[0][x][y - 1][z]
                   + f.Q[0][x][y][(z - 1 + p.Lz) % p.Lz]
                   + f.Q[0][x][y - 1][(z - 1 + p.Lz) % p.Lz]);
        Qnbminus = (f.Q[0][x - 1][y][z] + f.Q[0][x - 1][y - 1][z]
                    + f.Q[0][x - 1][y][(z - 1 + p.Lz) % p.Lz]
                    + f.Q[0][x - 1][y - 1][(z - 1 + p.Lz) % p.Lz]);

        f.Z[0][x][y][z]
            = f.Z[0][x][y][z] - 0.25 * (Qnbplus - Qnbminus) * p.dt / p.dx;

        // Zy (dQy/dy term)
        Qnbplus = (f.Q[1][x][y][z] + f.Q[1][x - 1][y][z]
                   + f.Q[1][x][y][(z - 1 + p.Lz) % p.Lz]
                   + f.Q[1][x - 1][y][(z - 1 + p.Lz) % p.Lz]);
        Qnbminus = (f.Q[1][x][y - 1][z] + f.Q[1][x - 1][y - 1][z]
                    + f.Q[1][x][y - 1][(z - 1 + p.Lz) % p.Lz]
                    + f.Q[1][x - 1][y - 1][(z - 1 + p.Lz) % p.Lz]);

        f.Z[1][x][y][z]
            = f.Z[1][x][y][z] - 0.25 * (Qnbplus - Qnbminus) * p.dt / p.dx;

        // Zz (dQz/dz term)
        Qnbplus = (f.Q[2][x][y][z] + f.Q[2][x - 1][y][z] + f.Q[2][x][y - 1][z]
                   + f.Q[2][x - 1][y - 1][z]);
        Qnbminus = (f.Q[2][x][y][(z - 1 + p.Lz) % p.Lz]
                    + f.Q[2][x - 1][y][(z - 1 + p.Lz) % p.Lz]
                    + f.Q[2][x][y - 1][(z - 1 + p.Lz) % p.Lz]
                    + f.Q[2][x - 1][y - 1][(z - 1 + p.Lz) % p.Lz]);

        f.Z[2][x][y][z]
            = f.Z[2][x][y][z] - 0.25 * (Qnbplus - Qnbminus) * p.dt / p.dx;
      }
    }
  }
  halo_field(f.E, p);

  free_field(p, dV);

#endif // SCALAR
}

/** Evolve the metric perturbations and conjugate momenta using
 * the field+fluid system as a source.
 *
 * Calculation of the stress energy tensor \f$ T^{ij} \f$ is performed in
 * stress_energy().
 *
 */
void evolve_uij(hydro_fields f, hydro_params p) {

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
	    + p.dt*16.0*M_PI*G*Tij[i][x][y][z];

	}
      }
    }
  }

  free_tensor(p, Tij);
}

