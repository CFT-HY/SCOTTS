/** @file energy.c
 *
 * Calculation of total energy and related quantities.
 */
#include "hydro.h"


/** Computes the total energy in the scalar field.
 *
 * Separately computes the kinetic, gradient and potential energy of
 * the scalar field and returns the total, summed over all sites.
 */
float field_energy(hydro_fields f, hydro_params p) {

  int x, y, z;

  float vol;

  float Etot = 0.0;

  float a = 0.0;
  float b = 0.0;
  float c = 0.0;

  vol = p.dx*p.dx*p.dx;

  for(x = 1; x <= p.slicex; x++) {
    for(y = 1; y <= p.slicey; y++) {
      for(z = 0; z < p.Lz; z++) {

	
	a += 0.5*f.pi_future[x][y][z]*f.pi_future[x][y][z]*vol;

	b += 0.5*((f.phi[x+1][y][z] - f.phi[x][y][z])/p.dx)
	  *((f.phi[x+1][y][z] - f.phi[x][y][z])/p.dx)*vol;
	
	b += 0.5*((f.phi[x][y+1][z] - f.phi[x][y][z])/p.dx)
	  *((f.phi[x][y+1][z] - f.phi[x][y][z])/p.dx)*vol;
	
	b += 0.5*((f.phi[x][y][(z+1)%p.Lz] - f.phi[x][y][z])/p.dx)
	  *((f.phi[x][y][(z+1)%p.Lz] - f.phi[x][y][z])/p.dx)*vol;

#ifndef SCALAR
	c += Vf(p, f.T[x][y][z], f.phi[x][y][z])*vol;
#else
	c += Vf(p, p.Tconst, f.phi[x][y][z])*vol;
#endif

      }
    }
  }
  
  float atot = reduce_sum(a, p);
  float btot = reduce_sum(b, p);
  float ctot = reduce_sum(c, p);

#ifdef SCALAR
  printf0(p,"Totals: momentum %g gradient %g potential %g\n", atot, btot, ctot);
#endif // SCALAR

  return a + b + c;
}



/** Computes the gradient energy in the scalar field.
 *
 * Note that this does _not_ currently sum over all sites, only those
 * on the current node.
 */
float gradient_energy_field(hydro_fields f, hydro_params p) {

  int x, y, z;

  float vol;

  float Etot = 0.0;


  vol = p.dx*p.dx*p.dx;

  for(x = 1; x <= p.slicex; x++) {
    for(y = 1; y <= p.slicey; y++) {
      for(z = 0; z < p.Lz; z++) {

	
	Etot += 0.5*((f.phi[x+1][y][z] - f.phi[x][y][z])/p.dx)
	  *((f.phi[x+1][y][z] - f.phi[x][y][z])/p.dx)*vol;
	
	Etot += 0.5*((f.phi[x][y+1][z] - f.phi[x][y][z])/p.dx)
	  *((f.phi[x][y+1][z] - f.phi[x][y][z])/p.dx)*vol;
	
	Etot += 0.5*((f.phi[x][y][(z+1)%p.Lz] - f.phi[x][y][z])/p.dx)
	  *((f.phi[x][y][(z+1)%p.Lz] - f.phi[x][y][z])/p.dx)*vol;
	
      }
    }
  }
  
  return Etot;
}

/** Computes the kinetic energy in the scalar field.
 *
 * Note that this does _not_ currently sum over all sites, only those
 * on the current node.
 */
float kinetic_energy_field(hydro_fields f, hydro_params p) {

  int x, y, z;

  float vol;

  float Etot = 0.0;


  vol = p.dx*p.dx*p.dx;

  for(x = 1; x <= p.slicex; x++) {
    for(y = 1; y <= p.slicey; y++) {
      for(z = 0; z < p.Lz; z++) {

	Etot += 0.5*f.pi_future[x][y][z]*f.pi_future[x][y][z]*vol;
	
      }
    }
  }
  
  return Etot;
}

/** Compute the total internal energy of the system.
 *
 * Total (field+fluid) energy. Borrowed directly from the
 * 1+1D spherical fortran code, which might explain the strange
 * way of performing the calculation.
 *
 * NB: This function currently does _not_ sum over all sites, only
 * those on the current node.
 */
float total_energy(hydro_fields f, hydro_params p) {

  int x, y, z;

  float vol;

  float Etot = 0.0;
  float restE = 0.0;
  float kinE = 0.0;
  float kinphi = 0.0;
  float grdphi = 0.0;

  vol = p.dx*p.dx*p.dx;

  for(x = 1; x <= p.slicex; x++) {
    for(y = 1; y <= p.slicey; y++) {
      for(z = 0; z < p.Lz; z++) {


#ifndef SCALAR	       
	// rest energy
	restE += (f.E[x][y][z]/f.W[x][y][z])*vol;
	
	// kinetic energy
	kinE += f.kappa[x][y][z]*(f.E[x][y][z]/f.W[x][y][z])
	  *(f.W[x][y][z]*f.W[x][y][z]-1.0)*vol;
#endif // SCALAR

		
	// momentum squared (scalar field kinetic energy)
	kinphi += 0.5*f.pi_future[x][y][z]*f.pi_future[x][y][z]*vol;
	
	// gradient term
	
	grdphi += 0.5*((f.phi[x+1][y][z] - f.phi[x][y][z])/p.dx)
	  *((f.phi[x+1][y][z] - f.phi[x][y][z])/p.dx)*vol;
	
	grdphi += 0.5*((f.phi[x][y+1][z] - f.phi[x][y][z])/p.dx)
	  *((f.phi[x][y+1][z] - f.phi[x][y][z])/p.dx)*vol;
	
	grdphi += 0.5*((f.phi[x][y][(z+1)%p.Lz] 
			  - f.phi[x][y][z])/p.dx)
	  *((f.phi[x][y][(z+1)%p.Lz] 
	     - f.phi[x][y][z])/p.dx)*vol;


      }
      
    }
  }

  Etot = (restE+kinE+kinphi+grdphi);
  
  return Etot;

}



/** Compute the total kinetic energy in the fluid.
 *
 * NB: This function does _not_ currently sum over all sites.
 */
float kinetic_energy_fluid(hydro_fields f, hydro_params p) {

  int x, y, z;

  float vol;

  float kinE = 0.0;

  vol = p.dx*p.dx*p.dx;

  for(x = 1; x <= p.slicex; x++) {
    for(y = 1; y <= p.slicey; y++) {
      for(z = 0; z < p.Lz; z++) {

#ifndef SCALAR
	// kinetic energy
	kinE += f.kappa[x][y][z]*(f.E[x][y][z]/f.W[x][y][z])
	  *(f.W[x][y][z]*f.W[x][y][z]-1.0)*vol;
#endif // SCALAR

      }
      
    }
  }
  
  return kinE;


}



/** Compute the rest energy in the fluid.
 *
 * NB: This function does _not_ currently sum over all sites.
 */
float rest_energy(hydro_fields f, hydro_params p) {

  int x, y, z;

  float vol;

  float restE = 0.0;

  vol = p.dx*p.dx*p.dx;

  for(x = 1; x <= p.slicex; x++) {
    for(y = 1; y <= p.slicey; y++) {
      for(z = 0; z < p.Lz; z++) {

#ifndef SCALAR
	// kinetic energy
	restE += (f.E[x][y][z]/f.W[x][y][z])*vol;
	  
#endif // SCALAR

      }
      
    }
  }

  
  return restE; // /vol;



}



/** Compute the total pressure (despite the name).
 *
 * NB: This function does _not_ currently sum over all sites.
 */
float avg_pressure(hydro_fields f, hydro_params p) {

  int x, y, z;

  float vol;

  float press;

  vol = p.dx*p.dx*p.dx;

  press = 0.0;

  for(x = 1; x <= p.slicex; x++) {
    for(y = 1; y <= p.slicey; y++) {
      for(z = 0; z < p.Lz; z++) {

#ifndef SCALAR
	// kinetic energy
	press += f.p[x][y][z]*vol;
#endif // SCALAR

      }
      
    }
  }

  return press;



}
   
			    

/** Compute total energy, i.e. 00 component of stress-energy.
 *
 * The zero-zero component of the stress energy tensor.  It's
 * basically just a different way of calculating total_energy(), and
 * therefore serves as a cross-check.
 */
float tzerozero(hydro_fields f, hydro_params p) {

  float total = 0.0;

  int x, y, z;

  float vol;

  vol = p.dx*p.dx*p.dx;

  for(x = 1; x <= p.slicex; x++) {
    for(y = 1; y <= p.slicey; y++) {
      for(z = 0; z < p.Lz; z++) {

	// d_mu phi d^nu phi
	// pi field if 00, otherwise grad mu, grad nu
	total += f.pi_future[x][y][z]*f.pi_future[x][y][z];

	// d_alpha phi d^alpha phi
	// (minus sign if 00, otherwise plus)
	total -= (0.5*f.pi_future[x][y][z]*f.pi_future[x][y][z]
		  - 0.125*((f.phi[x+1][y][z] - f.phi[x-1][y][z])/p.dx)
		  *((f.phi[x+1][y][z] - f.phi[x-1][y][z])/p.dx)
		  - 0.125*((f.phi[x][y+1][z] - f.phi[x][y-1][z])/p.dx)
		  *((f.phi[x][y+1][z] - f.phi[x][y-1][z])/p.dx)
		  - 0.125*((f.phi[x][y][(z+1)%p.Lz] 
			    - f.phi[x][y][(z-1+p.Lz)%p.Lz])/p.dx)
		  *((f.phi[x][y][(z+1)%p.Lz]
		     - f.phi[x][y][(z-1+p.Lz)%p.Lz])/p.dx));

#ifndef SCALAR
	// radiative fluid pressure
	// (minus sign if 00, otherwise plus)
	total -= (p.gdeg*f.T[x][y][z]*f.T[x][y][z]*f.T[x][y][z]*f.T[x][y][z]);
	
	// potential
	// (minus sign if 00, otherwise plus)
	total -= (-1.0*Vf(p, f.T[x][y][z], f.phi[x][y][z]));
#else
	total -= (-1.0*Vf(p, p.Tconst, f.phi[x][y][z]));
#endif // SCALAR



#ifndef SCALAR
	// fluid energy
	// (remember U_mu U_nu's at the end, no other sign)
	total += (4.0*p.gdeg
		  *f.T[x][y][z]*f.T[x][y][z]*f.T[x][y][z]*f.T[x][y][z]
		  - f.T[x][y][z]*VTf(p, f.T[x][y][z], f.phi[x][y][z]))
	  *f.W[x][y][z]*f.W[x][y][z];
	// U = (W, U1, U2, U3)
#endif // SCALAR

      }
    }
  }

  return total;


}



/** Compute sources of metric perturbations.
 *
 * Terms in the stress-energy tensor that are *LINEAR* in the metric,
 * and therefore source metric perturbations.
 */
void stress_energy(hydro_fields f, hydro_params p, float ****Tij) {

  int x, y, z;

  for(x = 1; x <= p.slicex; x++) {
    for(y = 1; y <= p.slicey; y++) {
      for(z = 0; z < p.Lz; z++) {
	

	Tij[CPT_11][x][y][z] = 0.0;
	Tij[CPT_21][x][y][z] = 0.0;
	Tij[CPT_31][x][y][z] = 0.0;
	Tij[CPT_22][x][y][z] = 0.0;
	Tij[CPT_32][x][y][z] = 0.0;
	Tij[CPT_33][x][y][z] = 0.0;

	if(p.gwsource != GW_FIELD) {
#ifndef SCALAR
	  // fluid bits
	  Tij[CPT_11][x][y][z] += (4.0*p.gdeg*f.T[x][y][z]
				   *f.T[x][y][z]*f.T[x][y][z]*f.T[x][y][z]
				   - f.T[x][y][z]*VTf(p, f.T[x][y][z], 
						      f.phi[x][y][z]))
	    *f.U[0][x][y][z]*f.U[0][x][y][z];
	  
	  Tij[CPT_21][x][y][z] += (4.0*p.gdeg*f.T[x][y][z]
				   *f.T[x][y][z]*f.T[x][y][z]*f.T[x][y][z]
				   - f.T[x][y][z]*VTf(p, f.T[x][y][z],
						      f.phi[x][y][z]))
	    *f.U[1][x][y][z]*f.U[0][x][y][z];
	  
	  Tij[CPT_31][x][y][z] += (4.0*p.gdeg*f.T[x][y][z]
				   *f.T[x][y][z]*f.T[x][y][z]*f.T[x][y][z]
				   - f.T[x][y][z]*VTf(p, f.T[x][y][z],
						      f.phi[x][y][z]))
	    *f.U[2][x][y][z]*f.U[0][x][y][z];
	  
	  Tij[CPT_22][x][y][z] += (4.0*p.gdeg*f.T[x][y][z]
				   *f.T[x][y][z]*f.T[x][y][z]*f.T[x][y][z]
				   - f.T[x][y][z]*VTf(p, f.T[x][y][z],
						      f.phi[x][y][z]))
	    *f.U[1][x][y][z]*f.U[1][x][y][z];
	  
	  Tij[CPT_32][x][y][z] += (4.0*p.gdeg*f.T[x][y][z]
				   *f.T[x][y][z]*f.T[x][y][z]*f.T[x][y][z]
				   - f.T[x][y][z]*VTf(p, f.T[x][y][z],
						      f.phi[x][y][z]))
	  *f.U[2][x][y][z]*f.U[1][x][y][z];
	  
	  Tij[CPT_33][x][y][z] += (4.0*p.gdeg*f.T[x][y][z]
				*f.T[x][y][z]*f.T[x][y][z]*f.T[x][y][z]
				   - f.T[x][y][z]*VTf(p, f.T[x][y][z],
						      f.phi[x][y][z]))
	    *f.U[2][x][y][z]*f.U[2][x][y][z];
#endif // SCALAR
	}


	
	if(p.gwsource != GW_FLUID) {
	  // Gradient bits
	  Tij[CPT_11][x][y][z] +=
	    0.25*((f.phi[x+1][y][z]
		   - f.phi[x-1][y][z])/p.dx)
	    *((f.phi[x+1][y][z] - f.phi[x-1][y][z])/p.dx);

	  Tij[CPT_21][x][y][z] += 
	    0.25*((f.phi[x+1][y][z]
		   - f.phi[x-1][y][z])/p.dx)
	    *((f.phi[x][y+1][z] - f.phi[x][y-1][z])/p.dx);
	  
	  Tij[CPT_31][x][y][z] +=
	    0.25*((f.phi[x][y][(z+1)%p.Lz]
		   - f.phi[x][y][(z-1+p.Lz)%p.Lz])/p.dx)
	    *((f.phi[x+1][y][z] - f.phi[x-1][y][z])/p.dx);
	  
	  Tij[CPT_22][x][y][z] += 
	    0.25*((f.phi[x][y+1][z]
		   - f.phi[x][y-1][z])/p.dx)
	    *((f.phi[x][y+1][z] - f.phi[x][y-1][z])/p.dx);
	  
	  Tij[CPT_32][x][y][z] +=
	    0.25*((f.phi[x][y][(z+1)%p.Lz]
		   - f.phi[x][y][(z-1+p.Lz)%p.Lz])/p.dx)
	    *((f.phi[x][y+1][z] - f.phi[x][y-1][z])/p.dx);
	  
	  Tij[CPT_33][x][y][z] +=
	    0.25*((f.phi[x][y][(z+1)%p.Lz]
		   - f.phi[x][y][(z-1+p.Lz)%p.Lz])/p.dx)
	    *((f.phi[x][y][(z+1)%p.Lz] - f.phi[x][y][(z-1+p.Lz)%p.Lz])/p.dx);
	}

      }
    }
  }

}


/** Compute energy density and store as a field.
 *
 * As for total_energy() but calculated on a per-lattice-site
 * basis and then stored in en.
 */
void energy_density(hydro_fields f, hydro_params p, float ***en) {

  int x, y, z;

  float vol;

  vol = p.dx*p.dx*p.dx;

  for(x = 1; x <= p.slicex; x++) {
    for(y = 1; y <= p.slicey; y++) {
      for(z = 0; z < p.Lz; z++) {

	en[x][y][z] = 0.0;


#ifndef SCALAR
	// rest energy
	en[x][y][z] += (f.E[x][y][z]/f.W[x][y][z])*vol;
	
	// kinetic energy
	en[x][y][z] += f.kappa[x][y][z]*(f.E[x][y][z]/f.W[x][y][z])
	  *(f.W[x][y][z]*f.W[x][y][z]-1.0)*vol;
#endif // SCALAR	

	// momentum squared (scalar field kinetic energy)
	en[x][y][z] += 0.5*f.pi_future[x][y][z]*f.pi_future[x][y][z]*vol;
	
	// gradient term
	
	en[x][y][z] += 0.125*((f.phi[x+1][y][z] - f.phi[x-1][y][z])/p.dx)
	  *((f.phi[x+1][y][z] - f.phi[x-1][y][z])/p.dx)*vol;
	
	en[x][y][z] += 0.125*((f.phi[x][y+1][z] - f.phi[x][y-1][z])/p.dx)
	  *((f.phi[x][y+1][z] - f.phi[x][y-1][z])/p.dx)*vol;
	
	en[x][y][z] += 0.125*((f.phi[x][y][(z+1)%p.Lz] 
			       - f.phi[x][y][(z-1+p.Lz)%p.Lz])/p.dx)
	  *((f.phi[x][y][(z+1)%p.Lz]
	     - f.phi[x][y][(z-1+p.Lz)%p.Lz])/p.dx)*vol;

	  

      }
      
    }
  }


}

/** Compute the total of b^2 on the lattice,  where b is the curl of the enthalpy
 * current, b = curl wU/\bar{w}.
 */
float get_btot(hydro_fields f, hydro_params p){
#ifndef SCALAR
  int x, y, z;
  float loc_tot_enth = 0;
  float mean_enth = 0;
  float ***Wnb = make_field(p);
  float b_tot = 0;
  
  for(x = 1; x <= p.slicex; x++) {
    for(y = 1; y <= p.slicey; y++) {
      for(z = 0; z < p.Lz; z++) {
	loc_tot_enth += f.p[x][y][z] + f.E[x][y][z]/f.W[x][y][z];
      }
    }
  }

  mean_enth = reduce_sum(loc_tot_enth, p)/(p.Lx*p.Ly*p.Lz);

  // Construct node centered W.
  for(x = 1; x <= p.slicex; x++) {
    for(y = 1; y <= p.slicey; y++) {
      for(z = 0; z < p.Lz; z++) {
	Wnb[x][y][z] = 0.125*(f.W[x-1][y][z]
			      + f.W[x][y][z]
			      + f.W[x][y-1][z]
			      + f.W[x-1][y-1][z]
			      + f.W[x][y][((z-1+p.Lz)%p.Lz)]
			      + f.W[x][y-1][((z+p.Lz-1)%p.Lz)]
			      + f.W[x-1][y][((z+p.Lz-1)%p.Lz)]
			      + f.W[x-1][y-1][((z+p.Lz-1)%p.Lz)]);
      }
    }
  }
  halo_field(Wnb, p);
  
  for(x = 1; x <= p.slicex; x++) {
    for(y = 1; y <= p.slicey; y++) {
      for(z = 0; z < p.Lz; z++) {
	b_tot += pow((f.Z[2][x][y+1][z]/Wnb[x][y+1][z]
		      - f.Z[2][x][y-1][z]/Wnb[x][y-1][z]
		      - f.Z[1][x][y][(z+1)%p.Lz]/Wnb[x][y][(z+1)%p.Lz]
		      + f.Z[1][x][y][(z-1+p.Lz)%p.Lz]
		      /Wnb[x][y][(z-1+p.Lz)%p.Lz]
		      )/(2*p.dx*mean_enth)
		     , 2);
	
	b_tot += pow((f.Z[0][x][y][(z+1)%p.Lz]/Wnb[x][y][(z+1)%p.Lz]
		      - f.Z[0][x][y][(z-1+p.Lz)%p.Lz]
		      /Wnb[x][y][(z-1+p.Lz)%p.Lz]
		      - f.Z[2][x+1][y][z]/Wnb[x+1][y][z]
		      + f.Z[2][x-1][y][z]/Wnb[x-1][y][z]
		      )/(2*p.dx*mean_enth)
		     , 2);

	b_tot += pow((f.Z[1][x+1][y][z]/Wnb[x+1][y][z]
		       - f.Z[1][x-1][y][z]/Wnb[x-1][y][z]
		       - f.Z[0][x][y+1][z]/Wnb[x][y+1][z]
		       + f.Z[0][x][y-1][z]/Wnb[x][y-1][z]
		       )/(2*p.dx*mean_enth)
		      , 2);	
      }
    }
  }

  free_field(p,Wnb);
  return b_tot;
#else
  return 0;

#endif //!SCALAR

}
			    
/** Compute the total of c^2 on the lattice,  where c is the divergence of the
 * enthalpy current, c = div wU/\bar{w}.
 */
float get_ctot(hydro_fields f, hydro_params p){
#ifndef SCALAR
  int x, y, z;
  float loc_tot_enth = 0;
  float mean_enth = 0;
  float ctot = 0;
  float ***Wnb = make_field(p);
  
  for(x = 1; x <= p.slicex; x++) {
    for(y = 1; y <= p.slicey; y++) {
      for(z = 0; z < p.Lz; z++) {
	loc_tot_enth += f.p[x][y][z] + f.E[x][y][z]/f.W[x][y][z];
      }
    }
  }

  mean_enth = reduce_sum(loc_tot_enth, p)/(p.Lx*p.Ly*p.Lz);

  // Construct node centered W.
  for(x = 1; x <= p.slicex; x++) {
    for(y = 1; y <= p.slicey; y++) {
      for(z = 0; z < p.Lz; z++) {
	Wnb[x][y][z] = 0.125*(f.W[x-1][y][z]
			      + f.W[x][y][z]
			      + f.W[x][y-1][z]
			      + f.W[x-1][y-1][z]
			      + f.W[x][y][((z-1+p.Lz)%p.Lz)]
			      + f.W[x][y-1][((z+p.Lz-1)%p.Lz)]
			      + f.W[x-1][y][((z+p.Lz-1)%p.Lz)]
			      + f.W[x-1][y-1][((z+p.Lz-1)%p.Lz)]);
      }
    }
  }
  halo_field(Wnb, p);
    
  for(x = 1; x <= p.slicex; x++) {
    for(y = 1; y <= p.slicey; y++) {
      for(z = 0; z < p.Lz; z++) {
	ctot += pow((f.Z[0][x+1][y][z]/Wnb[x+1][y][z]
		     - f.Z[0][x-1][y][z]/Wnb[x-1][y][z]
		     + f.Z[0][x][y+1][z]/Wnb[x][y+1][z]
		     - f.Z[0][x][y-1][z]/Wnb[x][y-1][z]
		     + f.Z[0][x][y][(z+1)%p.Lz]/Wnb[x][y][(z+1)%p.Lz]
		     - f.Z[0][x][y][(z-1+p.Lz)%p.Lz]
		     /Wnb[x][y][(z-1+p.Lz)%p.Lz]
		     )/(2*p.dx*mean_enth),2);
      }
    }
  }

  free_field(p,Wnb);

  return ctot;
#else
  return 0;

#endif //!SCALAR

}
