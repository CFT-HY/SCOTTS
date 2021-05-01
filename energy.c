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
 * on the current core.
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
 * on the current core.
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
 * those on the current core.
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
  float traceTij;
  
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

#ifdef TRACEFREE
	// If trace free compiler flag then remove the trace. In single
	// precision the trace of Tij can cause trace of udot to become large and leak into
	// hdot:
	traceTij = (Tij[CPT_11][x][y][z] + Tij[CPT_22][x][y][z] + Tij[CPT_33][x][y][z])/3.;

	Tij[CPT_11][x][y][z] -= traceTij;
	Tij[CPT_22][x][y][z] -= traceTij;
	Tij[CPT_33][x][y][z] -= traceTij;
#endif // TRACEFREE
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

/** Compute total vorticity of temperature current 
 * (curl J)^2 on local core. Allows calculation of enstrophy of J.
 *
 */
float get_curlJ_tot(hydro_fields f, hydro_params p){
#ifndef SCALAR
  int x, y, z;

  float ****J = make_vector(p);
  float temp;

  float curlJ_tot = 0;
  float vol=p.dx*p.dx*p.dx;

  // Construct temperature current (J) (centered at cell)

  for(x = 1; x <= p.slicex; x++) {
    for(y = 1; y <= p.slicey; y++) {
      for(z = 0; z < p.Lz; z++) {

	J[0][x][y][z] = 0.5*(f.V[0][x][y][z] + f.V[0][x+1][y][z]
				)*f.T[x][y][z]*f.W[x][y][z];
	J[1][x][y][z] = 0.5*(f.V[1][x][y][z] + f.V[1][x][y+1][z]
				)*f.T[x][y][z]*f.W[x][y][z];
	J[2][x][y][z] = 0.5*(f.V[2][x][y][z] + f.V[2][x][y][(z+1)%p.Lz]
				)*f.T[x][y][z]*f.W[x][y][z];
      }
    }
  }

  halo_field(J[0], p);
  halo_field(J[1], p);
  halo_field(J[2], p);


  // Construct (curl J)^2.
  // Use centered first-order difference so all components live in the same place,
  // and we avoid generation of spurious vorticity.
  for(x = 1; x <= p.slicex; x++) {
    for(y = 1; y <= p.slicey; y++) {
      for(z = 0; z < p.Lz; z++) {
	temp = (J[2][x][y+1][z] - J[2][x][y-1][z]
			  - J[1][x][y][(z+1)%p.Lz]
			  + J[1][x][y][(z-1+p.Lz)%p.Lz])/(2*p.dx);
	curlJ_tot += temp*temp*vol;
	
	temp = (J[0][x][y][(z+1)%p.Lz]
			  - J[0][x][y][(z-1+p.Lz)%p.Lz]
			  - J[2][x+1][y][z] + J[2][x-1][y][z])/(2*p.dx);
	curlJ_tot += temp*temp*vol;

	temp = (J[1][x+1][y][z] - J[1][x-1][y][z]
			  - J[0][x][y+1][z] + J[0][x][y-1][z])/(2*p.dx);
	curlJ_tot += temp*temp*vol;
      }
    }
  }
  free_vector(p, J);
  return curlJ_tot;
#else
  return 0;
#endif //!SCALAR   
}

/** Compute total divergence of temperature current on local core
 * (div J)^2.
 *
 */
float get_divJ_tot(hydro_fields f, hydro_params p){
#ifndef SCALAR
  int x, y, z;
  float ****J = make_vector(p);
  float divJ_tot = 0;
  float temp;
  float vol=p.dx*p.dx*p.dx;
  // Construct temperature current (centered at cell)

  for(x = 1; x <= p.slicex; x++) {
    for(y = 1; y <= p.slicey; y++) {
      for(z = 0; z < p.Lz; z++) {

	J[0][x][y][z] = 0.5*(f.V[0][x][y][z] + f.V[0][x+1][y][z]
				)*f.T[x][y][z]*f.W[x][y][z];
	J[1][x][y][z] = 0.5*(f.V[1][x][y][z] + f.V[1][x][y+1][z]
				)*f.T[x][y][z]*f.W[x][y][z];
	J[2][x][y][z] = 0.5*(f.V[2][x][y][z] + f.V[2][x][y][(z+1)%p.Lz]
				)*f.T[x][y][z]*f.W[x][y][z];
      }
    }
  }

  halo_field(J[0], p);
  halo_field(J[1], p);
  halo_field(J[2], p);

  // Construct (div J)^2
  // Use centered first-order difference so all components live in the same place,
  for(x = 1; x <= p.slicex; x++) {
    for(y = 1; y <= p.slicey; y++) {
      for(z = 0; z < p.Lz; z++) {
	temp = (J[0][x+1][y][z] - J[0][x-1][y][z]
		+ J[1][x][y+1][z] - J[1][x][y-1][z]
		+ J[2][x][y][(z+1)%p.Lz] - J[2][x][y][(z-1+p.Lz)%p.Lz])/(2*p.dx);
	divJ_tot += temp*temp*vol;
      }
    }
  }

  free_vector(p, J);
  return divJ_tot;
  
#else
  return 0;
#endif //!SCALAR
}
