/* energy.c
 *
 * Calculation of total energy and related quantities.
 */
#include "hydro.h"


/* double field_energy(hydro_fields f, hydro_params p)
 *
 * Total energy in the fields: kinetic, gradient and potential.
 */
double field_energy(hydro_fields f, hydro_params p) {

  int x, y, z;

  double vol;

  double Etot = 0.0;


  vol = p.dx*p.dx*p.dx;

  for(x = 1; x <= p.slicex; x++) {
    for(y = 1; y <= p.slicey; y++) {
      for(z = 0; z < p.Lz; z++) {

	
	Etot += 0.5*f.pifull[x][y][z]*f.pifull[x][y][z]*vol;

	Etot += 0.5*((f.phi[x+1][y][z] - f.phi[x][y][z])/p.dx)
	  *((f.phi[x+1][y][z] - f.phi[x][y][z])/p.dx)*vol;
	
	Etot += 0.5*((f.phi[x][y+1][z] - f.phi[x][y][z])/p.dx)
	  *((f.phi[x][y+1][z] - f.phi[x][y][z])/p.dx)*vol;
	
	Etot += 0.5*((f.phi[x][y][(z+1)%p.Lz] - f.phi[x][y][z])/p.dx)
	  *((f.phi[x][y][(z+1)%p.Lz] - f.phi[x][y][z])/p.dx)*vol;
	
	Etot += Vf(p, f.T[x][y][z], f.phi[x][y][z])*vol;

      }
    }
  }
  
  return Etot;
}



/* double gradient_energy(hydro_fields f, hydro_params p)
 *
 * Field gradient energy only.
 */
double gradient_energy(hydro_fields f, hydro_params p) {

  int x, y, z;

  double vol;

  double Etot = 0.0;


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



/* double total_energy(hydro_fields f, hydro_params p)
 *
 * Total (field+fluid) energy. Borrowed directly from the
 * 1+1D spherical fortran code, which might explain the strange
 * way of performing the calculation.
 */
double total_energy(hydro_fields f, hydro_params p) {

  int x, y, z;

  double vol;

  double Etot = 0.0;
  double restE = 0.0;
  double kinE = 0.0;
  double kinphi = 0.0;
  double grdphi = 0.0;

  vol = p.dx*p.dx*p.dx;

  for(x = 1; x <= p.slicex; x++) {
    for(y = 1; y <= p.slicey; y++) {
      for(z = 0; z < p.Lz; z++) {

	       
	// rest energy
	restE += (f.E[x][y][z]/f.W[x][y][z])*vol;
	
	// kinetic energy
	kinE += f.kappa[x][y][z]*(f.E[x][y][z]/f.W[x][y][z])
	  *(f.W[x][y][z]*f.W[x][y][z]-1.0)*vol;

		
	// momentum squared (scalar field kinetic energy)
	kinphi += 0.5*f.pifull[x][y][z]*f.pifull[x][y][z]*vol;
	
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



/* double kinetic_energy(hydro_fields f, hydro_params p)
 *
 * Fluid kinetic energy only.
 */
double kinetic_energy(hydro_fields f, hydro_params p) {

  int x, y, z;

  double vol;

  double Etot = 0.0;
  double restE = 0.0;
  double kinE = 0.0;
  double kinphi = 0.0;
  double grdphi = 0.0;

  vol = p.dx*p.dx*p.dx;

  for(x = 1; x <= p.slicex; x++) {
    for(y = 1; y <= p.slicey; y++) {
      for(z = 0; z < p.Lz; z++) {

	// kinetic energy
	kinE += f.kappa[x][y][z]*(f.E[x][y][z]/f.W[x][y][z])
	  *(f.W[x][y][z]*f.W[x][y][z]-1.0)*vol;

      }
      
    }
  }

  Etot = kinE;
  
  return Etot; // /vol;



}



/* double rest_energy(hydro_fields f, hydro_params p)
 *
 * Fluid rest energy only.
 */
double rest_energy(hydro_fields f, hydro_params p) {

  int x, y, z;

  double vol;

  double Etot = 0.0;
  double restE = 0.0;
  double kinE = 0.0;
  double kinphi = 0.0;
  double grdphi = 0.0;

  vol = p.dx*p.dx*p.dx;

  for(x = 1; x <= p.slicex; x++) {
    for(y = 1; y <= p.slicey; y++) {
      for(z = 0; z < p.Lz; z++) {

	// kinetic energy
	kinE += f.kappa[x][y][z]*(f.E[x][y][z]/f.W[x][y][z])
	  *vol;

      }
      
    }
  }

  Etot = kinE;
  
  return Etot; // /vol;



}
    
			    

/* double tzerozero(hydro_fields f, hydro_params p)
 *
 * The zero-zero component of the stress energy tensor.
 * It's basically just a different way of calculating total_energy,
 * (see above) and therefore serves as a cross-check.
 */
double tzerozero(hydro_fields f, hydro_params p) {

  double total = 0.0;

  int x, y, z;

  double vol;

  vol = p.dx*p.dx*p.dx;

  for(x = 1; x <= p.slicex; x++) {
    for(y = 1; y <= p.slicey; y++) {
      for(z = 0; z < p.Lz; z++) {

	// d_mu phi d^nu phi
	// pi field if 00, otherwise grad mu, grad nu
	total += f.pifull[x][y][z]*f.pifull[x][y][z];

	// d_alpha phi d^alpha phi
	// (minus sign if 00, otherwise plus)
	total -= (0.5*f.pifull[x][y][z]*f.pifull[x][y][z]
		  - 0.125*((f.phi[x+1][y][z] - f.phi[x-1][y][z])/p.dx)
		  *((f.phi[x+1][y][z] - f.phi[x-1][y][z])/p.dx)
		  - 0.125*((f.phi[x][y+1][z] - f.phi[x][y-1][z])/p.dx)
		  *((f.phi[x][y+1][z] - f.phi[x][y-1][z])/p.dx)
		  - 0.125*((f.phi[x][y][(z+1)%p.Lz] 
			    - f.phi[x][y][(z-1+p.Lz)%p.Lz])/p.dx)
		  *((f.phi[x][y][(z+1)%p.Lz]
		     - f.phi[x][y][(z-1+p.Lz)%p.Lz])/p.dx));

	// radiative fluid pressure
	// (minus sign if 00, otherwise plus)
	total -= (p.gdeg*f.T[x][y][z]*f.T[x][y][z]*f.T[x][y][z]*f.T[x][y][z]);
	
	// potential
	// (minus sign if 00, otherwise plus)
	total -= (-1.0*Vf(p, f.T[x][y][z], f.phi[x][y][z]));

	// fluid energy
	// (remember U_mu U_nu's at the end, no other sign)
	total += (4.0*p.gdeg
		  *f.T[x][y][z]*f.T[x][y][z]*f.T[x][y][z]*f.T[x][y][z]
		  - f.T[x][y][z]*VTf(p, f.T[x][y][z], f.phi[x][y][z]))
	  *f.W[x][y][z]*f.W[x][y][z];
	// U = (W, U1, U2, U3)

      }
    }
  }

  return total;


}



/* void stress_energy(hydro_fields f, hydro_params p, double ****Tij)
 *
 * Terms in the stress-energy tensor that are *LINEAR* in the metric,
 * and therefore source metric perturbations.
 */
void stress_energy(hydro_fields f, hydro_params p, double ****Tij) {

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


/* void energy_density(hydro_fields f, hydro_params p, double ***en)
 *
 * As for total_energy, see above, but calculated on a per-lattice-site
 * basis and then stored in en.
 */
void energy_density(hydro_fields f, hydro_params p, double ***en) {

  int x, y, z;

  double vol;

  vol = p.dx*p.dx*p.dx;

  for(x = 1; x <= p.slicex; x++) {
    for(y = 1; y <= p.slicey; y++) {
      for(z = 0; z < p.Lz; z++) {

	en[x][y][z] = 0.0;

	// rest energy
	en[x][y][z] += (f.E[x][y][z]/f.W[x][y][z])*vol;
	
	// kinetic energy
	en[x][y][z] += f.kappa[x][y][z]*(f.E[x][y][z]/f.W[x][y][z])
	  *(f.W[x][y][z]*f.W[x][y][z]-1.0)*vol;
	
	// momentum squared (scalar field kinetic energy)
	en[x][y][z] += 0.5*f.pifull[x][y][z]*f.pifull[x][y][z]*vol;
	
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
    
			    
