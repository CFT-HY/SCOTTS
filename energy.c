/** @file energy.c
 *
 * Calculation of total energy and related quantities.
 */
#include "hydro.h"



/** Calculates sum of various components of the total energy, split into
 *  symmetric and broken phase. Returns as an array of floats.
 *  
 *  Sets energies array to be as follows in order:
 *  rest energy, fluid kinetic energy, kinetic energy scalar, gradient energy scalar, potential
 *  energy scalar.
 *  Each energy component is a pair of floats, split into total in symmetric and
 *  broken phase, in that order. 
 *  Note that this function does not sum over all sites.
 */
void calculate_energies(hydro_fields f, hydro_params p, float *energies) {

  float rest_E_symm = 0;
  float rest_E_broken = 0;
  float kin_fluid_symm = 0;
  float kin_fluid_broken = 0;
  float kin_phi_symm = 0;
  float kin_phi_broken = 0;
  float grad_phi_symm = 0;
  float grad_phi_broken = 0;
  float pot_phi_symm = 0;
  float pot_phi_broken = 0;

  float vol = p.dx*p.dx*p.dx;
  
  int x, y, z;

  float phi_broken;
#ifdef BAG
  phi_broken = p.phi_0/2.;
#else
#ifdef SCALAR
  phi_broken =  (p.alpha*p.Tconst
		 + sqrt((p.alpha*p.Tconst)*(p.alpha*p.Tconst)
			- 4.0*p.lambda*p.gamma
			*(p.Tconst*p.Tconst - p.T0*p.T0))
		 )/(2.0*p.lambda);

#endif // SCALAR
#endif // BAG
  
  for(x = 1; x <= p.slicex; x++) {
    for(y = 1; y <= p.slicey; y++) {
      for(z = 0; z < p.Lz; z++) {
#if !defined(SCALAR) && !defined(BAG)
	phi_broken =  (p.alpha*f.T[x][y][z]
		       + sqrt((p.alpha*p.Tconst)*(p.alpha*f.T[x][y][z])
			      - 4.0*p.lambda*p.gamma
			      *(f.T[x][y][z]*f.T[x][y][z] - p.T0*p.T0))
		       )/(2.0*p.lambda);
#endif // !SCALAR && !BAG

	if (f.phi[x][y][z] < phi_broken/2){
#ifndef SCALAR
	  rest_E_symm += (f.E[x][y][z]/f.W[x][y][z])*vol;

	  kin_fluid_symm += f.kappa[x][y][z]*(f.E[x][y][z]/f.W[x][y][z])
	    *(f.W[x][y][z]*f.W[x][y][z]-1.0)*vol;

	  pot_phi_symm += Vf(p, f.T[x][y][z], f.phi[x][y][z])*vol;
#else
	  pot_phi_symm += Vf(p, p.Tconst, f.phi[x][y][z])*vol;
	    
#endif

	  grad_phi_symm += 0.5*((f.phi[x+1][y][z] - f.phi[x][y][z])/p.dx)
	    *((f.phi[x+1][y][z] - f.phi[x][y][z])/p.dx)*vol;
	
	  grad_phi_symm += 0.5*((f.phi[x][y+1][z] - f.phi[x][y][z])/p.dx)
	    *((f.phi[x][y+1][z] - f.phi[x][y][z])/p.dx)*vol;
	
	  grad_phi_symm += 0.5*((f.phi[x][y][(z+1)%p.Lz] 
			  - f.phi[x][y][z])/p.dx)
	    *((f.phi[x][y][(z+1)%p.Lz] 
	       - f.phi[x][y][z])/p.dx)*vol;


	  kin_phi_symm += 0.5*f.pi_future[x][y][z]*f.pi_future[x][y][z]*vol;
	}
	else {
#ifndef SCALAR
	  rest_E_broken += (f.E[x][y][z]/f.W[x][y][z])*vol;

	  kin_fluid_broken += f.kappa[x][y][z]*(f.E[x][y][z]/f.W[x][y][z])
	    *(f.W[x][y][z]*f.W[x][y][z]-1.0)*vol;

	  pot_phi_broken += Vf(p, f.T[x][y][z], f.phi[x][y][z])*vol;
#else
	  pot_phi_broken += Vf(p, p.Tconst, f.phi[x][y][z])*vol;
	  
#endif

	  
	  grad_phi_broken += 0.5*((f.phi[x+1][y][z] - f.phi[x][y][z])/p.dx)
	    *((f.phi[x+1][y][z] - f.phi[x][y][z])/p.dx)*vol;
	
	  grad_phi_broken += 0.5*((f.phi[x][y+1][z] - f.phi[x][y][z])/p.dx)
	    *((f.phi[x][y+1][z] - f.phi[x][y][z])/p.dx)*vol;
	
	  grad_phi_broken += 0.5*((f.phi[x][y][(z+1)%p.Lz] 
			  - f.phi[x][y][z])/p.dx)
	    *((f.phi[x][y][(z+1)%p.Lz] 
	       - f.phi[x][y][z])/p.dx)*vol;

	  

	  kin_phi_broken += 0.5*f.pi_future[x][y][z]*f.pi_future[x][y][z]*vol;
	}
      }
    }
  }
  energies[0] =  rest_E_symm;
  energies[1] =  rest_E_broken;
  energies[2] =  kin_fluid_symm;
  energies[3] =  kin_fluid_broken;
  energies[4] =  kin_phi_symm;
  energies[5] =  kin_phi_broken;
  energies[6] =  grad_phi_symm;
  energies[7] =  grad_phi_broken;
  energies[8] =  pot_phi_symm;
  energies[9] =  pot_phi_broken;  
}

/** Calculates sum the temperature over sites on the core, split into
 *  symmetric and broken phase. Returns as an array of floats, with symmetric
 *  first and broken second.
 *  
 *  Note that this function does not sum over all sites, only the ones on the
 *  current core.
 */
void calculate_T_sum(hydro_fields f, hydro_params p, float *T_sum) {

  float T_sum_symm = 0;
  float T_sum_broken = 0;

  float phi_broken;
  
  int x, y, z;
#ifndef SCALAR
#ifdef BAG
  phi_broken = p.phi_0;
#endif
  
  for(x = 1; x <= p.slicex; x++) {
    for(y = 1; y <= p.slicey; y++) {
      for(z = 0; z < p.Lz; z++) {

#ifndef BAG
	phi_broken =  (p.alpha*f.T[x][y][z]
		       + sqrt((p.alpha*p.Tconst)*(p.alpha*f.T[x][y][z])
			      - 4.0*p.lambda*p.gamma
			      *(f.T[x][y][z]*f.T[x][y][z] - p.T0*p.T0))
		       )/(2.0*p.lambda);
#endif

	if (f.phi[x][y][z] < phi_broken/2){
	  T_sum_symm += f.T[x][y][z];

	}
	else {
	  T_sum_broken += f.T[x][y][z];
	}
      }
    }
  }
#endif // !SCALAR
  T_sum[0] =  T_sum_symm;
  T_sum[1] =  T_sum_broken;
  
}

/** Compute the total pressure. Split this into symmetric and broken phases,
 * symmetric first.
 *
 * NB: This function does _not_ currently sum over all sites.
 */
void calculate_pressure_sum(hydro_fields f, hydro_params p, float *pressure_sum) {

  int x, y, z;

  float vol;

  float press_symm = 0;
  float press_broken = 0;
#ifndef SCALAR
  vol = p.dx*p.dx*p.dx;
  float phi_broken;
#ifdef BAG
  phi_broken = p.phi_0/2.;
#endif // BAG
  
  for(x = 1; x <= p.slicex; x++) {
    for(y = 1; y <= p.slicey; y++) {
      for(z = 0; z < p.Lz; z++) {
#ifndef BAG
	phi_broken =  (p.alpha*f.T[x][y][z]
		       + sqrt((p.alpha*p.Tconst)*(p.alpha*f.T[x][y][z])
			      - 4.0*p.lambda*p.gamma
			      *(f.T[x][y][z]*f.T[x][y][z] - p.T0*p.T0))
		       )/(2.0*p.lambda);
#endif // !BAG

	if (f.phi[x][y][z] < phi_broken/2){
	  press_symm += f.p[x][y][z]*vol;

	}
	else {
	  press_broken += f.p[x][y][z]*vol;
	}
      }
    }
  }
#endif // !SCALAR

  pressure_sum[0] = press_symm;
  pressure_sum[1] = press_broken;

}



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

