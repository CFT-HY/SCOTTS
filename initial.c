/** @file initial.c
 *
 * Initial conditions for hydro evolution.
 *
 * To create a system with bubbles, one would first call
 * initial_blank() *before* nucleating bubbles with nucleate_at()
 * or try_nucleate().
 *
 * Bubbles are generally nucleated as blobs of scalar field with the 
 * Gaussian ansatz. The value of the field at the centre of the blob is
 * is the broken-phase value of the scalar field,
 * \f[
 *  \phi_\text{min} = \frac{\alpha T_N + \sqrt{(\alpha T_N)^2 
 *                                      - 4 \lambda\gamma (T_N^2 - T_0^2)}}
 *                         {2\lambda}
 * \f]
 * 
 * The surface tension is
 * \f[
 *  \sigma = \frac{2\sqrt{2}}{81} \frac{\alpha^3}{\lambda^{5/2}}
 * \f]
 *
 * In the BAG model
 * \f[
 * \phi_\text{min} = \frac{\alpha + \sqrt{\alpha^2 
 *                                      - 4 \lambda\gamma}}
 *                         {2\lambda}
 * \f]
 * \f[
 *  \sigma=\frac{\left(\alpha \left(\alpha + \sqrt{\alpha^2 - 4 \gamma \lambda}
 *                                  \right) -2\gamma \lambda \right)^{3/2}}
 *              {24 \lambda^{5/2}}
 * \f]
 * The radius of the critical bubble is
 * \f[
 *  R_\text{crit} =  \frac{2 \sigma}{V(0,T_N)-V(\phi_b,T_N)}
 * \f]
 * This is scaled by `p.scale` to get the initial conditions for
 * a freshly-nucleated Gaussian blob. This scaling may be necessary
 * to get the bubbles to grow due to e.g. lattice artifacts.
 *
 *
 *
 * These expressions are used in various functions in this file.
 *
 * `phimin` and `Rcritical` can be specified in the parameters file.
 *
 * Contributors:
 * - 2010-2017 David Weir
 * - 2018-     Daniel Cutting
 */
#include "hydro.h"


/** Initialise the system to an empty box with no bubbles.
 *
 * Note that this must be run at the start of the simulation, prior to
 * (for example) nucleate_at(), even if one is just nucleating one
 * bubble.
 *
 * Initialises the scalar field `f.phi` and its conjugate momentum
 * `f.pi_future` to zero everywhere.
 *
 * If we are _not_ running a scalar-only simulation (i.e. `SCALAR` is
 * not defined) we also initialise the temperature to `p.Tconst`,
 * which is the temperature at which we nucleate bubbles; there is no
 * physical input to determine this temperature, it is set 'by hand'.
 *
 * The internal energy of the fluid is then set appropriately.
 *
 * The remaining quantities - the fluid 3-velocity `f.V` and momentum
 * density `f.Z` are set everywhere to zero, and the body-centred
 * gamma factor `f.W` is set to unity.
 *
 * Does not initialise the haloes - you must call halo_field() where
 * appropriate shortly after running this.
 */
void initial_blank(hydro_fields f, hydro_params p) {

  int x, y, z, i;

  for(x = 1; x <= p.slicex; x++) {
    for(y = 1; y <= p.slicey; y++) {
      for(z = 0; z < p.Lz; z++) {
	
	f.phi[x][y][z] = 0.0; 
	
	f.pi_future[x][y][z] = 0.0;

#ifndef SCALAR	


	f.T[x][y][z] = p.Tconst;


	// Note when the fluid is at rest, f.W = 1 and so
	// f.E is just the internal rest-energy i.e. \epsilon.
	
	f.E[x][y][z] = 3.0*p.gdeg*f.T[x][y][z]*f.T[x][y][z]
	  *f.T[x][y][z]*f.T[x][y][z]
	  + Vf(p, f.T[x][y][z], f.phi[x][y][z])
	  - f.T[x][y][z]*VTf(p, f.T[x][y][z], f.phi[x][y][z]);

	
	for(i = 0; i < 3; i++) {
	  f.Z[i][x][y][z] = 0.0;
	  f.V[i][x][y][z] = 0.0;
	}

	f.W[x][y][z] = 1.0;


#endif // SCALAR

      }
    }
  }
}


/** Compute the safe distance required around new bubbles.
 *
 * Note, returns the distance in number of lattice sites, rounded
 * (up or down).
 */
int safe_distance(hydro_fields f, hydro_params p) {

  double safe_lattice_distance = sqrt(2.0*p.R_scaled*p.R_scaled/
				     (p.dx*p.dx));
  
  return (int)round(safe_lattice_distance);
 	
}


/** Check whether a bubble can be nucleated at a given point.
 *
 * Are the conditions suitable for nucleation at `(x0, y0, z0)`?
 * Returns 1 if so, 0 if not.
 *
 * Two conditions are applied:
 *  - The distance between the bubbles must
 *    be at least the 'safe distance' computed by safe_distance().
 *  - There is no scalar radiation or other disturbances
 *    of a greater magnitude than threshold_phi (which is a hard-coded
 *    percentage of the field value at the centre of the Gaussian ansatz).
 */
int can_nucleate(hydro_fields f, hydro_params p, int x0, int y0, int z0) {

  int safedist = safe_distance(f,p);

  double threshold_phi = 0.00005*p.phimin;

  int direct_x, direct_y, direct_z;
  int wrap_x, wrap_y, wrap_z;
  int delta_x, delta_y, delta_z;

  int is_good = 1;
  int x, y, z;

  for(x = 1; x <= p.slicex; x++) {
    for(y = 1; y <= p.slicey; y++) {
      for(z = 0 ; z < p.Lz; z++) {

	// These are distances 'within the lattice'
	direct_x = abs(p.shiftx + x - 1 - x0);
	direct_y = abs(p.shifty + y - 1 - y0);
	direct_z = abs(z - z0);

	// These are distances which 'wrap around' the boundary
	wrap_x = (p.Lx - direct_x);
	wrap_y = (p.Ly - direct_y);
	wrap_z = (p.Lz - direct_z);

	// Check which of direct_ and wrap_ is shorter
	
	if(direct_x < wrap_x)
	  delta_x = direct_x;
	else
	  delta_x = wrap_x;

	if(direct_y < wrap_y)
	  delta_y = direct_y;
	else
	  delta_y = wrap_y;

	if(direct_z < wrap_z)
	  delta_z = direct_z;
	else
	  delta_z = wrap_z;

	// Apply two conditions:
	// 1. The distance between the bubbles must
	//    be at least the 'safe distance' computed by safe_distance.
	// 2. There is no scalar radiation or other disturbances
	//    of a greater magnitude than threshold_phi.
	if(((delta_x*delta_x + delta_y*delta_y + delta_z*delta_z)
	    < safedist*safedist)
	   && (fabs(f.phi[x][y][z]) > threshold_phi)) {

	  fprintf(stderr,
		  "Dead (%d, %d, %d) "
		  "(dx=%d, dy=%d, dz=%d, safe=%d, phi %lf)\n",
		  x0, y0, z0,
		  delta_x, delta_y, delta_z, safedist, f.phi[x][y][z]);

	  is_good = 0;

	  break;

	}
      }
      
      if(!is_good)
	break;
    }
    
    if(!is_good)
      break;
  }

  // Only allow nucleation if *all* the volumes look okay.  
  return reduce_and(is_good, p);
}


/** Nucleate one bubble at a given location.
 *
 * Nucleate a bubble at `(x0, y0, z0)`. Does _not_ check whether a
 * bubble 'can' be nucleated there, see can_nucleate() to test whether
 * this is the case. The Gaussian ansatz is added to the scalar field
 * `f.phi`, and so if one does not check that the area is empty, this
 * can end up overshooting the broken phase.
 */
void nucleate_at(hydro_fields f, hydro_params p, int x0, int y0, int z0) {

  int x, y, z, i;  

  double threshold_phi = 0.00005*p.phimin;
  
  printf0(p, "Nucleating at (%d, %d, %d)\n", x0, y0, z0);

  int direct_x, direct_y, direct_z;
  int wrap_x, wrap_y, wrap_z;
  int delta_x, delta_y, delta_z;

  double phival;

  for(x = 1; x <= p.slicex; x++) {
    for(y = 1; y <= p.slicey; y++) {
      for(z = 0; z < p.Lz; z++) {

	direct_x = abs(p.shiftx + x - 1 - x0);
	direct_y = abs(p.shifty + y - 1 - y0);
	direct_z = abs(z - z0);

	wrap_x = (p.Lx - direct_x);
	wrap_y = (p.Ly - direct_y);
	wrap_z = (p.Lz - direct_z);

	if(direct_x < wrap_x)
	  delta_x = direct_x;
	else
	  delta_x = wrap_x;

	if(direct_y < wrap_y)
	  delta_y = direct_y;
	else
	  delta_y = wrap_y;

	if(direct_z < wrap_z)
	  delta_z = direct_z;
	else
	  delta_z = wrap_z;

	// Gaussian ansatz
	phival = p.phimin*exp(-1.0*p.dx*p.dx*((double)(delta_x*delta_x 
						    + delta_y*delta_y
						    + delta_z*delta_z))
			    /(2.0*p.R_scaled*p.R_scaled));
	
	f.phi[x][y][z] += phival;


#ifndef SCALAR
	// Only bother computing the resulting fluid energy density
	// if the value of phi is 'not small'
	if(phival > threshold_phi) {

	  f.E[x][y][z] = 3.0*p.gdeg*f.T[x][y][z]*f.T[x][y][z]
	    *f.T[x][y][z]*f.T[x][y][z]
	    + Vf(p, f.T[x][y][z], f.phi[x][y][z])
	    - f.T[x][y][z]*VTf(p, f.T[x][y][z], f.phi[x][y][z]);

	}

	// Don't set up the energy density according to the Eqn of state?!
	//	f.E[x][y][z] = 0.0;

#endif // SCALAR
      }
    }
  }
#ifndef SCALAR
  halo_field(f.E,p);
#endif //!SCALAR
  halo_field(f.phi,p);
	
}


/** Full nucleation process including site selection.
 *
 * Attempt exactly one nucleation:
 * - choose a site at random uniformly out of the whole box
 * - check if nucleation is allowed there with can_nucleate()
 * - if so, nucleate with nucleate_at()
 * - if not, do nothing
 */
int try_nucleate(hydro_fields f, hydro_params p) {
  
  printf0(p, "Trying to nucleate a bubble (safe distance = %d)\n", 
	  safe_distance(f, p));

  int try_x = random() % p.Lx;
  int try_y = random() % p.Ly;
  int try_z = random() % p.Lz;


  // NB: number of attempts set to 1, but can always retry
  if(!can_nucleate(f, p, try_x, try_y, try_z)) {

    printf0(p, "Not allowed to nucleate at (%d,%d,%d)!\n",
	    try_x, try_y, try_z);

    return 0;

  } else {
    nucleate_at(f, p, try_x, try_y, try_z);
  }
  

  return 1;
}


/** How many bubbles to nucleate at the given timestep.
 *
 * Provides for random nucleation based a list of timesteps at which
 * nucleation should take place (which may be more than once per
 * timestep towards the end).  If p.nucleation set to NUC_OFF, never
 * nucleates more bubbles.
 *
 * Such a nucleation list can be generated by nucproc.py
 */
int bubbles_at_step(hydro_fields f, hydro_params p, double t_sim, int step) {
  if(p.initial!=INIT_BUBBLE) {
    return 0;
  } else if(p.nucleation == NUC_OFF) {
    return 0;
  } else if((p.nucleation == NUC_LIST) ||
	    (p.nucleation == NUC_FILE) ||
	    (p.nucleation == NUC_FILE_LOC)) {
    
    int j;

    int nuc_count = 0;

    for(j=0; j < p.n_nucsteps; j++) {
      if(p.nucsteps[j] == step) {
	nuc_count++;
	printf0(p, "Parameter file requires nucleation at t=%lf step=%d\n",
		t_sim, step);
      }
    }

    return nuc_count;

  }

  return 0;
}

/** Initialises a spherical overdensity of fluid in the origin of the simulation. 
 *
 * Keeps field in symmetric phase everywhere. Initialises a gaussian overdensity in T.
 */
void fluid_sphere(hydro_fields f, hydro_params p){
#ifdef SCALAR
  printf0(p, "Error - fluid sphere incompatible with SCALAR compiler flag, "
	  "exiting...");
  die(100);
#else
  int x, y, z;
  
  int direct_x, direct_y, direct_z;
  int wrap_x, wrap_y, wrap_z;
  int delta_x, delta_y, delta_z;
  
  double x0 = 0;
  int y0 = 0;
  int z0 = 0;
  
  for(x = 1; x <= p.slicex; x++) {
    for(y = 1; y <= p.slicey; y++) {
      for(z = 0; z < p.Lz; z++) {

	direct_x = abs(p.shiftx + x - 1 - x0);
	direct_y = abs(p.shifty + y - 1 - y0);
	direct_z = abs(z - z0);

	wrap_x = (p.Lx - direct_x);
	wrap_y = (p.Ly - direct_y);
	wrap_z = (p.Lz - direct_z);

	if(direct_x < wrap_x)
	  delta_x = direct_x;
	else
	  delta_x = wrap_x;

	if(direct_y < wrap_y)
	  delta_y = direct_y;
	else
	  delta_y = wrap_y;

	if(direct_z < wrap_z)
	  delta_z = direct_z;
	else
	  delta_z = wrap_z;

	// Symmetric
	f.phi[x][y][z] = 0;

	
	// Gaussian in T
	f.T[x][y][z] = p.Tconst;
	f.T[x][y][z] += (p.T_central - p.Tconst)*exp(-1.0*p.dx*p.dx*
						   ((double)(delta_x*delta_x 
							    + delta_y*delta_y
							    + delta_z*delta_z))
						     /(2.0*p.sphere_radius
						     *p.sphere_radius));
	
	f.E[x][y][z] = (3.0*p.gdeg*f.T[x][y][z]*f.T[x][y][z]
			*f.T[x][y][z]*f.T[x][y][z]
			+ Vf(p, f.T[x][y][z], f.phi[x][y][z])
			- f.T[x][y][z]*VTf(p, f.T[x][y][z], f.phi[x][y][z]));
	
      }
    }
  }
#endif 
}

