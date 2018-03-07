/** @file initial.c
 *
 * Initial conditions for hydro evolution.
 *
 * To create a system with bubbles, one would first call
 * initial_blank() *before* nucleating bubbles with nucleate_at()
 * or do_nucleate().
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
 * The radius of the critical bubble is
 * \f[
 *  R_\text{crit} = - \frac{2 \sigma}{V(\phi,T_N)}
 * \f]
 * This is scaled by `p.scale` to get the initial conditions for
 * a freshly-nucleated Gaussian blob. This scaling may be necessary
 * to get the bubbles to grow due to e.g. lattice artifacts.
 *
 * These expressions are used in various functions in this file.
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
 * `f.pifull` to zero everywhere.
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
	
	f.pifull[x][y][z] = 0.0;

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


  float surface_tension = (2.0*sqrt(2.0)/81.0)*(
						p.alpha*p.alpha*p.alpha
						/(p.lambda*p.lambda
						  *sqrt(p.lambda)));

  float phimin =  (p.alpha*p.Tconst
		   + sqrt((p.alpha*p.Tconst)*(p.alpha*p.Tconst)
			  - 4.0*p.lambda*p.gamma
			  *(p.Tconst*p.Tconst - p.T0*p.T0))
		   )/(2.0*p.lambda);
  
  float R_critical = -2.0*surface_tension/Vf(p, p.Tconst, phimin);
  
  float R_scaled = p.scale*R_critical;

  float safe_lattice_distance = sqrt(2.0*R_scaled*R_scaled/
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

  float phimin =  (p.alpha*p.Tconst
		   + sqrt((p.alpha*p.Tconst)*(p.alpha*p.Tconst)
			  - 4.0*p.lambda*p.gamma
			  *(p.Tconst*p.Tconst - p.T0*p.T0))
		   )/(2.0*p.lambda);

  float threshold_phi = 0.00005*phimin;

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


  float surface_tension = (2.0*sqrt(2.0)/81.0)*(
						p.alpha*p.alpha*p.alpha
						/(p.lambda*p.lambda
						  *sqrt(p.lambda)));

  float phimin =  (p.alpha*p.Tconst
		   + sqrt((p.alpha*p.Tconst)*(p.alpha*p.Tconst)
			  - 4.0*p.lambda*p.gamma
			  *(p.Tconst*p.Tconst - p.T0*p.T0))
		   )/(2.0*p.lambda);
  
  float R_critical = -2.0*surface_tension/Vf(p, p.Tconst, phimin);
  
  float R_scaled = p.scale*R_critical;

  float threshold_phi = 0.00005*phimin;
  
  printf0(p, "Nucleating at (%d, %d, %d)\n", x0, y0, z0);

  int direct_x, direct_y, direct_z;
  int wrap_x, wrap_y, wrap_z;
  int delta_x, delta_y, delta_z;

  float phival;

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
	phival = phimin*exp(-1.0*p.dx*p.dx*((float)(delta_x*delta_x 
						    + delta_y*delta_y
						    + delta_z*delta_z))
			    /(2.0*R_scaled*R_scaled));
	
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
 
}


/** Initialise the system with a shock tube.
 *
 * "Shock tube"-style initial conditions for testing the fluid.
 * No scalar field initialisation.
 */
void initial_3D(hydro_fields f, hydro_params p) {

  
  int x, y, z, i; 

  float surface_tension = 2.0*sqrt(2.0)/81.0*p.alpha*p.alpha*p.alpha
    /(p.lambda*p.lambda*sqrt(p.lambda));

  float phimin =  ( p.alpha*p.Tconst 
		     + sqrt((p.alpha*p.Tconst)*(p.alpha*p.Tconst)
			    - 4*p.lambda*p.gamma
			    *(p.Tconst*p.Tconst - p.T0*p.T0)) )/ (2*p.lambda);


  float R_critical = 2.0*surface_tension/(-1.0*Vf(p, p.Tconst, phimin));

  float R_scaled = p.scale*R_critical;

  printf0(p,
	  "** phimin %g, R_critical %g, R_scaled %g\n",
	  phimin, R_critical, R_scaled);
  
  float Er, El, rE;

  rE = 50.0;

  Er = 1.0/(1.0+rE);
  El = rE*Er;



  for(x=1;x<=p.slicex;x++) {
    for(y=1;y<=p.slicey;y++) {
      for(z=0;z<p.Lz;z++) {

	f.phi[x][y][z] = 0.0; // cstrab*(1.0 + 0.1*drand48());
	
	f.pifull[x][y][z] = 0.0;

#ifndef SCALAR	
	f.T[x][y][z] = 0.0; // p.Tconst;
	
	//  sqrt((x-p.Lx/2)*(x-p.Lx/2)+(y-p.Ly/2)*(y-p.Ly/2)) < 40)
	/*
	if( (x+p.shiftx-1 + y+p.shifty-1) < p.Lx/2
	    || (x+p.shiftx-1 + y+p.shifty-1) > 3*p.Lx/2) 
	  f.E[x][y][z] = El;
	else
	  f.E[x][y][z] = Er;
	*/

	f.E[x][y][z] = 6.0;

	// For debugging purposes
	// fprintf(stderr,"xc = %lf fi = %lf E = %lf\n", xc[x], phi[x],E[x]);
	
	f.Z[0][x][y][z] = 0.0;	
	f.Z[1][x][y][z] = 0.0;
	f.Z[2][x][y][z] = 0.0;
	
	
	f.W[x][y][z] = 1.0;
#endif // SCALAR
      }
    }
  }
}



/** Full nucleation process including site selection.
 *
 * Attempt exactly one nucleation:
 * - choose a site
 * - check if nucleation is allowed there
 * - if so, nucleate
 * - if not, do nothing
 */
int do_nucleate(hydro_fields f, hydro_params p) {
  
  printf0(p, "Trying to nucleate a bubble (safe distance = %d)\n", 
	  safe_distance(f,p));

  int tryx = random()%p.Lx;
  int tryy = random()%p.Ly;
  int tryz = random()%p.Lz;

  int attempts = 0;


  // NB: number of attempts set to 1, but can always retry
  if(!can_nucleate(f,p,tryx,tryy,tryz)) {

    printf0(p, "Not allowed to nucleate at (%d,%d,%d)!\n",
	    tryx, tryy, tryz);

    return 0;

  } else {
    nucleate_at(f, p, tryx, tryy, tryz);
    halo_field(f.phi, p);
  }
  

  return 1;
}


/** How many bubbles to nucleate at the given timestep.
 *
 * Provides for random nucleation based a list of times at which
 * nucleation should take place (which may be more than once per
 * timestep towards the end).  If p.nucleation set to NUC_OFF, never
 * nucleates more bubbles.
 *
 * Such a nucleation list can be generated by nucproc.py
 */
int bubbles_at_step(hydro_fields f, hydro_params p, float t, int step) {

  if(p.nucleation == NUC_OFF) {
    return 0;
  } else if((p.nucleation == NUC_LIST) || (p.nucleation == NUC_FILE)) {

    int j;

    int nuc_count = 0;

    for(j=0; j < p.n_nucsteps; j++) {
      if(p.nucsteps[j] == step) {
	nuc_count++;
	printf0(p, "Parameter file requires nucleation at t=%lf step=%d\n",
		t, step);
      }
    }

    return nuc_count;

  }

  return 0;
}



/** Initialise an invariant scaling profile for new bubbles.
 *
 * Based on my Q-ball code. Not currently in use.
 */
void init_profile(hydro_fields *f, hydro_params *p) {
  float xdummy, phidummy, vdummy;

  printf0(*p, "Loading invariant profile (may be slow)!\n");

  if(access("profile",R_OK) != 0) {
    printf0(*p, "Unable to access profile file, \"profile\"\n");
    die(100);
  }

  FILE *fp = fopen("profile","r");

  int lines = 0;

  while(!feof(fp)
	&& (fscanf(fp,"%f%f%f",&xdummy,&phidummy,&vdummy) == 3)) {
    lines++;
  }


  printf0(*p, "Invariant profile file has %d usable lines\n", lines);


  rewind(fp);

  int j;

  int i = 0;

  int imax;

  float dist;

  p->Linv = lines;
  f->x_inv = (float *) malloc(lines*sizeof(float));
  f->phi_inv = (float *) malloc(lines*sizeof(float));
  f->V_inv = (float *) malloc(lines*sizeof(float));

  while(i < lines) {
    if(fscanf(fp, "%f%f%f", &(f->x_inv[i]), &(f->phi_inv[i]),
	      &(f->V_inv[i])) != 3 && (!(p->rank))) {
      fprintf(stderr, "File inconsistent between first and second reads: "
	      "odd... giving up\n");
      die(100);
    }
    i++;
  }

  fclose(fp);

  printf0(*p, "Done loading invariant profile\n");


  // Not quite: we want ...???

}
