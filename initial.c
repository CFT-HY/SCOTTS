/* initial.c
 *
 * Initial conditions for hydro evolution.
 *
 * Also, nucleation.
 */
#include "hydro.h"


/* void initial_scalar_bubble(hydro_fields f, hydro_params p)
 *
 * Create a single isolated scalar bubble at the
 * centre of the box.
 */
void initial_scalar_bubble(hydro_fields f, hydro_params p) {

  int x, y, z, i;

  double sigmlo = 2.0*sqrt(2.0)/81.0*p.alpha*p.alpha*p.alpha
    /(p.lambda*p.lambda*sqrt(p.lambda));

  printf0(p,
	  "** Initial conditions magic:\n"
	  "** sigmlo %g\n", sigmlo);
  
  double phimin =  ( p.alpha*p.Tconst
		     + sqrt((p.alpha*p.Tconst)*(p.alpha*p.Tconst)
			    - 4*p.lambda*p.gamma
			    *(p.Tconst*p.Tconst - p.T0*p.T0)) )/(2*p.lambda);

  double cstrab = 1.0*phimin;

  double Rlapab = 2.0*sigmlo/(-1.0*Vf(p, p.Tconst, phimin));

  double Rtenab = p.scale*Rlapab;

  printf0(p,
	  "** phimin %g, cstrab %g, Rtenab %g\n",
	  phimin, cstrab, Rtenab);


  for(x=1;x<=p.slicex;x++) {
    for(y=1;y<=p.slicey;y++) {
      for(z=0;z<p.Lz;z++) {
    
	f.phi[x][y][z] = cstrab*exp(-1.0*p.dx*p.dx
				    *(((double)(p.shiftx+x-1-p.Lx/2))
				      *((double)(p.shiftx+x-1-p.Lx/2))
				      + ((double)(p.shifty+y-1-p.Ly/2))
				      *((double)(p.shifty+y-1-p.Ly/2))
				      + ((double)(z-p.Lz/2))
				      *((double)(z-p.Lz/2)))
				    /2.0/(Rtenab*Rtenab) );
	
	f.pifull[x][y][z] = 0.0;
	
	f.T[x][y][z] = p.Tconst;
	
	f.E[x][y][z] = 3.0*p.gdeg*f.T[x][y][z]*f.T[x][y][z]
	  *f.T[x][y][z]*f.T[x][y][z]
	  + Vf(p, f.T[x][y][z], f.phi[x][y][z])
	  - f.T[x][y][z]*VTf(p, f.T[x][y][z], f.phi[x][y][z]);
	
	f.Z[0][x][y][z] = 0.0;
	f.V[0][x][y][z] = 0.0;
	f.W[x][y][z] = 1.0;
      }
    }
  }
}




/* void initial_blank(hydro_fields f, hydro_params p)
 *
 * Empty box with no bubbles.
 */
void initial_blank(hydro_fields f, hydro_params p) {

  int x, y, z, i;

  double sigmlo = 2.0*sqrt(2.0)/81.0*p.alpha*p.alpha*p.alpha
    /(p.lambda*p.lambda*sqrt(p.lambda));

  printf0(p,
	  "** Initial conditions magic:\n"
	  "** sigmlo %g\n", sigmlo);

  double phimin =  ( p.alpha*p.Tconst
		     + sqrt((p.alpha*p.Tconst)*(p.alpha*p.Tconst)
			    - 4*p.lambda*p.gamma
			    *(p.Tconst*p.Tconst - p.T0*p.T0)) )/(2*p.lambda);
  
  double cstrab = 1.0*phimin;

  double Rlapab = 2.0*sigmlo/(-1.0*Vf(p, p.Tconst, phimin));

  double Rtenab = p.scale*Rlapab;

  printf0(p,
	  "** phimin %g, cstrab %g, Rtenab %g\n",
	  phimin, cstrab, Rtenab);


  for(x=1;x<=p.slicex;x++) {
    for(y=1;y<=p.slicey;y++) {
      for(z=0;z<p.Lz;z++) {
	
	f.phi[x][y][z] = 0.0; 
	
	f.pifull[x][y][z] = 0.0;
	
	f.T[x][y][z] = p.Tconst;
	
	f.E[x][y][z] = 3.0*p.gdeg*f.T[x][y][z]*f.T[x][y][z]
	  *f.T[x][y][z]*f.T[x][y][z]
	  + Vf(p, f.T[x][y][z], f.phi[x][y][z])
	  - f.T[x][y][z]*VTf(p, f.T[x][y][z], f.phi[x][y][z]);
	
	f.Z[0][x][y][z] = 0.0;
	f.V[0][x][y][z] = 0.0;
	f.W[x][y][z] = 1.0;
      }
    }
  }
}



/* int safe_distance(hydro_fields f, hydro_params p)
 *
 * Returns the safe distance required around a newly
 * nucleated bubble.
 */
int safe_distance(hydro_fields f, hydro_params p) {


  double sigmlo = 2.0*sqrt(2.0)/81.0*p.alpha*p.alpha*p.alpha
    /(p.lambda*p.lambda*sqrt(p.lambda));

  double phimin =  ( p.alpha*p.Tconst
		     + sqrt((p.alpha*p.Tconst)*(p.alpha*p.Tconst)
			    - 4*p.lambda*p.gamma
			    *(p.Tconst*p.Tconst - p.T0*p.T0)) )/ (2*p.lambda);

  double cstrab = 1.0*phimin;

  double Rlapab = 2.0*sigmlo/(-1.0*Vf(p, p.Tconst, phimin));
  
  double Rtenab = p.scale*Rlapab;

  double sigma = 1.0/(sqrt(p.dx*p.dx/2.0/(Rtenab*Rtenab) ));
  
  return (int)(round(sigma));
 	
}


/* int can_nucleate(hydro_fields f, hydro_params p, int x0, int y0, int z0)
 *
 * Are the conditions suitable for nucleation at (x0, y0, z0)?
 * Returns 1 if so.
 */
int can_nucleate(hydro_fields f, hydro_params p, int x0, int y0, int z0) {

  int safedist = safe_distance(f,p);
  double threshphi = 0.0001;

  int shortx, shorty, shortz;
  int longx, longy, longz;
  int dx, dy, dz;

  int is_good = 1;
  int x, y, z;

  for(x=1;x<=p.slicex;x++) {
    for(y=1;y<=p.slicey;y++) {
      for(z=0;z<p.Lz;z++) {

	shortx = abs(p.shiftx+x-1-x0);
	shorty = abs(p.shifty+y-1-y0);
	shortz = abs(z - z0);

	longx = (p.Lx - shortx);
	longy = (p.Ly - shorty);
	longz = (p.Lz - shortz);

	if(shortx < longx)
	  dx = shortx;
	else
	  dx = longx;

	if(shorty < longy)
	  dy = shorty;
	else
	  dy = longy;

	if(shortz < longz)
	  dz = shortz;
	else
	  dz = longz;

	if(((dx*dx + dy*dy + dz*dz) < safedist*safedist)
	   && (fabs(f.phi[x][y][z]) > threshphi)) {
	  fprintf(stderr,
		  "Dead (%d,%d,%d) (dx=%d, dy=%d, dz=%d, safe=%d, phi %lf)\n",
		  x0, y0, z0, dx, dy, dz, safedist, f.phi[x][y][z]);
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


/* nucleate_at(hydro_fields f, hydro_params p, int x0, int y0, int z0)
 *
 * Nucleate a bubble at (x0, y0, z0). Equivalent to initial_scalar_bubble
 * at a given point.
 */
void nucleate_at(hydro_fields f, hydro_params p, int x0, int y0, int z0) {

  int x, y, z, i;  
  
  double sigmlo = 2.0*sqrt(2.0)/81.0*p.alpha*p.alpha*p.alpha
    /(p.lambda*p.lambda*sqrt(p.lambda));



  double phimin =  ( p.alpha*p.Tconst
		     + sqrt((p.alpha*p.Tconst)*(p.alpha*p.Tconst)
			    - 4*p.lambda*p.gamma
			    *(p.Tconst*p.Tconst - p.T0*p.T0)) )/ (2*p.lambda);
  
  double cstrab = 1.0*phimin;
  
  double Rlapab = 2.0*sigmlo/(-1.0*Vf(p, p.Tconst, phimin));
  
  double Rtenab = p.scale*Rlapab;

  printf0(p, "Nucleating at (%d,%d,%d)\n",x0,y0,z0);

  int shortx, shorty, shortz;
  int longx, longy, longz;
  int dx, dy, dz;

  double phival;

  for(x=1;x<=p.slicex;x++) {
    for(y=1;y<=p.slicey;y++) {
      for(z=0;z<p.Lz;z++) {

	shortx = abs(p.shiftx+x-1-x0);
	shorty = abs(p.shifty+y-1-y0);
	shortz = abs(z - z0);

	longx = (p.Lx - shortx);
	longy = (p.Ly - shorty);
	longz = (p.Lz - shortz);

	if(shortx < longx)
	  dx = shortx;
	else
	  dx = longx;

	if(shorty < longy)
	  dy = shorty;
	else
	  dy = longy;

	if(shortz < longz)
	  dz = shortz;
	else
	  dz = longz;


	phival = cstrab*exp(-1.0*p.dx*p.dx*((double)((dx*dx) 
						     + (double)(dy*dy)
						     + (double)(dz*dz)))
			    /2.0/(Rtenab*Rtenab) );
	
	f.phi[x][y][z] += phival;

	if(phival > 1e-8) {

	  f.E[x][y][z] = 3.0*p.gdeg*f.T[x][y][z]*f.T[x][y][z]
	    *f.T[x][y][z]*f.T[x][y][z]
	    + Vf(p, f.T[x][y][z], f.phi[x][y][z])
	    - f.T[x][y][z]*VTf(p, f.T[x][y][z], f.phi[x][y][z]);

	}

	// Don't set up the energy density according to the Eqn of state?
	//	f.E[x][y][z] = 0.0;
      }
    }
  }
  
  //  fprintf(stderr,"Done\n");
}


/* void initial_3D(hydro_fields f, hydro_params p)
 *
 * "Shock tube"-style initial conditions for testing the fluid.
 * No scalar field initialisation.
 */
void initial_3D(hydro_fields f, hydro_params p) {
 
  int x, y, z, i; 

  double sigmlo = 2.0*sqrt(2.0)/81.0*p.alpha*p.alpha*p.alpha
    /(p.lambda*p.lambda*sqrt(p.lambda));

  double phimin =  ( p.alpha*p.Tconst 
		     + sqrt((p.alpha*p.Tconst)*(p.alpha*p.Tconst)
			    - 4*p.lambda*p.gamma
			    *(p.Tconst*p.Tconst - p.T0*p.T0)) )/ (2*p.lambda);

  double cstrab = 1.0*phimin;

  double Rlapab = 2.0*sigmlo/(-1.0*Vf(p, p.Tconst, phimin));

  double Rtenab = p.scale*Rlapab;

  printf0(p,
	  "** phimin %g, cstrab %g, Rtenab %g\n",
	  phimin, cstrab, Rtenab);
  
  double Er, El, rE;

  rE = 50.0;

  Er = 1.0/(1.0+rE);
  El = rE*Er;



  for(x=1;x<=p.slicex;x++) {
    for(y=1;y<=p.slicey;y++) {
      for(z=0;z<p.Lz;z++) {

	f.phi[x][y][z] = 0.0; // cstrab*(1.0 + 0.1*drand48());
	
	f.pifull[x][y][z] = 0.0;
	
	f.T[x][y][z] = 0.0; // p.Tconst;
	
	//  sqrt((x-p.Lx/2)*(x-p.Lx/2)+(y-p.Ly/2)*(y-p.Ly/2)) < 40)
	if( (x+p.shiftx-1 + y+p.shifty-1) < p.Lx/2
	    || (x+p.shiftx-1 + y+p.shifty-1) > 3*p.Lx/2) 
	  f.E[x][y][z] = El;
	else
	  f.E[x][y][z] = Er;
	
	
	// For debugging purposes
	// fprintf(stderr,"xc = %lf fi = %lf E = %lf\n", xc[x], phi[x],E[x]);
	
	f.Z[0][x][y][z] = 0.0;	
	f.Z[1][x][y][z] = 0.0;
	f.Z[2][x][y][z] = 0.0;
	
	
	f.W[x][y][z] = 1.0;
      }
    }
  }
}



/* int do_nucleate(hydro_fields f, hydro_params p) 
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


/* int should_nucleate(hydro_fields f, hydro_params p, double t, int step)
 *
 * Should an attempt to nucleate be carried out at the current step?
 *
 * Provides for either random nucleation based on an exponentially
 * growing nucleation rate (NUC_EXP) or else nucleation based
 * on a list of times at which nucleation should take place
 * (which may be more than once per timestep towards the end).
 *
 * Such a nucleation list can be generated by nucproc.py
 */
int should_nucleate(hydro_fields f, hydro_params p, double t, int step) {

  if(p.nucleation == NUC_EXP) {

    double end = ((double)p.steps)*p.dt;
    double prob = exp(p.beta*(t-end));

    srand48(1);

    double randin = drand48();


    if(randin < prob) {
      printf0(p, "Recommending nucleation at t=%lf, prob=%lf, random=%lf\n",
	      t, prob, randin);

    //    fprintf(stderr,
    //	    "Rank %d recommending nucleation at t %lf (probability=%lf)\n",
    //	    p.rank, t, prob);
      return 1;
    }
  } else {

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



/* void init_profile(hydro_fields *f, hydro_params *p)
 *
 * Initialise an invariant scaling profile for new bubbles.
 * Based on my Q-ball code. Not currently in use.
 */
void init_profile(hydro_fields *f, hydro_params *p) {
  double xdummy, phidummy, vdummy;

  printf0(*p, "Loading invariant profile (may be slow)!\n");

  if(access("profile",R_OK) != 0) {
    printf0(*p, "Unable to access profile file, \"profile\"\n");
    die(100);
  }

  FILE *fp = fopen("profile","r");

  int lines = 0;

  while(!feof(fp)
	&& (fscanf(fp,"%lf%lf%lf",&xdummy,&phidummy,&vdummy) == 3)) {
    lines++;
  }


  printf0(*p, "Invariant profile file has %d usable lines\n", lines);


  rewind(fp);

  int j;

  int i = 0;

  int imax;

  double dist;

  p->Linv = lines;
  f->x_inv = (double *) malloc(lines*sizeof(double));
  f->phi_inv = (double *) malloc(lines*sizeof(double));
  f->V_inv = (double *) malloc(lines*sizeof(double));

  while(i < lines) {
    if(fscanf(fp, "%lf%lf%lf", &(f->x_inv[i]), &(f->phi_inv[i]),
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
