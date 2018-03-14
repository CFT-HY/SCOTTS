/** @file notused.c
 *
 * Functions that are not currently used and have been taken
 * out of the other files during refactoring by David & Daniel.
 */


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
	
	f.pi_future[x][y][z] = 0.0;

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
	f.pi[x][y][z] = f.pi_future[x][y][z] - 0.5*p.dt
	  *(f.phi[x+1][y][z] + f.phi[x][y+1][z] 
	    + f.phi[x][y][((z+1)%p.Lz)] 
	    - 6.0*f.phi[x][y][z] 
	    + f.phi[x-1][y][z] + f.phi[x][y-1][z]
	    + f.phi[x][y][((z-1+p.Lz)%p.Lz)])/(p.dx*p.dx)
	  + 0.125*p.dt*p.dt
	  *(f.pi_future[x+1][y][z] - 2.0*f.pi_future[x][y][z]
	    + f.pi_future[x-1][y][z])/(p.dx*p.dx)
#ifndef SCALAR
	  + 0.5*p.dt*Vdf(p, f.T[x][y][z], f.phi[x][y][z] 
			 - 0.25*p.dt*f.pi_future[x][y][z]);
#else
	  + 0.5*p.dt*Vdf(p, p.Tconst, f.phi[x][y][z] 
			 - 0.25*p.dt*f.pi_future[x][y][z]);
#endif
      }
    }
  }

  halo_field(f.pi, p);
}
