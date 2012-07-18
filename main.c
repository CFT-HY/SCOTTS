/* main.c
 * 
 * Execution entrypoint.
 */
#include "hydro.h"


void alloc_fields(hydro_fields *f, hydro_params p) {
  // location labels set up in initial condition functions
  f->xe = (double *) malloc(p.N*sizeof(double));
  f->xc = (double *) malloc(p.N*sizeof(double));

  // fields below also set up in initial condition functions
  f->phi = (double *) malloc(p.N*sizeof(double));
  f->pifull = (double *) malloc(p.N*sizeof(double));
  f->T = (double *) malloc(p.N*sizeof(double));
  f->E = (double *) malloc(p.N*sizeof(double));
  f->W = (double *) malloc(p.N*sizeof(double));
  
  // pi gets initialised in backstep
  f->pi = (double *) malloc(p.N*sizeof(double));

  // phiold initialised by evolve_field
  f->phiold = (double *) malloc(p.N*sizeof(double));

  // equation of state (obtained by eq_of_state first time
  // and used in hydro)
  f->kappa = (double *) malloc(p.N*sizeof(double));



  // pressure (obtained by eq_of_state first time
  // and used in hydro)
  f->p = (double *) malloc(p.N*sizeof(double));

  f->V = (double **) malloc(3*sizeof(double *));

  f->V[0] = (double *) malloc(p.N*sizeof(double));
  f->V[1] = (double *) malloc(p.N*sizeof(double));
  f->V[2] = (double *) malloc(p.N*sizeof(double));
  
  f->Z = (double **) malloc(3*sizeof(double *));

  f->Z[0] = (double *) malloc(p.N*sizeof(double));
  f->Z[1] = (double *) malloc(p.N*sizeof(double));
  f->Z[2] = (double *) malloc(p.N*sizeof(double));

  f->U = (double **) malloc(3*sizeof(double *));

  f->U[0] = (double *) malloc(p.N*sizeof(double));
  f->U[1] = (double *) malloc(p.N*sizeof(double));
  f->U[2] = (double *) malloc(p.N*sizeof(double));

  f->F = (double **) malloc(3*sizeof(double *));

  f->F[0] = (double *) malloc(p.N*sizeof(double));
  f->F[1] = (double *) malloc(p.N*sizeof(double));
  f->F[2] = (double *) malloc(p.N*sizeof(double));

  // (calloc considered harmful)
}


void dump(double *field, hydro_params p) {
  int x;

  fprintf(stderr,"%g", field[0]);
  for(x=1;x<p.N;x++) {
    fprintf(stderr,", %g", field[x]);
  }
  fprintf(stderr,"\n");
}

void zero_fields(hydro_fields f, hydro_params p) {

  int x;

  // We don't need to do this because create_gaussian_bubble
  // initialises everything, but it keeps valgrind quiet.
  for(x=0;x<p.N;x++) {
    f.xe[x] = 0.0000;
    f.xc[x] = 0.0000;

    f.phi[x] = 0.0000;
    f.pifull[x] = 0.0000;
    f.T[x] = 0.0000;
    f.E[x] = 0.0000;
    f.Z[0][x] = 0.0000;
    f.Z[1][x] = 0.0000;
    f.Z[2][x] = 0.0000;
    f.U[0][x] = 0.0000;
    f.U[1][x] = 0.0000;
    f.U[2][x] = 0.0000;
    f.V[0][x] = 0.0000;
    f.V[1][x] = 0.0000;
    f.V[2][x] = 0.0000;
    f.W[x] = 0.0000;

    f.pi[x] = 0.0000;
    f.phiold[x] = 0.0000;
    f.kappa[x] = 0.0000;
    f.p[x] = 0.0000;
  }

}

void free_fields(hydro_fields *f) {
  free(f->xe);
  free(f->xc);
  free(f->phi);
  free(f->pifull);
  free(f->phiold);
  free(f->T);
  free(f->E);
  free(f->Z[0]);
  free(f->Z[1]);
  free(f->Z[2]);
  free(f->U[0]);
  free(f->U[1]);
  free(f->U[2]);
  free(f->V[0]);
  free(f->V[1]);
  free(f->V[2]);
  free(f->W);
  free(f->kappa);
  free(f->p);
  free(f->pi);

  free(f->F[0]);
  free(f->F[1]);
  free(f->F[2]);
  free(f->F);

  free(f->U);
  free(f->V);
 
}

int iix(int x, int y, int z, hydro_params p) {
  return ((z+p.Lz)%p.Lz)*p.Ly*p.Lx + ((y+p.Ly)%p.Ly)*p.Lx + ((x+p.Lx)%p.Lx);
}

int **init_nb(hydro_params p) {
  int x, y, z;

  int **nb;

  // set up nb lookup table
  nb = (int **) malloc(p.N*sizeof(int *));
  
  for(x=0; x<p.Lx; x++) {
    for(y=0; y<p.Ly; y++) {
      for(z=0; z<p.Lz; z++) {
    
	nb[iix(x, y, z, p)] = (int *)malloc(14*sizeof(int));

	nb[iix(x, y, z, p)][0] = iix(x+1, y, z, p);
	nb[iix(x, y, z, p)][1] = iix(x-1, y, z, p);
        nb[iix(x, y, z, p)][2] = iix(x, y+1, z, p);
	nb[iix(x, y, z, p)][3] = iix(x, y-1, z, p);
	nb[iix(x, y, z, p)][4] = iix(x, y, z+1, p);
	nb[iix(x, y, z, p)][5] = iix(x, y, z-1, p);

	/*
	 * Composed directions improve performance as soon as the system
	 * ceases to fit inside cache. At small volumes (and hence
	 * lower dimensionalities) it is up to 20% slower, but
	 * we are interested in good performance on large volumes:
	 * not likely that each node will be able to fit everything
	 * in cache.
	 */
      
	nb[iix(x, y, z, p)][DIR_02] = iix(x+1, y+1, z, p);
	nb[iix(x, y, z, p)][DIR_04] = iix(x+1, y, z+1, p);
        nb[iix(x, y, z, p)][DIR_24] = iix(x, y+1, z+1, p);
	nb[iix(x, y, z, p)][DIR_13] = iix(x-1, y-1, z, p);
	nb[iix(x, y, z, p)][DIR_15] = iix(x-1, y, z-1, p);
	nb[iix(x, y, z, p)][DIR_35] = iix(x, y-1, z-1, p);

	nb[iix(x, y, z, p)][DIR_024] = iix(x+1, y+1, z+1, p);
	nb[iix(x, y, z, p)][DIR_135] = iix(x-1, y-1, z-1, p);
      }
    }
  }

  return nb;
}

void free_nb(int **nb, hydro_params p) {
  int x;

  for(x=0; x<p.N; x++) {
    free(nb[x]);
  }
  
  free(nb);
}


int main(int argc, char *argv[])
{
  fprintf(stderr,"3D hydro port\n");

  if(argc != 2) {
    fprintf(stderr,"Usage: hydro <parameter file>\n");
    return 100;
  }

  mallopt( M_MMAP_MAX, 0 );  /* don't use mmap */
  /* HACK: don't release memory by calling sbrk */
  mallopt( M_TRIM_THRESHOLD, -1 );


  // Parse params from stdin
  hydro_params p = get_parameters(argv[1]);

  

  p.N = p.Lx*p.Ly*p.Lz;

  // Calculate terms in potential
  p.alpha = 1.0/sqrt(2.0*p.sigma*
			      pow(p.lcorr, 5.0)/3.0);

  p.gamma = (p.Lheat 
		      + 6.0*p.sigma/p.lcorr) 
    / (6.0*p.sigma*p.lcorr);

  p.lambda = 1.0/(3.0*p.sigma*pow(p.lcorr,3.0));

  // What did we find?
  fprintf(stderr,
	  "-- calculated potential terms\n" \
	  "-- alpha %g, gamma %g, lambda %g\n", \
	  p.alpha, p.gamma, p.lambda);

  // Make these user-modifiable eventually
  p.Tconst = 0.86;
  p.xxwall = -.5;

  // gdeg = 34.25
  p.a = 34.25*M_PI*M_PI/90.0;

  // Ugly but it's a straight conversion of the fortran
  p.T0 = sqrt(1.0-2.0*p.alpha*p.alpha
		       /9.0/p.gamma/p.lambda);


  // time iterate
  int step;
  double t = 0.0;

  // space iterate
  int x;

  double wmax;

  hydro_fields f;

  double initial_energy, current_energy;

  // Just runs malloc on all the fields therein
  alloc_fields(&f, p);

  fprintf(stderr, "- Allocated fields.\n");

  // For safety, set everything to zero
  zero_fields(f, p);

  fprintf(stderr, "- Zeroed fields.\n");

  // initial conditions
  /*  if(p.initial == INIT_BUBBLE) {
    create_gaussian_bubble(f, p);
  } else if(p.initial == INIT_SHOCK_TUBE) {
    create_shock_tube(f, p);
  } else {
    fprintf(stderr, "Unknown initial condition parameter!\n");
    exit(100);
  }
  */

  initial_scalar_bubble(f,p);
  //      initial_3D(f,p);
    // initial_step(f,p);

  fprintf(stderr, "Initial conditions done\n");

  int **nb;

  nb = init_nb(p);

  fprintf(stderr, "- Initialised neighbour lookup table.\n");

  initial_energy = total_energy(f, nb, p);

  fprintf(stderr, "Initial avg energy per site: %g\n", 
	  initial_energy);

  // Output headers
  // printf("Step\ttime\tenergy\twallps\n");

  // Output of some field at every time step
  //  FILE *phi_fh = fopen("field.dat","w");
#ifdef SILO
  if(p.silointerval > 0)
    silo_init(p);
#endif // SILO

#ifdef PAPI
  papi_init();
#endif // PAPI

  // Step back for leapfrog initial conds
  //  evolve_backstep(f, nb, p);

  // Record how many lattice sites in output file, to keep self contained
  //  fwrite(&p.N, sizeof(int), 1, phi_fh);


  for(step = 0; step < p.steps; step++) {

    if((p.silointerval > 0) && (step % p.silointerval == 0)) {
#ifdef SILO
      write_silo_step(f, p, step);
#endif // SILO
    }

    if((p.interval > 0) && (step % p.interval == 0)) {


      // Write time step, x coords and velocity field
      //      fwrite(&t, sizeof(double), 1, phi_fh);
      //      fwrite(f.E, sizeof(double), p.N, phi_fh);
      //      fwrite(f.p, sizeof(double), p.N, phi_fh);

      current_energy = total_energy(f, nb, p);

      
      fprintf(stderr,"%04d\t%6lf\t%6lf\t%6lf\t%6lf\n",
	     step, t,
	     current_energy, 
	     field_energy(f, nb, p), wmax);
      

      fflush(stdout);

      fprintf(stderr, "Energy violation: %lf%%\n",
	      100.0*fabs((current_energy-initial_energy)/initial_energy));
    }

    // Do field step
    evolve_field(f, nb, p);

    
    
    // Calculate EOS
    eq_of_state(f, p);



    // Do the hydro bits
    evolve_hydro(f, nb, p);

    if(step == p.steps - 1)
      for(x=0;x<p.Lx;x++) {
	fprintf(stdout,"%lf %.10lf %.10lf %.10lf %.10lf %.10lf %.10lf %.10lf\n",
		(x*p.dx),
		f.V[0][iix(x,x,0,p)],
		f.V[0][iix(x,x+1,0,p)],
		f.V[0][iix(x,x-1,0,p)],
		f.U[0][iix(x,x,0,p)], // f.Ux[nb[iix(x,x,0,p)][2]],
		f.phi[iix(x,x,0,p)], // f.Ux[nb[iix(x,x,0,p)][4]],
		f.T[iix(x,x,0,p)], // f.Ux[nb[nb[iix(x,x,0,p)][2]][4]],
		f.Z[0][iix(x,x,0,p)]);
      }


    // Advection of state variables
    advect_E(f, nb, p);
    // Advection of momentum
    advect_Z(f, nb, p);


    //    artificial_viscosity(f, nb, p);



    // Solve for T
    find_Ta(f, p);


    
    t += p.dt;

    wmax = get_gamma_max(f, p);


    
  } // main loop ends here



#ifdef PAPI
  papi_finalise();
#endif // PAPI
  

  //  fclose(phi_fh);

  // Clean up memory
  free_nb(nb, p);
  free_fields(&f);
  
  return 0;
}
