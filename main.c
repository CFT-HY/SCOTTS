/* main.c
 * 
 * Execution entrypoint.
 */
#include "hydro.h"


int main(int argc, char *argv[])
{

  hydro_params p;

#ifdef MPI

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &p.rank);
  MPI_Comm_size(MPI_COMM_WORLD, &p.size);

  if(!p.rank)
    fprintf(stderr, "3D hydro port with MPI\n");

#else // MPI

  fprintf(stderr, "3D hydro port (serial)\n");

  p.rank = 0;
  p.size = 1;

#endif // MPI

  if(argc != 2) {
    if(!p.rank)
      fprintf(stderr,"Usage: hydro <parameter file>\n");
    return 100;
  }


#ifdef HAVE_MALLOC_H

  // Inspired by Kari

  mallopt( M_MMAP_MAX, 0 );  /* don't use mmap */
  /* HACK: don't release memory by calling sbrk */
  mallopt( M_TRIM_THRESHOLD, -1 );

  fprintf(stderr, "Disabled sbrk\n");

#endif // HAVE_MALLOC_H


  // Parse params from stdin
  get_parameters(argv[1], &p);

  // How big is the system
  p.N = p.Lx*p.Ly*p.Lz;

  // Calculate terms in potential
  p.alpha = 1.0/sqrt(2.0*p.sigma*
			      pow(p.lcorr, 5.0)/3.0);

  p.gamma = (p.Lheat 
		      + 6.0*p.sigma/p.lcorr) 
    / (6.0*p.sigma*p.lcorr);

  p.lambda = 1.0/(3.0*p.sigma*pow(p.lcorr,3.0));

  // What did we find?
  if(!p.rank) 
    fprintf(stderr,
	    "-- calculated potential terms\n"	\
	    "-- alpha %g, gamma %g, lambda %g\n",	\
	    p.alpha, p.gamma, p.lambda);

  // Make these user-modifiable eventually
  p.Tconst = 0.86;
  p.xxwall = -.5;

  // gdeg = 34.25
  p.a = 34.25*M_PI*M_PI/90.0;

  // Ugly but it's a straight conversion of the fortran
  p.T0 = sqrt(1.0-2.0*p.alpha*p.alpha
		       /9.0/p.gamma/p.lambda);


  // Set up layout for any (or no) parallelism
  layout(&p);

  // Useless
  init_comms_time(&p);

  // time iterate
  int step;
  double t = 0.0;

  // space iterate
  int x;

  double wmax;

  hydro_fields f;

  double initial_energy, current_energy;

  int **nb;
  int **inverse;

  double current_field_energy, current_wmax;

  double cpu_time_used;

  clock_t start, end;


  nb = init_nb(&p);
  inverse = init_inverse(&p);

  if(!p.rank)
    fprintf(stderr, "- Initialised neighbour (and inverse) lookup table.\n");



  // Just runs malloc on all the fields therein
  alloc_fields(&f, p);

  if(!p.rank)
    fprintf(stderr, "- Allocated fields.\n");


  // For safety, set everything to zero
  zero_fields(f, p);

  if(!p.rank)
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

  initial_scalar_bubble(f, p, inverse);
  // initial_3D(f,p,inverse);
  // initial_step(f,p);


  // Communicate everything that neighbours need
  halo_field(f.phi, p);
  halo_field(f.pifull, p);
  halo_field(f.E, p);
  halo_field(f.Z[0], p);
  halo_field(f.Z[1], p);
  halo_field(f.Z[2], p);
  halo_field(f.W, p);


  if(!p.rank)
    fprintf(stderr, "Initial conditions done\n");


  initial_energy = reduce_sum(total_energy(f, nb, p), p);

  if(!p.rank)
    fprintf(stderr, "Initial avg energy per site: %g\n", 
	    initial_energy);


#ifdef SILO

#ifdef MPI
#error SILO and MPI do not play well yet
#endif // MPI

  if(p.silointerval > 0)
    silo_init(p);

#endif // SILO



#ifdef PAPI

  papi_init();

#endif // PAPI



  // Step back for leapfrog initial conds (doubtful we need this)
  // evolve_backstep(f, nb, p);
  // NB: if we later enable this, careful to halo the necessary stuff

  cpu_time_used = 0.0;


#ifdef MPI
  // Mostly so that the timing is remotely fair
  MPI_Barrier(MPI_COMM_WORLD);

#endif // MPI


  start = clock();

  for(step = 0; step < p.steps; step++) {

    if((p.silointerval > 0) && (step % p.silointerval == 0)) {

#ifdef SILO

      write_silo_step(f, p, step);

#endif // SILO

    }

    if((p.interval > 0) && (step % p.interval == 0)) {

      current_energy = reduce_sum(total_energy(f, nb, p), p);
      current_field_energy = reduce_sum(field_energy(f, nb, p), p);
      current_wmax = reduce_max(get_gamma_max(f, p), p);
      
      if(!p.rank)
	fprintf(stderr,"%04d\t%6lf\t%6lf\t%6lf\t%6lf\n",
		step, t,
		current_energy,
		current_field_energy,
		current_wmax);
      
      /*
      if(!p.rank)
	fprintf(stderr, "Energy violation: %lf%%\n",
		100.0*fabs((current_energy-initial_energy)/initial_energy));
      */

    }

    // Do field step
    evolve_field(f, nb, p);
    
    // Calculate EOS
    eq_of_state(f, p);

    // Do the hydro bits
    evolve_hydro(f, nb, p);

   
    // Dump a whole lot of stuff (currently just use Silo for that)    
    /*
    if(step == p.steps - 1)
      for(x=0;x<p.N;x++) {
	if(inverse[x][0] == inverse[x][1]) {
	  fprintf(stdout,"%lf %.10lf %.10lf %.10lf %.10lf %.10lf\n",
		  (inverse[x][0]*p.dx),
		  f.V[0][x],
		  f.U[0][x],
		  f.phi[x],
		  f.T[x],
		  f.Z[0][x]);
	  
	}
      }    
    */

    // Advection of state variables
    advect_E(f, nb, p);
    // Advection of momentum
    advect_Z(f, nb, p);

    // Don't bother with art viscosity, yet
    // artificial_viscosity(f, nb, p);

    // Solve for T
    find_Ta(f, p);
    
    t += p.dt;
    
  } // main loop ends here


  end = clock();


#ifdef PAPI

  papi_finalise();

#endif // PAPI
  

#ifdef MPI

  MPI_Finalize();

#endif // MPI


  cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;

  fprintf(stderr, "CPU time in main loop was %lf, of which %lf was comms\n",
	  cpu_time_used, get_comms_time(&p));


  // Clean up memory
  free_nb(nb, &p);
  free_inverse(inverse, &p);
  free_fields(&f);
  
  return 0;

}





void alloc_fields(hydro_fields *f, hydro_params p) {

  // fields below also set up in initial condition functions
  f->phi = (double *) malloc(p.fieldN*sizeof(double));
  f->pifull = (double *) malloc(p.fieldN*sizeof(double));
  f->T = (double *) malloc(p.fieldN*sizeof(double));
  f->E = (double *) malloc(p.fieldN*sizeof(double));
  f->W = (double *) malloc(p.fieldN*sizeof(double));
  
  // pi gets initialised in backstep
  f->pi = (double *) malloc(p.fieldN*sizeof(double));

  // phiold initialised by evolve_field
  f->phiold = (double *) malloc(p.fieldN*sizeof(double));

  // equation of state (obtained by eq_of_state first time
  // and used in hydro)
  f->kappa = (double *) malloc(p.fieldN*sizeof(double));



  // pressure (obtained by eq_of_state first time
  // and used in hydro)
  f->p = (double *) malloc(p.fieldN*sizeof(double));

  f->V = (double **) malloc(3*sizeof(double *));

  f->V[0] = (double *) malloc(p.fieldN*sizeof(double));
  f->V[1] = (double *) malloc(p.fieldN*sizeof(double));
  f->V[2] = (double *) malloc(p.fieldN*sizeof(double));
  
  f->Z = (double **) malloc(3*sizeof(double *));

  f->Z[0] = (double *) malloc(p.fieldN*sizeof(double));
  f->Z[1] = (double *) malloc(p.fieldN*sizeof(double));
  f->Z[2] = (double *) malloc(p.fieldN*sizeof(double));

  f->U = (double **) malloc(3*sizeof(double *));

  f->U[0] = (double *) malloc(p.fieldN*sizeof(double));
  f->U[1] = (double *) malloc(p.fieldN*sizeof(double));
  f->U[2] = (double *) malloc(p.fieldN*sizeof(double));

  f->F = (double **) malloc(3*sizeof(double *));

  f->F[0] = (double *) malloc(p.fieldN*sizeof(double));
  f->F[1] = (double *) malloc(p.fieldN*sizeof(double));
  f->F[2] = (double *) malloc(p.fieldN*sizeof(double));

  // (calloc considered harmful)
}




void zero_fields(hydro_fields f, hydro_params p) {

  int x;

  // We don't need to do this because create_gaussian_bubble
  // initialises everything, but it keeps valgrind quiet.
  for(x=0;x<p.fieldN;x++) {

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

  free(f->phi);
  free(f->pifull);
  free(f->T);
  free(f->E);
  free(f->W);

  free(f->pi);

  free(f->phiold);

  free(f->kappa);
  free(f->p);

  free(f->V[0]);
  free(f->V[1]);
  free(f->V[2]);

  free(f->V);

  free(f->Z[0]);
  free(f->Z[1]);
  free(f->Z[2]);

  free(f->Z);

  free(f->U[0]);
  free(f->U[1]);
  free(f->U[2]);

  free(f->U);

  free(f->F[0]);
  free(f->F[1]);
  free(f->F[2]);

  free(f->F);



 
}

