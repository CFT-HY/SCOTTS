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

  printf0(p, "3D hydro port with MPI\n");

#else // MPI

  printf0(p, "3D hydro port (serial)\n");

  p.rank = 0;
  p.size = 1;

#endif // MPI

  printf0(p, "Built: %s %s\n", __DATE__, __TIME__);

  if(argc != 2) {
    printf0(p, "Usage: hydro <parameter file>\n");
    return 100;
  }


#ifdef HAVE_MALLOC_H

  // Inspired by Kari

  mallopt( M_MMAP_MAX, 0 );  /* don't use mmap */
  /* HACK: don't release memory by calling sbrk */
  mallopt( M_TRIM_THRESHOLD, -1 );

  printf0(p, "Disabled sbrk\n");


#endif // HAVE_MALLOC_H


  // Parse parameters from the filename specified on the command line
  get_parameters(argv[1], &p);

  // Seed - make sure everyone gets the same one (if necessary)
  srandom(p.seed);

  // How big is the system
  p.N = p.Lx*p.Ly*p.Lz;


  //  p.scale = 0.3;

  // Calculate terms in potential
  p.alpha = 1.0/sqrt(2.0*p.sigma*pow(p.lcorr, 5.0)/3.0);

  p.gamma = (p.Lheat + 6.0*p.sigma/p.lcorr)/(6.0*p.sigma*p.lcorr);

  p.lambda = 1.0/(3.0*p.sigma*pow(p.lcorr,3.0));


  // What did we find?
  printf0(p,
	  "-- calculated potential terms\n"
	  "-- alpha %g, gamma %g, lambda %g\n",
	  p.alpha, p.gamma, p.lambda);


  // Make these user-modifiable eventually
  p.Tconst = 0.86;

  // gdeg = 34.25
  p.a = 34.25*M_PI*M_PI/90.0;

  // Ugly but it's a straight conversion of the fortran
  p.T0 = sqrt(1.0-2.0*p.alpha*p.alpha/9.0/p.gamma/p.lambda);

  printf0(p,
	  "-- T0 %g\n", p.T0);

  // Set up layout for any (or no) parallelism
  layout(&p);


  // Start a counter for communications time
  init_comms_time(&p);


  int bcount = 0;

  // time iterates: count steps and t
  int step;
  double t = 0.0;

  // space iterate
  int x;

  double wmax;

  hydro_fields f;

  double initial_energy, current_energy;

  double current_field_energy, current_wmax;

  double cpu_time_used;

  clock_t start, end;

  int still_nucleate = 1;

  int step_start = 0;

  double gwen;

  // Just runs malloc on all the fields therein (see alloc.c)
  alloc_fields(&f, p);

  // Managed to get this far, so we probably have enough memory
  printf0(p, "- Allocated fields.\n");

  if(p.checkpointinterval > 0  && usable_checkpoint(f, p)) {
    printf0(p, "Found a usable checkpoint file\n");

    step_start = load_checkpoint(f, p);

    printf0(p, "Checkpoint load done.\n");
    printf0(p, "WARNING: Bubble count reset after restart!\n");

  } else {
    // For safety, set everything to zero
    zero_fields(f, p);

    printf0(p, "- Zeroed fields.\n");

    // initial conditions -- one or more bubbles, or blank and then nucleate?
    /*  if(p.initial == INIT_BUBBLE) {
	create_gaussian_bubble(f, p);
	} else if(p.initial == INIT_SHOCK_TUBE) {
	create_shock_tube(f, p);
	} else {
	fprintf(stderr, "Unknown initial condition parameter!\n");
	die(100);
	}
    */

    initial_blank(f, p);

    // You get one bubble for free
    //  initial_scalar_bubble(f, p);

    for(step=0;step<p.bubbles;step++) {
      start = clock();
      still_nucleate = do_nucleate(f, p);
      end = clock();
      if(!p.rank)
	fprintf(stderr,"Nucleation attempt took %lf\n",
		((double) (end - start)) / CLOCKS_PER_SEC);
      
      bcount += still_nucleate;
      
      // Each is an independent attempt, do not disable!
      //      if(!still_nucleate)
      //	break;

    }

    //  still_nucleate = 0;
    // initial_3D(f,p);
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
  }

  initial_energy = reduce_sum(total_energy(f, p), p);

  printf0(p, "Initial avg energy per site: %g\n", 
	  initial_energy/((double)p.N));



#ifdef SILO

  if(p.silointerval > 0)
    silo_init(p);

#endif // SILO



#ifdef PAPI

  papi_init();

#endif // PAPI



  // Step back for leapfrog initial conds (doubtful we need this)
  // evolve_backstep(f, p);
  // NB: if we later enable this, careful to halo the necessary stuff

  cpu_time_used = 0.0;


#ifdef MPI

  // Mostly so that the timing is remotely fair
  MPI_Barrier(MPI_COMM_WORLD);

#endif // MPI


  printf0(p, "Starting main loop.\n");


  start = clock();

  double t00;

  for(step = step_start; step < p.steps; step++) {


    // Disabled: always try!
    //     if(still_nucleate && should_nucleate(f, p, t, step)) {

    if(should_nucleate(f, p, t, step)) {
      still_nucleate = do_nucleate(f, p);

      bcount += still_nucleate;
    }


    if((p.checkpointinterval > 0) && (step % p.checkpointinterval == 0) \
       && (step != step_start)) {
      int wall_start = (int)time(NULL);
      checkpoint(f, p, step);
      printf0(p, "Checkpointing took %d seconds of walltime\n", ((int)time(NULL))-wall_start);
    }


    if((p.silointerval > 0) && (step % p.silointerval == 0)) {

#ifdef SILO

      write_silo_step(f, p, step);

#endif // SILO

    }

    if((p.interval > 0) && (step % p.interval == 0)) {

      histo_field(f.phi, p, step);



#ifdef FFT
      
      fft_vel(f, p, step);
      gwen = fft_tensor(f, p, step, current_energy);
    
#endif // FFT

      current_energy = reduce_sum(total_energy(f, p), p);
      current_field_energy = reduce_sum(field_energy(f, p), p);
      current_wmax = reduce_max(get_gamma_max(f, p), p);
      
      if(!p.rank) {
	printf("%04d\t%6lf\t%6lf\t%6lf\t%6lf\t%6lf\t%d\n",
	       step,
	       t,
	       current_energy,
	       current_field_energy,
	       current_wmax,
	       gwen,
	       bcount);
      }
      
      /*
	printf0(p, "Energy violation: %lf%%\n",
		100.0*fabs((current_energy-initial_energy)/initial_energy));
      */

    }

    // Do field step
    evolve_field(f, p);
    
    // Calculate EOS
    eq_of_state(f, p);

    // Do the hydro bits
    evolve_hydro(f, p);

    // Evolve metric perturbations
    evolve_uij(f, p);

    // Dump a whole lot of stuff (currently just use Silo for that)    
    if(step == p.steps - 1) {

#ifdef FFT
      
      fft_tensor(f,p,step,current_energy);
    
#endif // FFT
    

    /*
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

    }

    // Advection of state variables
    advect_E(f, p);
    // Advection of momentum
    advect_Z(f, p);

    // Don't bother with art viscosity, yet
    // artificial_viscosity(f, nb, p);

    // Solve for T
    find_Ta(f, p);
    
    t += p.dt;
    
  } // main loop ends here

  if(p.checkpointinterval > 0) {
    printf0(p, "Final checkpoint...\n");
    checkpoint(f, p, step);
  }
  
  current_energy = reduce_sum(total_energy(f, p), p);
  current_field_energy = reduce_sum(field_energy(f, p), p);
  current_wmax = reduce_max(get_gamma_max(f, p), p);
  
  //  printf0(p, "Final state:\n");
  if(!p.rank) {
    printf("%04d\t%6lf\t%6lf\t%6lf\t%6lf\t%d\n",
	   step,
	   t,
	   current_energy,
	   current_field_energy,
	   current_wmax,
	   bcount);
  }

  end = clock();



#ifdef FFT
      
  fft_tensor(f,p,step,current_energy);
    
#endif // FFT



#ifdef PAPI

  papi_finalise();

#endif // PAPI
  

#ifdef MPI

  MPI_Finalize();

#endif // MPI


  cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;

  printf0(p, "On master node, CPU time in main loop was %lfs,\n"
	  "of which %lfs was comms\n",
	  cpu_time_used, get_comms_time(&p));
  

  // Clean up memory
  free_fields(&f, p);
  
  return 0;

}
