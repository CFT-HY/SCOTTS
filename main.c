/* main.c
 * 
 * Execution entrypoint.
 */
#include "hydro.h"


int main(int argc, char *argv[]) {

  hydro_params p;

#ifdef MPI

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &p.rank);
  MPI_Comm_size(MPI_COMM_WORLD, &p.size);

  printf0(p, "3D hydro port with MPI\n");

#else // not MPI

  printf0(p, "3D hydro port (serial)\n");

  p.rank = 0;
  p.size = 1;

#endif // MPI

  // Note that this is when main.c was rebuilt, not the whole thing!
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

  // Scale factor (fixed to 1 for now, no expansion)
  p.a = 1.0;

  // Calculate potential terms, no longer used (potential directly supplied)
  /*
  p.alpha = 1.0/sqrt(2.0*p.sigma*pow(p.lcorr, 5.0)/3.0);
  p.gamma = (p.Lheat + 6.0*p.sigma/p.lcorr)/(6.0*p.sigma*p.lcorr);
  p.lambda = 1.0/(3.0*p.sigma*pow(p.lcorr,3.0));

  printf0(p,
	  "-- calculated potential terms\n"
	  "-- alpha %g, gamma %g, lambda %g\n",
	  p.alpha, p.gamma, p.lambda);

  p.T0 = sqrt(1.0-2.0*p.alpha*p.alpha/9.0/p.gamma/p.lambda);

  printf0(p,
	  "-- T0 %g\n", p.T0);
  */


  
  // Degrees of freedom, still hardcoded
  p.gdeg = p.gstar*M_PI*M_PI/90.0;


  // Set up layout for any (or no) parallelism
  layout(&p);

  // Start a counter for communications time
  init_comms_time(&p);

  // Initialised bubble count
  int bcount = 0;

  // On this timestep, do we continue to nucleate?
  int still_nucleate = 1;

  // How many bubbles to nucleate on a given timestep
  int howmany;

  // time iterates: count steps and t separately
  int step;
  float t = 0.0;

  // Struct that stores the fields
  hydro_fields f;

  // Storage of measurements of average stress-energy tensor (not used)
  //  float cpts[TENSOR_CPTS];

  float initial_energy, current_energy;
  float initial_field_energy, current_field_energy;
  float current_kinetic, current_gradient_energy, current_rest;
  float current_avgpress;
  float current_veltot;
  float gwen;

  // Timing counters
  float cpu_time_used;

  clock_t start, end;

  // What step does the for loop start on? For checkpoint restarts...
  int step_start = 0;


  int i;

  // Just runs malloc on all the fields therein (see alloc.c)
  alloc_fields(&f, p);

  // Managed to get this far, so we probably have enough memory
  printf0(p, "- Allocated fields.\n");

  // Restore checkpoint
  if(usable_checkpoint(f, p)) {
    printf0(p, "Found a usable checkpoint file\n");

    step_start = load_checkpoint(f, p);

    printf0(p, "Checkpoint load done.\n");
    printf0(p, "WARNING: Bubble count reset after restart!\n");

  } else {
    // For safety, set everything to zero
    zero_fields(f, p);

    printf0(p, "- Zeroed fields.\n");

    if(p.initial==INIT_PS){
#ifndef SCALAR
      // If initial is INIT_PS then initialise velocity power spec.
      // Then no bubbles nucleated.

	initial_blank(f,p);
	init_ps(f, p, f.Z);
	norm_power(f, p, f.Z);

	/*
	  memcpy(f.V[0][0][0], f.U[0][0][0], (p.slicex+2)*(p.slicey+2)
	  *(p.Lz)*sizeof(float));
	  memcpy(f.V[1][0][0], f.U[1][0][0], (p.slicex+2)*(p.slicey+2)
	  *(p.Lz)*sizeof(float));
	  memcpy(f.V[2][0][0], f.U[2][0][0], (p.slicex+2)*(p.slicey+2)
	  *(p.Lz)*sizeof(float));
	  */
	// init_ps(f, p, f.V);
	//    init_ps(f, p, f.Z);
#else
	
	printf0(p,"INIT_PS initial conditions invalid with SCALAR compiler flag.",
		"Exiting... \n");
	die(100);
	
#endif // !SCALAR	

	
    } else if(p.initial==INIT_BUBBLE){
      //Bubble initial conditions:
      if(p.bubbles > 1){
	// Instead, start off with an empty box
	initial_blank(f, p);
      
	start = clock();
	printf0(p, "Nucleating first bubble\n");
	nucleate_at(f,p,0,0,0);
	halo_field(f.phi,p);
      
	
	end = clock();
	if(!p.rank)
	  fprintf(stderr,"Nucleation attempt took %lf\n",
		  ((float) (end - start)) / CLOCKS_PER_SEC);
	bcount+=1;
	
	// p.bubbles is how many bubbles to make at the start (usually 1)
	for(step=1;step<p.bubbles;step++) {
	  start = clock();
	  still_nucleate = try_nucleate(f, p);
	  end = clock();
	  if(!p.rank)
	    fprintf(stderr,"Nucleation attempt took %lf\n",
		    ((float) (end - start)) / CLOCKS_PER_SEC);
	  
	  bcount += still_nucleate;
      
	  // Each is an independent attempt, do not disable!
	  //      if(!still_nucleate)
	  //	break;
	  
	}
      } else if(p.bubbles == 1){
	// One bubble only (not normally used)
	initial_blank(f, p);
	
	printf0(p, "Nucleating just one bubble\n");
	nucleate_at(f,p,0,0,0);
	halo_field(f.phi,p);
      } else{
	// Empty system
	initial_blank(f, p);
      }
    }else if(p.initial==INIT_SHOCK_TUBE){
      fprintf(stderr,"Sorry, shocktube is not implemented! Exiting... \n");
      die(100);
    }else{
      fprintf(stderr,"Invalid initial condition option! Exiting... \n");
      die(100);
    }
	
  
      
    // (don't) turn off nucleation after initial stage
    //  still_nucleate = 0;



    // Communicate everything that neighbours need
    halo_field(f.phi, p);
    halo_field(f.pi_future, p);
#ifndef SCALAR
    halo_field(f.E, p);
    halo_field(f.Z[0], p);
    halo_field(f.Z[1], p);
    halo_field(f.Z[2], p);
    halo_field(f.W, p);
#endif // SCALAR

    if(!p.rank)
      fprintf(stderr, "Initial conditions done\n");
  }

  initial_energy = reduce_sum(total_energy(f, p), p);
  initial_field_energy = reduce_sum(field_energy(f, p), p);

  printf0(p, "Initial avg energy per site: %g\n", 
	  initial_energy/((float)p.N));


#ifdef SILO

  // Are we using silo for visualisation?
  if(p.silointerval > 0 || p.silosliceinterval>0 )
    silo_init(p);

#endif // SILO



#ifdef PAPI

  // Are we using PAPI performance counters?
  papi_init();

#endif // PAPI

  // Reset cpu time
  cpu_time_used = 0.0;

  // Step back for leapfrog initial conds (doubtful we need this)
  // evolve_backstep(f, p);
  // NB: if we later enable this, careful to halo the necessary stuff



#ifdef MPI

  // Mostly so that the timing is remotely fair
  MPI_Barrier(MPI_COMM_WORLD);

#endif // MPI

  // For making sound files of the pressure
  /*
  FILE *press;


  if(!p.rank) {
    press = fopen("pressure.dat","a");
  }
  */

  printf0(p, "Starting main loop.\n");

  // Wall time measurement
  start = clock();



  for(step = step_start; step < p.steps; step++) {

    /*
    if(!p.rank) {
      fprintf(press, "%d %g\n", step, f.p[1][1][1]);
    }
    */



    // For reproducibility, reseed with timestep
    // so bubbles nucleate in same places (locations given by random())
    srandom(p.seed + step);

    // How many bubbles do we (try to) nucleate this timestep?
    // (Always 0 if initial condition not set to "bubble")
    howmany = bubbles_at_step(f, p, t, step);

    i = 0;

    while(i < howmany) {

      still_nucleate = try_nucleate(f, p);

      bcount += still_nucleate;
      i++;
    }


    // Checkpoint if necessary
    if((p.checkpointinterval > 0) && (step % p.checkpointinterval == 0) \
       && (step != step_start)) {
      int wall_start = (int)time(NULL);
      checkpoint(f, p, step);
      printf0(p, "Checkpointing took %d seconds of walltime\n",
	      ((int)time(NULL))-wall_start);
    }


#ifdef SILO
    // Write visualisation stuff if necessary
    if((p.silointerval > 0) && (step % p.silointerval == 0)) {
      write_silo_step(f, p, step);
    }

    if((p.silosliceinterval > 0) && (step % p.silosliceinterval == 0)) {
      write_silo_slice_step(f, p, step);
    }
#endif // SILO




#ifdef FFT
    if((p.fftinterval > 0) && (step % p.fftinterval == 0)) {

      histo_field(f.phi, p, step);


      if(p.uetcstart >= 0 && step >= p.uetcstart) {
	if(step == p.uetcstart) {
	  init_uetc(f, p);
	}

	// Calculate UETCs (not in use)
	fft_uetc(f, p, step);
      }      

      // Power spectrum of scalar field
      fft_field(f, p, f.phi, step);

#ifndef SCALAR
      // Velocity power spectrum
      fft_vel(f, p, step, f.V);
#endif //!SCALAR

      // Gravitational wave power spectrum (returns GW energy)
      gwen = fft_tensor(f, p, step, 0.0);

      // Average size of stress-energy tensor components (not used)
      /*
      didj(cpts, f, p);
      printf0(p, "cpts: %g %g %g %g %g %g\n",
	      cpts[CPT_11],
	      cpts[CPT_21],
	      cpts[CPT_31],
	      cpts[CPT_22],
	      cpts[CPT_32],
	      cpts[CPT_33]);
      */
    }
#endif // FFT


    // Measurements
    if((p.interval > 0) && (step % p.interval == 0)) {

      current_energy = reduce_sum(total_energy(f, p), p);
      current_kinetic = reduce_sum(kinetic_energy(f, p), p);
      current_rest = reduce_sum(rest_energy(f, p), p);
      current_avgpress = reduce_sum(avg_pressure(f, p), p);
      current_field_energy = reduce_sum(field_energy(f, p), p);
      current_gradient_energy = reduce_sum(gradient_energy(f, p), p);
      current_veltot = reduce_sum(get_veltot(f, p), p)
	/((float)(p.Lx*p.Ly*p.Lz));
      
      if(!p.rank) {
	printf("%04d\t%6lf\t%6lf\t%6lf\t%6lf\t%6lf\t%6lf\t%6lf\t%d\t%6lf\t%6lf\n",
	       step,
	       t,
	       current_energy,
	       current_kinetic,
	       current_field_energy,
	       current_gradient_energy,
	       current_veltot,
	       gwen,
	       bcount,
	       current_rest,
	       current_avgpress);
      }

      // Statement of energy violation (not shown; better to use KE)
      /*
      fprintf(stderr, "Energy violation: %.10lf%%, %lf%%\n",
              100.0*fabs((current_energy-initial_energy)
			 /initial_energy),
              100.0*fabs((current_field_energy-initial_field_energy)
			 /initial_field_energy));
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

    // Advection of state variables
    advect_E(f, p);
    // Advection of momentum
    advect_Z(f, p);

    // Don't bother with art viscosity, yet
    //    artificial_viscosity(f, nb, p);

    // Solve for T
    find_Ta(f, p);
    
    t += p.dt;


    // On the last step, do some extra GW power spectra FFTs
#ifdef FFT
    if(step == p.steps - 1) {

#ifndef SCALAR
      fft_vel(f, p, step, f.V);
#endif // !SCALAR

      fft_tensor(f,p,step,current_energy);
    }
#endif // FFT
    
  } // main loop ends here


  if(p.checkpointinterval > 0) {
    printf0(p, "Final checkpoint...\n");
    checkpoint(f, p, step);
  }


  current_energy = reduce_sum(total_energy(f, p), p);
  current_kinetic = reduce_sum(kinetic_energy(f, p), p);
  current_rest = reduce_sum(rest_energy(f, p), p);
  current_avgpress = reduce_sum(avg_pressure(f, p), p);
  current_field_energy = reduce_sum(field_energy(f, p), p);
  current_gradient_energy = reduce_sum(gradient_energy(f, p), p);
  current_veltot = reduce_sum(get_veltot(f, p), p)
    /((float)(p.Lx*p.Ly*p.Lz));
      
  if(!p.rank) {
    printf("%04d\t%6lf\t%6lf\t%6lf\t%6lf\t%6lf\t%6lf\t%6lf\t%d\t%6lf\t%6lf\n",
	   step,
	   t,
	   current_energy,
	   current_kinetic,
	   current_field_energy,
	   current_gradient_energy,
	   current_veltot,
	   gwen,
	   bcount,
	   current_rest,
	   current_avgpress);
  }

  
  // End time, for walltime calculation
  end = clock();


#ifdef FFT

#ifndef SCALAR
  fft_vel(f,p,step,f.V);
#endif // !SCALAR

  fft_tensor(f,p,step,current_energy);
    
#endif // FFT



#ifdef PAPI

  papi_finalise();

#endif // PAPI
  

#ifdef MPI

  MPI_Finalize();

#endif // MPI


  // Time spent running
  cpu_time_used = ((float) (end - start)) / CLOCKS_PER_SEC;

  printf0(p, "On master node, CPU time in main loop was %lfs,\n"
	  "of which %lfs was comms\n",
	  cpu_time_used, get_comms_time(&p));

  // Close pressure measurement file (not used)
  /*
  if(!p.rank) {
    fclose(press);
  }
  */

  // Clean up memory
  free_fields(&f, p);
  
  return 0;

}
