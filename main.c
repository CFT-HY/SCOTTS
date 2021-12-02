/* main.c
 *
 * Execution entrypoint.
 */
#include "hydro.h"

int main(int argc, char* argv[])
{

    hydro_params p;

#ifdef USE_MPI

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &p.rank);
    MPI_Comm_size(MPI_COMM_WORLD, &p.size);

    printf0(p, "3D hydro port with MPI\n");

#else // not USE_MPI

    printf0(p, "3D hydro port (serial)\n");

    p.rank = 0;
    p.size = 1;

#endif // MPI

    // Note that this is when main.c was rebuilt, not the whole thing!
    printf0(p, "Built: %s %s\n", __DATE__, __TIME__);

    if (argc != 2) {
        printf0(p, "Usage: hydro <parameter file>\n");
        return 100;
    }

#ifdef HAVE_MALLOC_H

    // Inspired by Kari

    mallopt(M_MMAP_MAX, 0); /* don't use mmap */
    /* HACK: don't release memory by calling sbrk */
    mallopt(M_TRIM_THRESHOLD, -1);

    printf0(p, "Disabled sbrk\n");

#endif // HAVE_MALLOC_H

    // Parse parameters from the filename specified on the command line
    get_parameters(argv[1], &p);

    // Seed - for initps, each core throws independent random numbers
    if (p.initial == INIT_PS) {

        int stride = (p.size > p.Lx) ? ((int)(p.size / p.Lx)) : 1;

        srandom(p.seed + (int)(p.rank / stride));
        srand48(p.seed + (int)(p.rank / stride));

    }
    // Seed - make sure everyone gets the same one (if necessary)
    else {
        srandom(p.seed);
        srand48(p.seed);
    }

    // How big is the system
    p.N = p.Lx * p.Ly * p.Lz;

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
    p.gdeg = p.gstar * M_PI * M_PI / 90.0;

    // Set up layout for any (or no) parallelism
    layout(&p);

    // Start a counter for communications time
    init_comms_time(&p);

    // Initialised bubble count
    int bcount = 0;

    // On this timestep, do we continue to nucleate?
    int still_nucleate = 1;

    // Track failed number of attempts at nucleating a bubble
    int current_attempts = 0;

    // How many bubbles to nucleate on a given timestep
    int howmany;

    // time iterates: count steps and time separately
    int step;
    float t_sim = 0.0;

    // Struct that stores the fields
    hydro_fields f;

#ifdef FFT
    // Struct that stores fields in momentum space and helper fields for FFTs.
    fft_fields fft_f;
#endif // FFT

    // Counter to specify which order to perform advection
    int adv_order = 0;

    // Storage of measurements of average stress-energy tensor (not used)
    //  float cpts[TENSOR_CPTS];

    float initial_energy, initial_field_energy;
    float gwen = 0;

    // Timing counters
    float cpu_time_used;

    clock_t start, end;

    // What step does the for loop start on? For checkpoint restarts...
    int step_start = 0;

    int i;

    // Just runs malloc on all the fields therein (see alloc.c)
    alloc_fields(&f, p);

    printf0(p, "- Allocated fields.\n");

#ifdef FFT
    fft_init(p, &fft_f);

    // Managed to get this far, so we probably have enough memory
    printf0(p, "- Initialised fft.\n");
#endif // FFT

    // Restore checkpoint
    if (usable_checkpoint(f, p)) {
        printf0(p, "Found a usable checkpoint file\n");

        step_start = load_checkpoint(f, p);

        printf0(p, "Checkpoint load done.\n");
        printf0(p, "WARNING: Bubble count reset after restart!\n");

    } else {
        // For safety, set everything to zero
        zero_fields(f, p);

        printf0(p, "- Zeroed fields.\n");

        if (p.initial == INIT_PS) {
#if defined(FFT) && (BAG) && !defined(SCALAR)
            // If initial is INIT_PS then initialise velocity power spec.
            // Then no bubbles nucleated.
            start = clock();
            initial_blank(f, p);
            init_ps(f, p, f.U);

            end = clock();

            if (!p.rank) {
                fprintf(stderr, "Init ps initialisation took %lf\n",
                        ((float)(end - start)) / CLOCKS_PER_SEC);
            }

#else

            printf0(p, "INIT_PS initial conditions require FFT and BAG flags,"
                       " invalid with SCALAR flag."
                       " Exiting... \n");
            die(100);

#endif // FFT && BAG && !SCALAR

        } else if (p.initial == INIT_BUBBLE) {
            if (p.nucleation == NUC_FILE_LOC) {
                initial_blank(f, p);
                if (p.bubbles > 0) {
                    printf0(
                        p,
                        "Using nucleation location file: \n"
                        "Ignoring initial bubbles parameter (bubbles %d). \n",
                        p.bubbles);
                }
            }
            // Bubble initial conditions:
            else if (p.bubbles > 1) {
                // Instead, start off with an empty box
                initial_blank(f, p);

                start = clock();
                printf0(p, "Nucleating first bubble\n");
                nucleate_at(f, p, 0, 0, 0);

                end = clock();
                if (!p.rank)
                    fprintf(stderr, "Nucleation attempt took %lf\n",
                            ((float)(end - start)) / CLOCKS_PER_SEC);
                bcount += 1;

                // p.bubbles is how many bubbles to make at the start (usually 1)
                for (step = 1; step < p.bubbles; step++) {
                    start = clock();
                    current_attempts = 0;
                    while (current_attempts < p.maxattempts) {
                        still_nucleate = try_nucleate(f, p);
                        if (still_nucleate) {
                            break;
                        } else {
                            current_attempts++;
                        }
                    }
                    end = clock();
                    if (!p.rank)
                        fprintf(stderr, "Nucleation attempt took %lf\n",
                                ((float)(end - start)) / CLOCKS_PER_SEC);

                    bcount += still_nucleate;

                    // Each is an independent attempt, do not disable!
                    //      if(!still_nucleate)
                    //	break;
                }
            } else if (p.bubbles == 1) {
                // One bubble only (not normally used)
                initial_blank(f, p);

                printf0(p, "Nucleating just one bubble\n");
                nucleate_at(f, p, 0, 0, 0);
                bcount += 1;
            } else {
                // Empty system
                initial_blank(f, p);
            }
        } else if (p.initial == INIT_SHOCK_TUBE) {
            fprintf(stderr,
                    "Sorry, shocktube is not implemented! Exiting... \n");
            die(100);
        } else if (p.initial == INIT_FLUID_SPHERE) {
            initial_blank(f, p);
            // Create a sphere of fluid at origin of the simulation box.
            fluid_sphere(f, p);
        } else {
            fprintf(stderr, "Invalid initial condition option! Exiting... \n");
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

        if (!p.rank)
            fprintf(stderr, "Initial conditions done\n");
    }

    initial_energy = reduce_sum(total_energy(f, p), p);
    initial_field_energy = reduce_sum(field_energy(f, p), p);

    printf0(p, "Initial avg energy per site: %g\n",
            initial_energy / ((float)p.N));

#ifdef SILO

    // Are we using silo for visualisation?
    if (p.silointerval > 0 || p.silosliceinterval > 0)
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

#ifdef USE_MPI

    // Mostly so that the timing is remotely fair
    MPI_Barrier(MPI_COMM_WORLD);

#endif // USE_MPI

    // For making sound files of the pressure
    /*
    FILE *press;


    if(!p.rank) {
      press = fopen("pressure.dat","a");
    }
    */

    // Spit out to stdout global headers
    write_global_headers(f, p);

    printf0(p, "Starting main loop.\n");

    // Wall time measurement
    start = clock();

    // Next index in list of bubble locations for NUC_LOC_FILE option.
    int bub_loc_ind = 0;

    for (step = step_start; step < p.steps; step++) {

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
        howmany = bubbles_at_step(f, p, t_sim, step);

        i = 0;

        while (i < howmany) {
            if (p.nucleation == NUC_FILE_LOC) {
                printf0(p,
                        "Nucleating a bubble on step %d at (%d, %d, %d)"
                        " without any checks.\n",
                        step, p.nuclocs[bub_loc_ind][0],
                        p.nuclocs[bub_loc_ind][1], p.nuclocs[bub_loc_ind][2]);

                nucleate_at(f, p, p.nuclocs[bub_loc_ind][0],
                            p.nuclocs[bub_loc_ind][1],
                            p.nuclocs[bub_loc_ind][2]);
                bub_loc_ind++;
                bcount++;
                i++;
            } else {
                current_attempts = 0;
                while (current_attempts < p.maxattempts) {
                    still_nucleate = try_nucleate(f, p);
                    if (still_nucleate) {
                        break;
                    } else {
                        current_attempts++;
                    }
                }
                bcount += still_nucleate;
            }
        }

        // Checkpoint if necessary
        if ((p.checkpointinterval > 0) && (step % p.checkpointinterval == 0)
            && (step != step_start)) {
            int wall_start = (int)time(NULL);
            checkpoint(f, p, step);
            printf0(p, "Checkpointing took %d seconds of walltime\n",
                    ((int)time(NULL)) - wall_start);
        }

#ifdef SILO
        // Write visualisation stuff if necessary
        if ((p.silointerval > 0) && (step % p.silointerval == 0)) {
            write_silo_step(f, p, step);
        }

        if ((p.silosliceinterval > 0) && (step % p.silosliceinterval == 0)) {
            write_silo_slice_step(f, p, step, p.siloslicecoord);
        }
#endif // SILO

#ifdef FFT

        // Initialise uetcs if it is time.
        if (p.uetcstart >= 0 && step == p.uetcstart) {
            init_uetc(f, p, fft_f);
        }

        if ((p.fftinterval > 0) && (step % p.fftinterval == 0)) {
            if (p.uetcscalar == 1) {
                long long N_broken = reduce_sum(get_N_broken(f, p), p);
                if (N_broken / p.N > p.uetcbrokenthresh) {
                    p.uetcstart = step;
                    printf0(p,
                            "Broken phase fraction threshold exceeded"
                            " starting uetcs on step %d\n",
                            step);
                    init_uetc(f, p, fft_f);
                    p.uetcscalar = 0;
                }
            }
            gwen = output_ps_uetcs(f, p, fft_f, step);
        }
#endif // FFT

        // Measurements
        if ((p.interval > 0) && (step % p.interval == 0)) {

	  write_globals(f, p, gwen, bcount, t_sim, step);

	  // Statement of energy violation (not shown; better to use KE)
	  /*
	    fprintf(stderr, "Energy violation: %.10lf%%, %lf%%\n",
	    100.0*fabs((current_energy-initial_energy)
	    /initial_energy),
	    100.0*fabs((current_field_energy-initial_field_energy)
	    /initial_field_energy));
	  */
        }

	// Store lorentz factor for dW/dt terms
	evolve_hydro_storeWold(f, p);
	
	// Do field step.
	// If initial is INIT_PS, field is set in broken phase everywhere
	// and won't evolve if compiled with BAG potential.
	if(p.initial != INIT_PS)
	  evolve_field(f, p);
	eq_of_state(f, p);
	// Do the hydro bits
	if(p.initial != INIT_PS)
	  evolve_hydro_fieldfluid(f, p);
	evolve_hydro_pressureacceleration(f, p);
	evolve_hydro_velocities(f, p);
	if ((p.kq > 0) || (p.kl > 0))
	  evolve_hydro_artviscosity(f, p);
	evolve_hydro_pressurework(f, p);
	// Advection of state variables
	//uncomment to advect half step and then reverse order.
	// In that case don't perform advect_E and advect_Z. 
	//advect_halfsteps(f, p); 
	advect_E(f, p, adv_order);
	advect_Z(f, p, adv_order);
	adv_order +=1;

	// Evolve dW/dt terms.
	evolve_hydro_boostfactor(f, p);
	//dump_max_min(f, p);
	// Solve for T.
	if(p.initial != INIT_PS)
	  find_Ta(f, p);

    
	//printf0(p,"Evolved hydro \n");
	//dump_max_min(f, p);

	// Evolve metric perturbations
	if (step >= p.metricstart) 
	  evolve_uij(f, p);


	t_sim += p.dt;

    } // main loop ends here


    if(p.checkpointinterval > 0) {
      printf0(p, "Final checkpoint...\n");
      checkpoint(f, p, step);
    }

    // Do ffts one last time for final timestep:
#ifdef FFT

    gwen = output_ps_uetcs(f, p, fft_f, step);
    fft_finalise(p, &fft_f);

#endif // FFT

    // Write globals one last time.
    write_globals(f, p, gwen, bcount, t_sim, step);

#ifdef PAPI

    papi_finalise();

#endif // PAPI

#ifdef USE_MPI

    MPI_Finalize();

#endif // USE_MPI

    // End time, for walltime calculation
    end = clock();

    // Time spent running
    cpu_time_used = ((float)(end - start)) / CLOCKS_PER_SEC;

    printf0(p,
            "On master node, CPU time in main loop was %lfs,\n"
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
