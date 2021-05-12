/** @file output.c
 *
 * Simple things, some of which we could in principle calculate offline, but
 * that are nice to know (energy-related quantites are in energy.c).
 */
#include "hydro.h"



/** Output to std out the global headers.
 */
void write_global_headers(hydro_fields f, hydro_params p){
  if(!p.rank) {
    printf("step,time,rest E sym,rest E broken,kin fluid sym,kin fluid broken,"
	   "kin phi sym,kin phi broken,grad phi sym,grad phi broken,"
	   "pot phi sym,pot phi broken,pressure sym,pressure broken,"
	   "T sym,T broken,V^2 tot sym,V^2 tot broken,"
	   "GW Energy,Nb,s_max,"
	   "gamma_max,curl J,div J,N broken,N links\n");
  }
}

/** Write out all global quantities to std out. Takes gwen calculated
 * from fft_tensor, number of bubbles at current time, and the
 * simulation time and timestep.
 *
 * Summary of global outputs in order:
 *
 * step: timestep at output
 * time: time at output
 * rest E sym: sum over rest energy in symmetric phase
 * rest E broken: sum over rest energy in broken phase
 * kin fluid sym: sum over kinetic energy in fluid in symmetric phase
 * kin fluid broken: sum over kinetic energy in fluid in broken phase
 * kin phi sym: sum over kinetic energy in scalar field in symmetric phase
 * kin phi broken: sum over kinetic energy in scalar field in broken phase
 * grad phi sym: sum over gradient energy in scalar field in symmetric phase.
 * grad phi broken: sum over gradient energy in scalar field in broken phase.
 * pot phi sym: sum over potential energy in scalar field in symmetric phase. 
 * pot phi broken: sum over potential energy in scalar field in broken phase.
 * pressure sym: sum over pressure in symmetric phase.
 * pressure broken: sum over pressure in broken phase.
 * T sym: sum over temperature at each site in symmetric phase 
 *        (divide by (p.N - N_broken) to get average temperature in symmetric phase)
 * T broken: sum over temperature at each site in broken phase
 *           (divide by (N_broken) to get average temperature in broken phase)
 * V^2 tot sym: sum over V^2 in symmetric phase. Can be used to find vrms
 * V^2 tot broken: sum over V^2 in broken phase. Can be used to find vrms
 * GW Energy: Average gravitational wave energy. (In units with G = 1).
 * Nb: Number of bubbles in the simulation.
 * s_max: Maximum s parameter in Crank-Nicholson update (see evolve_scalar)
 * gamma_max: Maximum gamma factor of fluid velocity in simulation.
 * curl J: sum over (curl J)^2 for all sites.
 * div J: sum over (div J)^2 for all sites.
 * N broken: Number of sites in broken phase.
 * N links: Number of faces over which the scalar field changes between
 *          symmetric and broken phase.
 */
void write_globals(hydro_fields f, hydro_params p, float gwen, int bcount,
		    float t_sim, int step){

  float s_max;
  float gamma_max;
  float curlJ_tot;
  float divJ_tot;
  long long N_broken;
  long long N_links;

  float energies[10];
  float T_sum[2];
  float pressure_sum[2];
  float Vsq_sum[2];

  calculate_energies(f, p, energies);
  calculate_T_sum(f, p, T_sum);
  calculate_pressure_sum(f, p, pressure_sum);
  calculate_Vsq_sum(f, p, Vsq_sum);

  reduce_sum_array(energies, p, 10);
  reduce_sum_array(T_sum, p, 2);
  reduce_sum_array(pressure_sum, p, 2);
  reduce_sum_array(Vsq_sum, p, 2);

  
  N_broken = reduce_sum(get_N_broken(f,p), p);
  N_links = reduce_sum(get_broken_links(f,p), p);

  s_max = reduce_max(get_s_max(f, p), p);
  gamma_max = reduce_max(get_gamma_max(f, p), p);
  curlJ_tot =  reduce_sum(get_curlJ_tot(f,p),p);
  divJ_tot = reduce_sum(get_divJ_tot(f,p), p);
  N_broken = reduce_sum(get_N_broken(f,p), p);
  N_links = reduce_sum(get_broken_links(f,p), p);



  
  if(!p.rank) {
    printf("%06d,%6lf,%6lf,%6lf,%6lf,%6lf,%6lf,%6lf,"
	   "%6lf,%6lf,%6lf,%6lf,%6lf,%6lf,"
	   "%6lf,%6lf,%6lf,%6lf,%6lf,"
	   "%d,%6lf,%6lf,%6lf,%6lf,%lli,%lli\n",
	   step,
	   t_sim,
	   energies[0],
	   energies[1],
	   energies[2],
	   energies[3],
	   energies[4],
	   energies[5],
	   energies[6],
	   energies[7],
	   energies[8],
	   energies[9],
	   pressure_sum[0],
	   pressure_sum[1],
	   T_sum[0],
	   T_sum[1],
	   Vsq_sum[0],
	   Vsq_sum[1],
	   gwen,
	   bcount,
	   s_max,
	   gamma_max,
	   curlJ_tot,
	   divJ_tot,
	   N_broken,
	   N_links);
  }

}

/** Returns the largest zone-centred gamma factor found anywhere
 * in the simulation box.
 */
float get_gamma_max(hydro_fields f, hydro_params p) {
#ifndef SCALAR

  int x, y, z, xmax;

  float gmax = f.W[0][0][0];

  float gtest;

  // Just search for maxmimum
  for(x = 1; x <= p.slicex; x++) {
    for(y = 1; y <= p.slicey; y++) {
      for(z = 0; z < p.Lz; z++) {
	gtest = f.W[x][y][z];
	if(gtest > gmax) {
	  gmax = gtest;
	}

      }
    }
  }

  return gmax;

#else
  return 0.0;
#endif
}

/** Returns the largest zone-centred damping factor s found anywhere
 * in the simulation box.
 */
float get_s_max(hydro_fields f, hydro_params p) {
  int x, y, z;
  float smax = 0;
  float stest;

  for(x = 1; x <= p.slicex; x++) {
    for(y = 1; y <= p.slicey; y++) {
      for(z = 0; z < p.Lz; z++) {
#ifdef DIMENSIONLESS
#ifndef SCALAR
	stest = 0.5*p.dt*(p.C*f.phi[x][y][z]*f.phi[x][y][z]
		       /f.T[x][y][z])*f.W[x][y][z];
#else
	stest = 0.5*p.dt*(p.C*f.phi[x][y][z]*f.phi[x][y][z]
		       /p.Tconst);
#endif //!SCALAR
#else
#ifndef SCALAR
	stest = 0.5*p.dt*p.C*f.W[x][y][z];
#else
	stest = 0.5*p.dt*p.C;
#endif //!SCALAR
#endif //DIMENSIONLESS
	if(stest>smax) {
	  smax = stest;
	}
      }
    }
  }

  return smax;

}

/** Calculate sum over V^2 within this core. Populate array with result split
 *  into symmetric and broken phase.
 */
void calculate_Vsq_sum(hydro_fields f, hydro_params p, float *V_sq_sum){
  
  int x, y, z;

  float Vsq_symm = 0;
  float Vsq_broken = 0;
  float phi_broken;
  
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
	  Vsq_symm += (f.V[0][x][y][z]*f.V[0][x][y][z]
		       + f.V[1][x][y][z]*f.V[1][x][y][z]
		       + f.V[2][x][y][z]*f.V[2][x][y][z]);

	}
	else {
	  Vsq_broken += (f.V[0][x][y][z]*f.V[0][x][y][z]
			 + f.V[1][x][y][z]*f.V[1][x][y][z]
			 + f.V[2][x][y][z]*f.V[2][x][y][z]);
	}
      }
    }
  }
#endif // !SCALAR

  V_sq_sum[0] = Vsq_symm;
  V_sq_sum[1] = Vsq_broken;

  
}



/** Find number of sites in the broken phase. Allows for calculation
 * of volume of broken phase, and therefore wall speed.
 */
long long get_N_broken(hydro_fields f, hydro_params p){
  int x, y, z;
  long long N_broken = 0;
  float phi_broken;

  for(x = 1; x <= p.slicex; x++) {
    for(y = 1; y <= p.slicey; y++) {
      for(z = 0; z < p.Lz; z++) {
#ifdef BAG
	if(f.phi[x][y][z] > p.phi_0/2.){
	  N_broken += 1;
	}
#else
	phi_broken =  (p.alpha*f.T[x][y][z]
		       + sqrt((p.alpha*p.Tconst)*(p.alpha*f.T[x][y][z])
			      - 4.0*p.lambda*p.gamma
			      *(f.T[x][y][z]*f.T[x][y][z] - p.T0*p.T0))
		       )/(2.0*p.lambda);

	if(f.phi[x][y][z] > phi_broken/2.){
	  N_broken += 1;
	}
#endif
      }
    }
  }

  return N_broken;

}

/** Find number of links between broken and symmetric phase.  Allows
 * for calculation of surface area of bubbles, and therefore wall
 * speed.
 */
long long get_broken_links(hydro_fields f, hydro_params p){
  int x, y, z;
  long long N_links = 0;
  float phi_broken;

  for(x = 1; x <= p.slicex; x++) {
    for(y = 1; y <= p.slicey; y++) {
      for(z = 0; z < p.Lz; z++) {
#ifdef BAG
	if((f.phi[x+1][y][z] - p.phi_0/2.)*(f.phi[x][y][z] - p.phi_0/2.) < 0){
	  N_links += 1;
	}
	if((f.phi[x][y+1][z] - p.phi_0/2.)*(f.phi[x][y][z] - p.phi_0/2.) < 0){
	  N_links += 1;
	}
	if((f.phi[x][y][(z+1)%p.Lz] - p.phi_0/2.)
	   *(f.phi[x][y][z] - p.phi_0/2.) < 0){
	  N_links += 1;
	}
#else
	phi_broken =  (p.alpha*f.T[x][y][z]
		       + sqrt((p.alpha*p.Tconst)*(p.alpha*f.T[x][y][z])
			      - 4.0*p.lambda*p.gamma
			      *(f.T[x][y][z]*f.T[x][y][z] - p.T0*p.T0))
		       )/(2.0*p.lambda);
	if((f.phi[x+1][y][z] - phi_broken/2.)
	   *(f.phi[x][y][z] - phi_broken/2.) < 0){
	  N_links += 1;
	}
	if((f.phi[x][y+1][z] - phi_broken/2.)
	   *(f.phi[x][y][z] - phi_broken/2.) < 0){
	  N_links += 1;
	}
	if((f.phi[x][y][(z+1)%p.Lz] - phi_broken/2.)
	   *(f.phi[x][y][z] - phi_broken/2.) < 0){
	  N_links += 1;
	}
#endif // BAG
      }
    }
  }

  return N_links;
}

/** Dumps a field to stderr. Expects field to have N entries.
 * For debugging purposes...
 */
void dump(float *field, hydro_params p) {
  int x;

  fprintf(stderr,"%g", field[0]);
  for(x=1;x<p.N;x++) {
    fprintf(stderr,", %g", field[x]);
  }
  fprintf(stderr,"\n");
}



/** Calculate a histogram of the field, and store in a file
 * labelled by the timestep.
 */
void histo_field(float ***field, hydro_params p, int step) {
  int x, y, z;

  float fmax = 0.0;
  float fmin = 0.0;

  float ftest;

  float start = clock();


  fmax = field[0][0][0];
  fmin = field[0][0][0];

  for(x = 1; x <= p.slicex; x++) {
    for(y = 1; y <= p.slicey; y++) {
      for(z = 0; z < p.Lz; z++) {

	ftest = field[x][y][z];
	if(ftest > fmax) {
	  fmax = ftest;
	} else if (ftest < fmin) {
	  fmin = ftest;
	}

      }
    }
  }

  float overall_max, overall_min;

  overall_max = reduce_max(fmax, p);
  overall_min = reduce_min(fmin, p);

  /* Avoid fencepost errors by making overall_max and overall_min
     slightly larger (and smaller, respectively) than the largest (and
     smallest) field values found in the lattice. */
  const double bin_epsilon = 0.01*(overall_max - overall_min);
  overall_min = overall_min - bin_epsilon;
  overall_max = overall_max + bin_epsilon;



  int nbins = 100;

  float df = (overall_max - overall_min)/((float)nbins);

  if(isnan(df)) {
    printf0(p, "Max: %lf, Min: %lf, df: %lf\n",
	    overall_max, overall_min, df);
    printf0(p, "Not doing histogram: bin width not a number\n");

    return;
  }


  if(fabs(df) < 0.001) {
    printf0(p, "Max: %lf, Min: %lf, df: %lf\n",
	    overall_max, overall_min, df);
    printf0(p, "Not doing histogram: bin width too small "
	    "(overflow would occur)\n");

    return;
  }

  int *count = (int *)malloc(nbins*sizeof(int));

  int i;

  for(i=0; i<nbins; i++) {
    count[i] = 0;
  }


  for(x = 1; x <= p.slicex; x++) {
    for(y = 1; y <= p.slicey; y++) {
      for(z = 0; z < p.Lz; z++) {
	count[(int)floor((field[x][y][z] - overall_min)/df)]++;
      }
    }
  }

  printf0(p, "Doing histogram, %d bins in [%lf,%lf]\n",
	  nbins, overall_min, overall_max);


  int ntemp;

  for(i=0; i<nbins; i++) {
    ntemp = count[i];
    count[i] = reduce_sum_int(ntemp, p);
  }

  if(!p.rank) {

    FILE *fp;

    char histodest[200];

    sprintf(histodest,"histo-%d.txt",step);

    fp = fopen(histodest, "w");

    for(i=0; i<nbins; i++) {
      fprintf(fp, "%lf %d\n", overall_min + df*((float)i), count[i]);
    }

    fclose(fp);
  }

  float end = clock();

  printf0(p, "Histogram stuff took %lf\n",
	  ((float) (end - start)) / CLOCKS_PER_SEC);


  free(count);
}


/** *Average* components of the stress-energy tensor for the scalar
 * field to leading order in the metric perturbation.
 */
void didj(float *cpts, hydro_fields f, hydro_params p) {

  float cpts_here[TENSOR_CPTS];

  int x, y, z;
  int i;

  for(i=0; i<TENSOR_CPTS; i++) {
    cpts_here[i] = 0.0;
  }

  for(x = 1; x <= p.slicex; x++) {
    for(y = 1; y <= p.slicey; y++) {
      for(z = 0; z < p.Lz; z++) {

          cpts_here[CPT_11] +=
            0.25*((f.phi[x+1][y][z] - f.phi[x-1][y][z])/p.dx)
            *((f.phi[x+1][y][z] - f.phi[x-1][y][z])/p.dx);

          cpts_here[CPT_21] +=
            0.25*((f.phi[x+1][y][z] - f.phi[x-1][y][z])/p.dx)
            *((f.phi[x][y+1][z] - f.phi[x][y-1][z])/p.dx);

          cpts_here[CPT_31] +=
            0.25*((f.phi[x][y][(z+1)%p.Lz]
		   - f.phi[x][y][(z-1+p.Lz)%p.Lz])/p.dx)
            *((f.phi[x+1][y][z] - f.phi[x-1][y][z])/p.dx);

          cpts_here[CPT_22] +=
            0.25*((f.phi[x][y+1][z]
		   - f.phi[x][y-1][z])/p.dx)
            *((f.phi[x][y+1][z] - f.phi[x][y-1][z])/p.dx);

          cpts_here[CPT_32] +=
            0.25*((f.phi[x][y][(z+1)%p.Lz]
		   - f.phi[x][y][(z-1+p.Lz)%p.Lz])/p.dx)
            *((f.phi[x][y+1][z] - f.phi[x][y-1][z])/p.dx);

          cpts_here[CPT_33] +=
            0.25*((f.phi[x][y][(z+1)%p.Lz]
		   - f.phi[x][y][(z-1+p.Lz)%p.Lz])/p.dx)
            *((f.phi[x][y][(z+1)%p.Lz] - f.phi[x][y][(z-1+p.Lz)%p.Lz])/p.dx);

      }
    }
  }


  for(i=0; i<TENSOR_CPTS; i++) {
    cpts[i] = reduce_sum(cpts_here[i], p)/((float)(p.Lx*p.Ly*p.Lz));
  }

}

/** Returns the maximum and its position of a specified field on the lattice.
 *
 * For debugging purposes.
 * Only ensures that the master node has the correct value.
 * Additional boolean int to determine if we want to find maximum
 * of absolute value or just maximum.
 */
value_loc find_max_loc(float ***field, hydro_params p, int abs_max){
  float field_max= -999999999999;
  int x, y, z;
  int max_loc[3] = {-1,-1,-1};
  value_rank out;
  value_loc max_and_loc;

  for(x = 1; x <= p.slicex; x++) {
    for(y = 1; y <= p.slicey; y++) {
      for(z = 0; z < p.Lz; z++) {
	if(fabs(field[x][y][z]) > field_max && abs_max){
	  field_max = fabs(field[x][y][z]);
	  max_loc[0] = x + p.shiftx - 1;
	  max_loc[1] = y + p.shifty - 1;
	  max_loc[2] = z;
	}
	else if(field[x][y][z] > field_max && !abs_max){
	  field_max = field[x][y][z];
	  max_loc[0] = x + p.shiftx - 1;
	  max_loc[1] = y + p.shifty - 1;
	  max_loc[2] = z;
	}
      }
    }
  }

  out = reduce_maxloc(field_max, p);
  if(!out.rank){
  }
  else if(p.rank == out.rank){
    MPI_Send(max_loc, 3, MPI_INT, 0, p.rank, MPI_COMM_WORLD);
  }
  else if(!p.rank){
    MPI_Recv(max_loc, 3, MPI_INT, out.rank, out.rank, MPI_COMM_WORLD,
	     MPI_STATUS_IGNORE);
  }

  max_and_loc.value = out.value;
  max_and_loc.loc[0] = max_loc[0];
  max_and_loc.loc[1] = max_loc[1];
  max_and_loc.loc[2] = max_loc[2];

  return max_and_loc;
}

/** Returns the minimum and its position of a specified field on the lattice.
 *
 * For debugging purposes.
 * Only ensures that the master node has the correct value.
 * Additional boolean int to determine if we want to find minimum
 * of absolute value or just minimum.
 */
value_loc find_min_loc(float ***field, hydro_params p, int abs_min){
  float field_min= +999999999999;
  int x, y, z;
  int min_loc[3] = {-1,-1,-1};
  value_rank out;
  value_loc min_and_loc;

  for(x = 1; x <= p.slicex; x++) {
    for(y = 1; y <= p.slicey; y++) {
      for(z = 0; z < p.Lz; z++) {
	if(fabs(field[x][y][z]) < field_min && abs_min){
	  field_min = fabs(field[x][y][z]);
	  min_loc[0] = x + p.shiftx - 1;
	  min_loc[1] = y + p.shifty - 1;
	  min_loc[2] = z;
	}
	else if(field[x][y][z] < field_min && !abs_min){
	  field_min = field[x][y][z];
	  min_loc[0] = x + p.shiftx - 1;
	  min_loc[1] = y + p.shifty - 1;
	  min_loc[2] = z;
	}
      }
    }
  }

  out = reduce_minloc(field_min, p);
  if(!out.rank){
  }
  else if(p.rank == out.rank){
    MPI_Send(min_loc, 3, MPI_INT, 0, p.rank, MPI_COMM_WORLD);
  }
  else if(!p.rank){
    MPI_Recv(min_loc, 3, MPI_INT, out.rank, out.rank, MPI_COMM_WORLD,
	     MPI_STATUS_IGNORE);
  }

  min_and_loc.value = out.value;
  min_and_loc.loc[0] = min_loc[0];
  min_and_loc.loc[1] = min_loc[1];
  min_and_loc.loc[2] = min_loc[2];

  return min_and_loc;
}

/** Dumps maximum and minimum values of all fields and their locations
 * on the lattice.
 *
 * For debugging purposes only, very costly.
 *
 */
void dump_max_min(hydro_fields f, hydro_params p){

  value_loc val_loc;

  val_loc = find_max_loc(f.phi, p, 0);

  printf0(p,"Max of phi = %g at (%d, %d, %d) \n", val_loc.value,
	  val_loc.loc[0], val_loc.loc[1], val_loc.loc[2]);

  val_loc = find_min_loc(f.phi, p, 0);

  printf0(p,"Min of phi = %g at (%d, %d, %d) \n", val_loc.value,
	  val_loc.loc[0], val_loc.loc[1], val_loc.loc[2]);

  val_loc = find_max_loc(f.pi, p, 0);

  printf0(p,"Max of pi = %g at (%d, %d, %d) \n", val_loc.value,
	  val_loc.loc[0], val_loc.loc[1], val_loc.loc[2]);

  val_loc = find_min_loc(f.pi, p, 0);

  printf0(p,"Min of pi = %g at (%d, %d, %d) \n", val_loc.value,
	  val_loc.loc[0], val_loc.loc[1], val_loc.loc[2]);

#ifndef SCALAR

  val_loc = find_max_loc(f.T, p, 0);

  printf0(p,"Max of T = %g at (%d, %d, %d) \n", val_loc.value,
	  val_loc.loc[0], val_loc.loc[1], val_loc.loc[2]);

  val_loc = find_min_loc(f.T, p, 0);

  printf0(p,"Min of T = %g at (%d, %d, %d) \n", val_loc.value,
	  val_loc.loc[0], val_loc.loc[1], val_loc.loc[2]);

  val_loc = find_max_loc(f.E, p, 0);

  printf0(p,"Max of E = %g at (%d, %d, %d) \n", val_loc.value,
	  val_loc.loc[0], val_loc.loc[1], val_loc.loc[2]);

  val_loc = find_min_loc(f.E, p, 0);

  printf0(p,"Min of E = %g at (%d, %d, %d) \n", val_loc.value,
	  val_loc.loc[0], val_loc.loc[1], val_loc.loc[2]);

  val_loc = find_max_loc(f.W, p, 0);

  printf0(p,"Max of W = %g at (%d, %d, %d) \n", val_loc.value,
	  val_loc.loc[0], val_loc.loc[1], val_loc.loc[2]);

  val_loc = find_min_loc(f.W, p, 0);

  printf0(p,"Min of W = %g at (%d, %d, %d) \n", val_loc.value,
	  val_loc.loc[0], val_loc.loc[1], val_loc.loc[2]);

  val_loc = find_max_loc(f.kappa, p, 0);

  printf0(p,"Max of kappa = %g at (%d, %d, %d) \n", val_loc.value,
	  val_loc.loc[0], val_loc.loc[1], val_loc.loc[2]);

  val_loc = find_min_loc(f.kappa, p, 0);

  printf0(p,"Min of kappa = %g at (%d, %d, %d) \n", val_loc.value,
	  val_loc.loc[0], val_loc.loc[1], val_loc.loc[2]);

  val_loc = find_max_loc(f.p, p, 0);

  printf0(p,"Max of p = %g at (%d, %d, %d) \n", val_loc.value,
	  val_loc.loc[0], val_loc.loc[1], val_loc.loc[2]);

  val_loc = find_min_loc(f.p, p, 0);

  printf0(p,"Min of p = %g at (%d, %d, %d) \n", val_loc.value,
	  val_loc.loc[0], val_loc.loc[1], val_loc.loc[2]);

  val_loc = find_max_loc(f.V[0], p, 1);

  printf0(p,"Max of |V[0]| = %g at (%d, %d, %d) \n", val_loc.value,
	  val_loc.loc[0], val_loc.loc[1], val_loc.loc[2]);

  val_loc = find_max_loc(f.V[1], p, 1);

  printf0(p,"Max of |V[1]| = %g at (%d, %d, %d) \n", val_loc.value,
	  val_loc.loc[0], val_loc.loc[1], val_loc.loc[2]);

  val_loc = find_max_loc(f.V[2], p, 1);

  printf0(p,"Max of |V[2]| = %g at (%d, %d, %d) \n", val_loc.value,
	  val_loc.loc[0], val_loc.loc[1], val_loc.loc[2]);

  val_loc = find_min_loc(f.V[0], p, 1);

  printf0(p,"Min of |V[0]| = %g at (%d, %d, %d) \n", val_loc.value,
	  val_loc.loc[0], val_loc.loc[1], val_loc.loc[2]);

  val_loc = find_min_loc(f.V[1], p, 1);

  printf0(p,"Min of |V[1]| = %g at (%d, %d, %d) \n", val_loc.value,
	  val_loc.loc[0], val_loc.loc[1], val_loc.loc[2]);

  val_loc = find_min_loc(f.V[2], p, 1);

  printf0(p,"Min of |V[2]| = %g at (%d, %d, %d) \n", val_loc.value,
	  val_loc.loc[0], val_loc.loc[1], val_loc.loc[2]);

  val_loc = find_max_loc(f.Z[0], p, 1);

  printf0(p,"Max of |Z[0]| = %g at (%d, %d, %d) \n", val_loc.value,
	  val_loc.loc[0], val_loc.loc[1], val_loc.loc[2]);

  val_loc = find_max_loc(f.Z[1], p, 1);

  printf0(p,"Max of |Z[1]| = %g at (%d, %d, %d) \n", val_loc.value,
	  val_loc.loc[0], val_loc.loc[1], val_loc.loc[2]);

  val_loc = find_max_loc(f.Z[2], p, 1);

  printf0(p,"Max of |Z[2]| = %g at (%d, %d, %d) \n", val_loc.value,
	  val_loc.loc[0], val_loc.loc[1], val_loc.loc[2]);

  val_loc = find_min_loc(f.Z[0], p, 1);

  printf0(p,"Min of |Z[0]| = %g at (%d, %d, %d) \n", val_loc.value,
	  val_loc.loc[0], val_loc.loc[1], val_loc.loc[2]);

  val_loc = find_min_loc(f.Z[1], p, 1);

  printf0(p,"Min of |Z[1]| = %g at (%d, %d, %d) \n", val_loc.value,
	  val_loc.loc[0], val_loc.loc[1], val_loc.loc[2]);

  val_loc = find_min_loc(f.Z[2], p, 1);

  printf0(p,"Min of |Z[2]| = %g at (%d, %d, %d) \n", val_loc.value,
	  val_loc.loc[0], val_loc.loc[1], val_loc.loc[2]);

  val_loc = find_max_loc(f.U[0], p, 1);

  printf0(p,"Max of |U[0]| = %g at (%d, %d, %d) \n", val_loc.value,
	  val_loc.loc[0], val_loc.loc[1], val_loc.loc[2]);

  val_loc = find_max_loc(f.U[1], p, 1);

  printf0(p,"Max of |U[1]| = %g at (%d, %d, %d) \n", val_loc.value,
	  val_loc.loc[0], val_loc.loc[1], val_loc.loc[2]);

  val_loc = find_max_loc(f.U[2], p, 1);

  printf0(p,"Max of |U[2]| = %g at (%d, %d, %d) \n", val_loc.value,
	  val_loc.loc[0], val_loc.loc[1], val_loc.loc[2]);

  val_loc = find_min_loc(f.U[0], p, 1);

  printf0(p,"Min of |U[0]| = %g at (%d, %d, %d) \n", val_loc.value,
	  val_loc.loc[0], val_loc.loc[1], val_loc.loc[2]);

  val_loc = find_min_loc(f.U[1], p, 1);

  printf0(p,"Min of |U[1]| = %g at (%d, %d, %d) \n", val_loc.value,
	  val_loc.loc[0], val_loc.loc[1], val_loc.loc[2]);

  val_loc = find_min_loc(f.U[2], p, 1);

  printf0(p,"Min of |U[2]| = %g at (%d, %d, %d) \n", val_loc.value,
	  val_loc.loc[0], val_loc.loc[1], val_loc.loc[2]);

#endif //!SCALAR

}
