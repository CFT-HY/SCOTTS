/* output.c
 * 
 * Simple things that we could in principle calculate offline, but
 * that are nice to know (energy-related quantites are in energy.c).
 */
#include "hydro.h"


/* float get_gamma_max(hydro_fields f, hydro_params p)
 *
 * Returns the largest zone-centred gamma factor found anywhere
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

/* float get_s_max(hydro_fields f, hydro_params p)
 *
 * Returns the largest zone-centred damping factor s found anywhere
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
		       /f.Tconst);
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
  
/* float get_veltot(hydro_fields f, hydro_params p)
 *
 * The sum of the fluid (3-)velocity everywhere. A strange quantity
 * on its own, but allows calculation of average fluid velocity.
 */
float get_veltot(hydro_fields f, hydro_params p) {

#ifndef SCALAR
  int x, y, z, xmax;

  
  float veltot = 0.0;

  // Just search for maxmimum
  for(x = 1; x <= p.slicex; x++) {
    for(y = 1; y <= p.slicey; y++) {
      for(z = 0; z < p.Lz; z++) {
	veltot += sqrt(f.V[0][x][y][z]*f.V[0][x][y][z]
	  + f.V[1][x][y][z]*f.V[1][x][y][z]
	  + f.V[2][x][y][z]*f.V[2][x][y][z]);

      }
    }
  }
  
  return veltot;
#else
  return 0.0;
#endif
}


/* void dump(float *field, hydro_params p) 
 *
 * Dumps a field to stderr. Expects field to have N entries.
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



/* void histo_field(float ***field, hydro_params p, int step)
 *
 * Calculate a histogram of the field, and store in a file
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


/* void didj(float *cpts, hydro_fields f, hydro_params p)
 *
 * *Average* components of the stress-energy tensor for the scalar
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

