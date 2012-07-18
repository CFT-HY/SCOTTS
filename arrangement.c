#include "hydro.h"

#ifdef MPI
#include <mpi.h>

// Inspired by the contents of setup_layout_generic.c
void layout(hydro_params p) {

  int rank, size;

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  static int n_primes = 8;
  static int prime[] = {2,3,5,7,11,13,17,19};

  int n;
  int nfactors[n_primes];
  int i = size;

  fprintf(stderr,"Putting %d sites on %d nodes\n", p.N, size);

  if( (p.Lx*p.Ly % size) > 0 ) {
    fprintf(stderr,"Number of (%d) nodes not a factor of Lx*Ly!\n", size);
    exit(-99);
  }

  for (n=0; n<n_primes; n++) {
    nfactors[n] = 0;
    while (i % prime[n] == 0) {
      i /= prime[n];
      nfactors[n]++;
    }
  }

  if(i != 1) {
    fprintf(stderr,"Can\'t figure out how to put %dx%d%d on %d nodes\n", p.Lx, p.Ly, p.Lz, size);
    exit(-100);
  }

  p.slicex = p.Lx;
  p.slicey = p.Ly;

  for(n=0;n<n_primes;n++)
    fprintf(stderr,"nfactors %d = %d\n", prime[n],nfactors[n]);

  for (n=n_primes-1; n>=0; n--) {
    for(i=0; i<nfactors[n]; i++) {
      
      if((p.slicex > p.slicey) && (p.slicex % prime[n] == 0)) {
	fprintf(stderr,"Dividing x (%d) by %d\n", p.slicex, prime[n]);
	p.slicex = p.slicex/prime[n];
      } else if(p.slicey % prime[n] == 0) {
	fprintf(stderr,"Dividing y (%d) by %d\n", p.slicey, prime[n]);
	p.slicey = p.slicey/prime[n];
      } else {
	fprintf(stderr,"Should not get here\n");
      }
    }
  }

  fprintf(stderr,"slicing: x=%d y=%d\n", p.slicex, p.slicey);
  fprintf(stderr,"meaning %d nodes in x, %d nodes in y\n", p.Lx/p.slicex, p.Ly/p.slicey);
  
  p.totalN = p.N;
  // number of sites + halos + double halos
  p.fieldN = p.N + p.Lz*(2*(p.slicex + p.slicey) + 4);
}




void generate_offsets(hydro_params p) {

  p.N = p.slicex*p.slicey*p.Lz;

  // Starting point in the box
  p.inner = (p.slicex-2)*(p.slicey-2)*p.Lz;

  // Sender boxes
  p.inner_xM = p.inner;
  p.inner_xP = p.inner_xM + (p.Lz*(p.slicey-1));
  p.inner_yM = p.inner_xP + (p.Lz*(p.slicey-1));
  p.inner_yP = p.inner_yM + (p.Lz*(p.slicex-1));

  // Double sender boxes
  p.inner_xMyM = p.inner_yP + (p.Lz*(p.slicex-1));
  p.inner_xMyP = p.inner_xMyM + p.Lz;
  p.inner_xPyM = p.inner_xMyP + p.Lz;
  p.inner_xPyP = p.inner_xPyM + p.Lz;

  fprintf(stderr,
	  "Checking: number of sites in box, %d, equals final sender box, %d\n",
	  p.N,
	  p.inner_xPyP + p.Lz);

  // Single haloes
  p.offset_xM = p.offset_xPyP + (p.Lz*p.slicex);
  p.offset_xP = p.offset_xM + (p.Lz*p.slicey);
  p.offset_yM = p.offset_xP + (p.Lz*p.slicey);
  p.offset_yP = p.offset_yM + (p.Lz*p.slicex);

  // Double haloes
  p.offset_xMyM = p.offset_yP + (p.Lz*p.slicex);
  p.offset_xMyP = p.offset_xMyM + p.Lz;
  p.offset_xPyM = p.offset_xMyP + p.Lz;
  p.offset_xPyP = p.offset_xPyM + p.Lz;

  fprintf(stderr,
	  "Checking: number of sites including halo, %d, equals final offset, %d\n",
	  p.fieldN,
	  p.offset_xPyP + p.Lz);

}

int local_iix(int x, int y, int z, hydro_params p) {

  if (x < -1 || x > p.slicex || y < -1 || y > p.slicey) {
    // Can't handle
    fprintf(stderr,"Unable to handle local location (%d,%d,%d)\n", x, y, z);
    exit(-45);

  } else if (x == -1) {
    if(y == -1) {
      return p.offset_xMyM + ((z+p.Lz)%p.Lz);
    } else if (y == p.slicey) {
      return p.offset_xMyP + ((z+p.Lz)%p.Lz);
    } else {
      return p.offset_xM + y*p.Lz + ((z+p.Lz)%p.Lz);
    }
    
  } else if (x == p.slicex) {
    if(y == -1) {
      return p.offset_xPyM + ((z+p.Lz)%p.Lz);
    } else if (y == p.slicey) {
      return p.offset_xPyP + ((z+p.Lz)%p.Lz);
    } else {
      return p.offset_xP + y*p.Lz + ((z+p.Lz)%p.Lz);
    }
  } else if (y == -1) {
    return p.offset_yM + x*p.Lz + ((z+p.Lz)%p.Lz);
  } else if (y == p.slicey) {
    return p.offset_yP + x*p.Lz + ((z+p.Lz)%p.Lz);



  }  else if (x == 0) {
    if(y == 0) {
      return p.inner_xMyM + ((z+p.Lz)%p.Lz);
    } else if (y == p.slicey-1) {
      return p.inner_xMyP + ((z+p.Lz)%p.Lz);
    } else {
      return p.inner_xM + (y-1)*p.Lz + ((z+p.Lz)%p.Lz);
    }
  } else if (x == p.slicex-1) {
    if(y == 0) {
      return p.inner_xPyM + ((z+p.Lz)%p.Lz);
    } else if (y == p.slicey-1) {
      return p.inner_xPyP + ((z+p.Lz)%p.Lz);
    } else {
      return p.inner_xP + (y-1)*p.Lz + ((z+p.Lz)%p.Lz);
    }
  } else if (y == 0) {
    return p.inner_yM + (x-1)*p.Lz + ((z+p.Lz)%p.Lz);
  } else if (y == p.slicey-1) {
    return p.inner_yP + (x-1)*p.Lz + ((z+p.Lz)%p.Lz);
  } else {
    // In the bulk
    return (x-1)*(p.slicey-2)*p.Lz + (y-1)*p.Lz + z;
  }
}

int **init_nb(hydro_params p) {
  int x, y, z;

  int **nb;

  // set up nb lookup table                                                                        
  // each node does the same thing
  nb = (int **) malloc(p.N*sizeof(int *));


  for(x=0; x<p.slicex; x++) {
    for(y=0; y<p.slicey; y++) {
      for(z=0; z<p.Lz; z++) {

        nb[local_iix(x,y,z,p)] = (int *)malloc(14*sizeof(int));

        nb[local_iix(x,y,z,p)][0] = local_iix(x+1, y, z, p);
        nb[local_iix(x,y,z,p)][1] = local_iix(x-1, y, z, p);
        nb[local_iix(x,y,z,p)][2] = local_iix(x, y+1, z, p);
        nb[local_iix(x,y,z,p)][3] = local_iix(x, y-1, z, p);
        nb[local_iix(x,y,z,p)][4] = local_iix(x, y, z+1, p);
        nb[local_iix(x,y,z,p)][5] = local_iix(x, y, z-1, p);

        nb[iix(x, y, z, p)][DIR_02] = local_iix(x+1, y+1, z, p);
        nb[iix(x, y, z, p)][DIR_04] = local_iix(x+1, y, z+1, p);
        nb[iix(x, y, z, p)][DIR_24] = local_iix(x, y+1, z+1, p);
        nb[iix(x, y, z, p)][DIR_13] = local_iix(x-1, y-1, z, p);
        nb[iix(x, y, z, p)][DIR_15] = local_iix(x-1, y, z-1, p);
        nb[iix(x, y, z, p)][DIR_35] = local_iix(x, y-1, z-1, p);

        nb[iix(x, y, z, p)][DIR_024] = local_iix(x+1, y+1, z+1, p);
        nb[iix(x, y, z, p)][DIR_135] = local_iix(x-1, y-1, z-1, p);

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


void halo_field(double *field, hydro_params p) {
  // set up

  // <<

  // >>

  // /\ /\

  // \/ \/

  // [ ']

  // [, ]

  // [' ]

  // [ ,]
}


#else // MPI


int iix(int x, int y, int z, hydro_params p) {
  return ((x+p.Lx)%p.Lx)*p.Ly*p.Lz + ((y+p.Ly)%p.Ly)*p.Lz + ((z+p.Lz)%p.Lz);
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

int get_z(int n, hydro_params p) {
  return (n % p.Lz);
}

int get_y(int n, hydro_params p) {
  n = n - (n % p.Lz);
  n = n / p.Lz;
  return (n % p.Ly);
}

int get_x(int n, hydro_params p) {
  n = n - (n % p.Lz);
  n = n / p.Lz;
  n = n - (n % p.Ly);
  return n / p.Ly;
}


void layout(hydro_params p) {
  fprintf(stderr,"Layout: doing nothing: no parallelism\n");
}

#endif // MPI
