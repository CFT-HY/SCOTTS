#include "hydro.h"

#ifdef MPI
#include <mpi.h>

// Inspired by the contents of setup_layout_generic.c
void layout(hydro_params *p) {


  static int n_primes = 8;
  static int prime[] = {2,3,5,7,11,13,17,19};

  int n;
  int nfactors[n_primes];
  int i = p->size;

  if(!p->rank)
    fprintf(stderr,"Putting %d sites on %d nodes\n", p->N, p->size);

  if( (p->Lx*p->Ly % p->size) > 0 ) {
    if(!p->rank)
      fprintf(stderr,"Number of (%d) nodes not a factor of Lx*Ly!\n", p->size);
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
    if(!p->rank)
      fprintf(stderr,"Can\'t figure out how to put %dx%d%d on %d nodes\n", p->Lx, p->Ly, p->Lz, p->size);
    exit(-100);
  }

  p->slicex = p->Lx;
  p->slicey = p->Ly;

  for(n=0;n<n_primes;n++)
    if(!p->rank)
      fprintf(stderr,"nfactors %d = %d\n", prime[n],nfactors[n]);

  for (n=n_primes-1; n>=0; n--) {
    for(i=0; i<nfactors[n]; i++) {
      
      if((p->slicex > p->slicey) && (p->slicex % prime[n] == 0)) {
	if(!p->rank)
	  fprintf(stderr,"Dividing x (%d) by %d\n", p->slicex, prime[n]);
	p->slicex = p->slicex/prime[n];
      } else if(p->slicey % prime[n] == 0) {
	if(!p->rank)
	  fprintf(stderr,"Dividing y (%d) by %d\n", p->slicey, prime[n]);
	p->slicey = p->slicey/prime[n];
      } else {
	if(!p->rank)
	  fprintf(stderr,"Should not get here\n");
	exit(100);
      }
    }
  }


  int nx = p->Lx/p->slicex;
  int ny = p->Ly/p->slicey;

  if(!p->rank) {
    fprintf(stderr,"slicing: x=%d y=%d\n", p->slicex, p->slicey);
    fprintf(stderr,"meaning %d nodes in x, %d nodes in y\n", nx, ny);
  }  

  p->totalN = p->N;
  // number of sites + halos + double halos
  p->fieldN = (p->slicex*p->slicey*p->Lz) + p->Lz*(2*(p->slicex + p->slicey) + 4);
  if(!p->rank)
    fprintf(stderr,"Each node will need %d = %d + %d entries to store halo\n",
	    p->fieldN,
	    p->totalN,
	    p->Lz*(2*(p->slicex + p->slicey) + 4));


  p->myposx = p->rank % nx;
  p->myposy = (p->rank - p->myposx)/nx;

  p->rank_xM = p->myposy*nx + ((p->myposx-1+nx)%nx);
  p->rank_xP = p->myposy*nx + ((p->myposx+1+nx)%nx);
  p->rank_yM = ((p->myposy - 1 + ny)%ny)*nx + p->myposx;
  p->rank_yP = ((p->myposy + 1 + ny)%ny)*nx + p->myposx;

  p->rank_xMyM = ((p->myposy - 1 + ny)%ny)*nx + ((p->myposx-1+nx)%nx);
  p->rank_xMyP = ((p->myposy + 1 + ny)%ny)*nx + ((p->myposx-1+nx)%nx);
  p->rank_xPyM = ((p->myposy - 1 + ny)%ny)*nx + ((p->myposx+1+nx)%nx);
  p->rank_xPyP = ((p->myposy + 1 + ny)%ny)*nx + ((p->myposx+1+nx)%nx);

  fprintf(stderr, "Node of rank %d, node position is (%d, %d)\n", p-> rank, p->myposx, p->myposy);
  fprintf(stderr,
	  "Node of rank %d, neighbours are < %d, > %d, /\\ %d, \\/ %d, |_ %d |^ %d ^| %d _| %d\n",
	  p->rank, p->rank_xM, p->rank_xP, p->rank_yM, p->rank_yP, p->rank_xMyM, p->rank_xMyP,
	  p->rank_xPyM, p->rank_xPyP);
}




void generate_offsets(hydro_params *p) {

  p->N = p->slicex*p->slicey*p->Lz;
  if(!p->rank)
    fprintf(stderr, "Resizing p->N to %d\n", p->N);

  // Starting point in the box
  p->inner = (p->slicex-2)*(p->slicey-2)*p->Lz;

  // Sender boxes
  p->inner_xM = p->inner;
  p->inner_xP = p->inner_xM + (p->Lz*(p->slicey-2));
  p->inner_yM = p->inner_xP + (p->Lz*(p->slicey-2));
  p->inner_yP = p->inner_yM + (p->Lz*(p->slicex-2));

  // Double sender boxes
  p->inner_xMyM = p->inner_yP + (p->Lz*(p->slicex-2));
  p->inner_xMyP = p->inner_xMyM + p->Lz;
  p->inner_xPyM = p->inner_xMyP + p->Lz;
  p->inner_xPyP = p->inner_xPyM + p->Lz;

  if(!p->rank)
    fprintf(stderr,
	    "Checking: number of sites in box, %d, equals final sender box, %d\n",
	    p->N,
	    p->inner_xPyP + p->Lz);
  
  // Single haloes
  p->offset_xM = p->inner_xPyP + p->Lz;
  p->offset_xP = p->offset_xM + (p->Lz*(p->slicey-2));
  p->offset_yM = p->offset_xP + (p->Lz*(p->slicey-2));
  p->offset_yP = p->offset_yM + (p->Lz*(p->slicex-2));

  // Endposts
  p->offset_xMyeM = p->offset_yP + (p->Lz*(p->slicex-2));
  p->offset_xMyeP = p->offset_xMyeM + p->Lz;
  p->offset_xPyeM = p->offset_xMyeP + p->Lz;
  p->offset_xPyeP = p->offset_xPyeM + p->Lz;
  p->offset_yMxeM = p->offset_xPyeP + p->Lz;
  p->offset_yMxeP = p->offset_yMxeM + p->Lz;
  p->offset_yPxeM = p->offset_yMxeP + p->Lz;
  p->offset_yPxeP = p->offset_yPxeM + p->Lz;

  // Double haloes
  p->offset_xMyM = p->offset_yPxeP + p->Lz;
  p->offset_xMyP = p->offset_xMyM + p->Lz;
  p->offset_xPyM = p->offset_xMyP + p->Lz;
  p->offset_xPyP = p->offset_xPyM + p->Lz;

  if(!p->rank)
    fprintf(stderr,
	    "Checking: number of sites including halo, %d, equals final offset, %d\n",
	    p->fieldN,
	    p->offset_xPyP + p->Lz);

}

int local_iix(int x, int y, int z, hydro_params p) {

  if (x < -1 || x > p.slicex || y < -1 || y > p.slicey) {
    // Can't handle
    fprintf(stderr,"Unable to handle local location (%d,%d,%d)\n", x, y, z);
    exit(-45);
  }


  // Deal with offsets (locations of halo material)

  // 1. Off left hand end of box
  if (x == -1) {

    // 1A Bottom left, double halo
    if(y == -1) {
      return p.offset_xMyM + ((z+p.Lz)%p.Lz);
    }

    // 1B Top left, double halo
    if (y == p.slicey) {
      return p.offset_xMyP + ((z+p.Lz)%p.Lz);
    }

    // 1C Bottom left, one up, endpost
    if (y == 0) {
      return p.offset_xMyeM + ((z+p.Lz)%p.Lz);
    }

    // 1D Top left, one down, endpost
    if (y == p.slicey-1) {
      return p.offset_xMyeP + ((z+p.Lz)%p.Lz);
    }

    // 1E Left, but not endpost
    return p.offset_xM + (y-1)*p.Lz + ((z+p.Lz)%p.Lz);
  }
    
  
  // 2. Off right hand end of box
  if (x == p.slicex) {

    // 2A Bottom right, double halo
    if(y == -1) {
      return p.offset_xPyM + ((z+p.Lz)%p.Lz);
    }

    // 2B Top right, double halo
    if (y == p.slicey) {
      return p.offset_xPyP + ((z+p.Lz)%p.Lz);
    } 

    // 2C Bottom right, one up, endpost
    if (y == 0) {
      return p.offset_xPyeM + ((z+p.Lz)%p.Lz);
    }

    // 2D Top right, one down, endpost
    if (y == p.slicey-1) {
      return p.offset_xPyeP + ((z+p.Lz)%p.Lz);
    }

    // 2E Right, but not endpost
    return p.offset_xP + (y-1)*p.Lz + ((z+p.Lz)%p.Lz);
   
  }

  // 3. Off bottom of box
  if (y == -1) {

    // 3A Bottom left, but one along, endpost
    if(x == 0) {
      return p.offset_yMxeM + ((z+p.Lz)%p.Lz);
    }

    // 3B Bottom right, one along, endpost
    if (x == p.slicex - 1) {
      return p.offset_yMxeP + ((z+p.Lz)%p.Lz);
    } 

    // 3C Bottom, but not end post
    return p.offset_yM + (x-1)*p.Lz + ((z+p.Lz)%p.Lz);

  }

  // 4. Off top of box
  if (y == p.slicey) {

    // 4A Top left, but one along, endpost
    if(x == 0) {
      return p.offset_yPxeM + ((z+p.Lz)%p.Lz);
    }

    // 4B Top right, but one along, endpost
    if (x == p.slicex - 1) {
      return p.offset_yPxeP + ((z+p.Lz)%p.Lz);
    } 

    // 4C Top, but not endpost
    return p.offset_yP + (x-1)*p.Lz + ((z+p.Lz)%p.Lz);
  }



  // Deal with inner edges of box

  // 1. Left hand side of box
  if (x == 0) {

    // 1A Bottom left corner, double halo/endpost
    if(y == 0) {
      return p.inner_xMyM + ((z+p.Lz)%p.Lz);
    }

    // 1B Top left corner, double halo/endpost
    if (y == p.slicey-1) {
      return p.inner_xMyP + ((z+p.Lz)%p.Lz);
    }

    // 1C Left, but not endpost
    return p.inner_xM + (y-1)*p.Lz + ((z+p.Lz)%p.Lz);

  }

  // 2. Right hand side of box
  if (x == p.slicex-1) {

    // 2A Bottom right corner, double halo/endpost
    if(y == 0) {
      return p.inner_xPyM + ((z+p.Lz)%p.Lz);
    }

    // 2B Top right corner, double halo/endpost
    if (y == p.slicey-1) {
      return p.inner_xPyP + ((z+p.Lz)%p.Lz);
    }
    
    // 2C Right, but not endpost
    return p.inner_xP + (y-1)*p.Lz + ((z+p.Lz)%p.Lz);
  }

  // 3. Bottom of box, but not bottom right or bottom left
  if (y == 0) {
    return p.inner_yM + (x-1)*p.Lz + ((z+p.Lz)%p.Lz);
  }

  // 4. Top of box, but not top right or top left
  if (y == p.slicey-1) {
    return p.inner_yP + (x-1)*p.Lz + ((z+p.Lz)%p.Lz);
  } 

  // 5. In the bulk
  return (x-1)*(p.slicey-2)*p.Lz + (y-1)*p.Lz + ((z+p.Lz)%p.Lz);
}


int **init_nb(hydro_params *p) {
  int x, y, z;

  int **nb;

  generate_offsets(p);

  p->inverse = (int **) malloc(p->N*sizeof(int *));

  for(x=0; x<p->slicex; x++) {
    for(y=0; y<p->slicey; y++) {
      for(z=0; z<p->Lz; z++) {
	
	//	fprintf(stderr,"Allocating %d\n", local_iix(x,y,z,*p));
        p->inverse[local_iix(x,y,z,*p)] = (int *)malloc(3*sizeof(int));
	p->inverse[local_iix(x,y,z,*p)][0] = x+(p->slicex*p->myposx);
	//	fprintf(stderr,"Entry %d is %d\n", local_iix(x,y,z,*p), p->inverse[local_iix(x,y,z,*p)][0]);
	p->inverse[local_iix(x,y,z,*p)][1] = y+(p->slicey*p->myposy);
	p->inverse[local_iix(x,y,z,*p)][2] = z;
      }
    }
  }

  //  fprintf(stderr,"Done allocating");
  //  fprintf(stderr,"Checking 0: %d\n", p->inverse[local_iix(0,0,0,*p)][0]);

  // set up nb lookup table                                                                        
  // each node does the same thing
  nb = (int **) malloc(p->N*sizeof(int *));


  for(x=0; x<p->slicex; x++) {
    for(y=0; y<p->slicey; y++) {
      for(z=0; z<p->Lz; z++) {

        nb[local_iix(x,y,z,*p)] = (int *)malloc(14*sizeof(int));

        nb[local_iix(x,y,z,*p)][0] = local_iix(x+1, y, z, *p);
        nb[local_iix(x,y,z,*p)][1] = local_iix(x-1, y, z, *p);
        nb[local_iix(x,y,z,*p)][2] = local_iix(x, y+1, z, *p);
        nb[local_iix(x,y,z,*p)][3] = local_iix(x, y-1, z, *p);
        nb[local_iix(x,y,z,*p)][4] = local_iix(x, y, z+1, *p);
        nb[local_iix(x,y,z,*p)][5] = local_iix(x, y, z-1, *p);

        nb[local_iix(x, y, z, *p)][DIR_02] = local_iix(x+1, y+1, z, *p);
        nb[local_iix(x, y, z, *p)][DIR_04] = local_iix(x+1, y, z+1, *p);
        nb[local_iix(x, y, z, *p)][DIR_24] = local_iix(x, y+1, z+1, *p);
        nb[local_iix(x, y, z, *p)][DIR_13] = local_iix(x-1, y-1, z, *p);
        nb[local_iix(x, y, z, *p)][DIR_15] = local_iix(x-1, y, z-1, *p);
        nb[local_iix(x, y, z, *p)][DIR_35] = local_iix(x, y-1, z-1, *p);

        nb[local_iix(x, y, z, *p)][DIR_024] = local_iix(x+1, y+1, z+1, *p);
        nb[local_iix(x, y, z, *p)][DIR_135] = local_iix(x-1, y-1, z-1, *p);

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
  MPI_Request request;
  MPI_Status stat;

  //  fprintf(stderr, "Haloing field\n");

  /* <<
   *
   * SEND LEFT
   */

  // First do the inner (slicey-2) bits
  MPI_Sendrecv(&field[p.inner_xM],
	       (p.slicey-2)*p.Lz, MPI_DOUBLE, p.rank_xM, 0,
	       &field[p.offset_xP],
	       (p.slicey-2)*p.Lz, MPI_DOUBLE, p.rank_xP, 0,
	       MPI_COMM_WORLD, &stat);

  // Then the bottom left corner gets haloed to the bottom endpost
  MPI_Sendrecv(&field[p.inner_xMyM],
	       p.Lz, MPI_DOUBLE, p.rank_xM, 1,
	       &field[p.offset_xPyeM],
	       p.Lz, MPI_DOUBLE, p.rank_xP, 1,
	       MPI_COMM_WORLD, &stat);

  // Then the top left corner gets haloed to the top endpost
  MPI_Sendrecv(&field[p.inner_xMyP],
	       p.Lz, MPI_DOUBLE, p.rank_xM, 2,
	       &field[p.offset_xPyeP],
	       p.Lz, MPI_DOUBLE, p.rank_xP, 2,
	       MPI_COMM_WORLD, &stat);


  /* >>
   *
   * SEND RIGHT
   */


  // First the inner (slicey-2) bits
  MPI_Sendrecv(&field[p.inner_xP],
	       (p.slicey-2)*p.Lz, MPI_DOUBLE, p.rank_xP, 10,
	       &field[p.offset_xM],
	       (p.slicey-2)*p.Lz, MPI_DOUBLE, p.rank_xM, 10,
	       MPI_COMM_WORLD, &stat);

  // Then the bottom right corner gets haloed to the bottom endpost
  MPI_Sendrecv(&field[p.inner_xPyM],
	       p.Lz, MPI_DOUBLE, p.rank_xP, 11,
	       &field[p.offset_xMyeM],
	       p.Lz, MPI_DOUBLE, p.rank_xM, 11,
	       MPI_COMM_WORLD, &stat);


  // Then the top right corner gets haloed to the top endpost
  MPI_Sendrecv(&field[p.inner_xPyP],
	       p.Lz, MPI_DOUBLE, p.rank_xP, 12,
	       &field[p.offset_xMyeP],
	       p.Lz, MPI_DOUBLE, p.rank_xM, 12,
	       MPI_COMM_WORLD, &stat);

  /* /\ /\
   *
   * SEND UP
   */

  // First the inner (slicex-2) bits
  MPI_Sendrecv(&field[p.inner_yP],
	       (p.slicex-2)*p.Lz, MPI_DOUBLE, p.rank_yP, 100,
	       &field[p.offset_yM],
	       (p.slicex-2)*p.Lz, MPI_DOUBLE, p.rank_yM, 100,
	       MPI_COMM_WORLD, &stat);

  // Then the top left corner goes to the left endpost
  MPI_Sendrecv(&field[p.inner_xMyP],
	       p.Lz, MPI_DOUBLE, p.rank_yP, 101,
	       &field[p.offset_yMxeM],
	       p.Lz, MPI_DOUBLE, p.rank_yM, 101,
	       MPI_COMM_WORLD, &stat);


  // Then the top right corner gets haloed to the right endpost
  MPI_Sendrecv(&field[p.inner_xPyP],
	       p.Lz, MPI_DOUBLE, p.rank_yP, 102,
	       &field[p.offset_yMxeP],
	       p.Lz, MPI_DOUBLE, p.rank_yM, 102,
	       MPI_COMM_WORLD, &stat);


  /* \/ \/
   *
   * SEND DOWN
   */

  // First the inner (slicex-2) bits
  MPI_Sendrecv(&field[p.inner_yM],
	       (p.slicex-2)*p.Lz, MPI_DOUBLE, p.rank_yM, 1000,
	       &field[p.offset_yP],
	       (p.slicex-2)*p.Lz, MPI_DOUBLE, p.rank_yP, 1000,
	       MPI_COMM_WORLD, &stat);

  // Then the top left corner goes to the left endpost
  MPI_Sendrecv(&field[p.inner_xMyM],
	       p.Lz, MPI_DOUBLE, p.rank_yM, 1001,
	       &field[p.offset_yPxeM],
	       p.Lz, MPI_DOUBLE, p.rank_yP, 1001,
	       MPI_COMM_WORLD, &stat);


  // Then the top right corner gets haloed to the right endpost
  MPI_Sendrecv(&field[p.inner_xPyM],
	       p.Lz, MPI_DOUBLE, p.rank_yM, 1002,
	       &field[p.offset_yPxeP],
	       p.Lz, MPI_DOUBLE, p.rank_yP, 1002,
	       MPI_COMM_WORLD, &stat);


  /* [ ']
   *
   * SEND UP AND RIGHT
   */
  MPI_Sendrecv(&field[p.inner_xPyP],
	       p.Lz, MPI_DOUBLE, p.rank_xPyP, 10001,
	       &field[p.offset_xMyM],
	       p.Lz, MPI_DOUBLE, p.rank_xMyM, 10001,
	       MPI_COMM_WORLD, &stat);


  /* [, ]
   *
   * SEND DOWN AND LEFT
   */
  MPI_Sendrecv(&field[p.inner_xMyM],
	       p.Lz, MPI_DOUBLE, p.rank_xMyM, 10002,
	       &field[p.offset_xPyP],
	       p.Lz, MPI_DOUBLE, p.rank_xPyP, 10002,
	       MPI_COMM_WORLD, &stat);


  /* [' ]
   *
   * SEND UP AND LEFT
   */
  MPI_Sendrecv(&field[p.inner_xMyP],
	       p.Lz, MPI_DOUBLE, p.rank_xMyP, 10003,
	       &field[p.offset_xPyM],
	       p.Lz, MPI_DOUBLE, p.rank_xPyM, 10003,
	       MPI_COMM_WORLD, &stat);


  /* [ ,]
   *
   * SEND DOWN AND RIGHT
   */
  MPI_Sendrecv(&field[p.inner_xPyM],
	       p.Lz, MPI_DOUBLE, p.rank_xPyM, 10004,
	       &field[p.offset_xMyP],
	       p.Lz, MPI_DOUBLE, p.rank_xMyP, 10004,
	       MPI_COMM_WORLD, &stat);

}

    
int get_z(int n, hydro_params p) {

  return p.inverse[n][2];
}

int get_y(int n, hydro_params p) {
  return p.inverse[n][1];
}

int get_x(int n, hydro_params p) {

  //  fprintf(stderr,"Requesting %d\n", n);
  return p.inverse[n][0];
}

double reduce_sum(double result, hydro_params p) {
  double total = 0.0;
  MPI_Allreduce(&result,&total,1,MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  return total;
}

double reduce_max(double result, hydro_params p) {
  double total = 0.0;
  MPI_Allreduce(&result,&total,1,MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  return total;
}

#else // MPI




int iix(int x, int y, int z, hydro_params p) {
  return ((z+p.Lz)%p.Lz)*p.Ly*p.Lx + ((y+p.Ly)%p.Ly)*p.Lx + ((x+p.Lx)%p.Lx);
}

int **init_nb(hydro_params *p) {
  int x, y, z;

  int **nb;

  // set up nb lookup table
  nb = (int **) malloc(p->N*sizeof(int *));
  
  for(x=0; x<p->Lx; x++) {
    for(y=0; y<p->Ly; y++) {
      for(z=0; z<p->Lz; z++) {
    
        nb[iix(x, y, z, *p)] = (int *)malloc(14*sizeof(int));

        nb[iix(x, y, z, *p)][0] = iix(x+1, y, z, *p);
        nb[iix(x, y, z, *p)][1] = iix(x-1, y, z, *p);
        nb[iix(x, y, z, *p)][2] = iix(x, y+1, z, *p);
        nb[iix(x, y, z, *p)][3] = iix(x, y-1, z, *p);
        nb[iix(x, y, z, *p)][4] = iix(x, y, z+1, *p);
        nb[iix(x, y, z, *p)][5] = iix(x, y, z-1, *p);

        /*
         * Composed directions improve performance as soon as the system
         * ceases to fit inside cache. At small volumes (and hence
         * lower dimensionalities) it is up to 20% slower, but
         * we are interested in good performance on large volumes:
         * not likely that each node will be able to fit everything
         * in cache.
         */
      
        nb[iix(x, y, z, *p)][DIR_02] = iix(x+1, y+1, z, *p);
        nb[iix(x, y, z, *p)][DIR_04] = iix(x+1, y, z+1, *p);
        nb[iix(x, y, z, *p)][DIR_24] = iix(x, y+1, z+1, *p);
        nb[iix(x, y, z, *p)][DIR_13] = iix(x-1, y-1, z, *p);
        nb[iix(x, y, z, *p)][DIR_15] = iix(x-1, y, z-1, *p);
        nb[iix(x, y, z, *p)][DIR_35] = iix(x, y-1, z-1, *p);

        nb[iix(x, y, z, *p)][DIR_024] = iix(x+1, y+1, z+1, *p);
        nb[iix(x, y, z, *p)][DIR_135] = iix(x-1, y-1, z-1, *p);
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


void layout(hydro_params *p) {
  fprintf(stderr,"Layout: doing nothing: no parallelism\n");
  p->fieldN = p->N;
}

void halo_field(double *field, hydro_params p) {
  // do nothing
}

double reduce_sum(double result, hydro_params p) {
  // do nothing
  return result;
}

double reduce_max(double result, hydro_params p) {
  // do nothing
  return result;
}

#endif // MPI
