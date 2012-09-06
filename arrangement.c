/* arrangement.c
 *
 * How the physical geometry is laid out in the system, and (where necessary)
 * communication between nodes (when MPI is used).
 */
#include "hydro.h"



double comms_time;

#ifdef MPI



/*
 * Inspired by the contents of setup_layout_generic.c in Kari's code
 */
void layout(hydro_params *p) {


  static int n_primes = 8;
  static int prime[] = {2,3,5,7,11,13,17,19};

  int n;
  int nfactors[n_primes];
  int i = p->size;

  int nx, ny;

  if(!p->rank)
    fprintf(stderr,"Putting %d sites on %d nodes\n", p->N, p->size);


  // You gave us a very stupid number of nodes to work with
  if( (p->Lx*p->Ly % p->size) > 0 ) {
    if(!p->rank)
      fprintf(stderr,"Number of (%d) nodes not a factor of Lx*Ly!\n", p->size);

    exit(-99);
  }

  // Work out prime factors of the number of nodes
  for (n=0; n<n_primes; n++) {
    nfactors[n] = 0;
    while (i % prime[n] == 0) {
      i /= prime[n];
      nfactors[n]++;
    }
  }

  // You gave us a stupid number of nodes to work with
  if(i != 1) {
    if(!p->rank)
      fprintf(stderr, "Can\'t figure out how to put %dx%d%d on %d nodes\n",
	      p->Lx, p->Ly, p->Lz, p->size);

    exit(-100);
  }

  // Show factorisation
  if(!p->rank) {
    for(n=0;n<n_primes;n++)
      fprintf(stderr,"nfactors %d = %d\n", prime[n], nfactors[n]);
  }


  p->slicex = p->Lx;
  p->slicey = p->Ly;

  // Keep cutting the biggest direction until we use up all the prime factors
  for (n=n_primes-1; n>=0; n--) {

    for(i=0; i<nfactors[n]; i++) {
      
      // Longest direction that we can still chop...
      if((p->slicex > p->slicey) && (p->slicex % prime[n] == 0)) {

	if(!p->rank)
	  fprintf(stderr,"Dividing x (%d) by %d\n", p->slicex, prime[n]);

	p->slicex = p->slicex/prime[n];
      } else if(p->slicey % prime[n] == 0) {

	if(!p->rank)
	  fprintf(stderr,"Dividing y (%d) by %d\n", p->slicey, prime[n]);

	p->slicey = p->slicey/prime[n];
      } else {
	// Oops!
	if(!p->rank &&  (p->slicex % prime[n] != 0)) {
	  fprintf(stderr,"Should not get here -- still can\'t divide by x!\n");
	  exit(100);
	} else {
	  if(!p->rank) {
	    fprintf(stderr,"Dividing x direction...\n");
	  }
	  p->slicex = p->slicex/prime[n];
	}

      }
    }
  }


  // Number of nodes in a given direction
  nx = p->Lx/p->slicex;
  ny = p->Ly/p->slicey;

  // Report our findings
  if(!p->rank) {
    fprintf(stderr,"slicing: x=%d y=%d\n", p->slicex, p->slicey);
    fprintf(stderr,"meaning %d nodes in x, %d nodes in y\n", nx, ny);
  }  


  // sites on a node = number of sites + halos + double halos
  p->fieldN = ((p->slicex+2)*(p->slicey+2)*p->Lz);

  if(!p->rank)
    fprintf(stderr,"Each node will need %d = %d + %d entries to store halo\n",
	    p->fieldN,
	    p->N,
	    p->Lz*(2*(p->slicex + p->slicey) + 4));


  // If we knew it worked, then MPI_Getcart might be useful here

  // Where am I
  p->myposx = p->rank % nx;
  p->myposy = (p->rank - p->myposx)/nx;

  // Who are my neighbours (having a full lookup table seems silly)
  p->rank_xM = p->myposy*nx + ((p->myposx-1+nx)%nx);
  p->rank_xP = p->myposy*nx + ((p->myposx+1+nx)%nx);
  p->rank_yM = ((p->myposy - 1 + ny)%ny)*nx + p->myposx;
  p->rank_yP = ((p->myposy + 1 + ny)%ny)*nx + p->myposx;

  p->rank_xMyM = ((p->myposy - 1 + ny)%ny)*nx + ((p->myposx-1+nx)%nx);
  p->rank_xMyP = ((p->myposy + 1 + ny)%ny)*nx + ((p->myposx-1+nx)%nx);
  p->rank_xPyM = ((p->myposy - 1 + ny)%ny)*nx + ((p->myposx+1+nx)%nx);
  p->rank_xPyP = ((p->myposy + 1 + ny)%ny)*nx + ((p->myposx+1+nx)%nx);

  fprintf(stderr, "Node of rank %d, node position is (%d, %d)\n",
	  p-> rank, p->myposx, p->myposy);
  fprintf(stderr,
	  "Node %d: Neighbours are < %d, > %d, /\\ %d, \\/ %d, " \
	  "|_ %d |^ %d ^| %d _| %d\n",
	  p->rank, p->rank_xM, p->rank_xP, p->rank_yM, p->rank_yP,
	  p->rank_xMyM, p->rank_xMyP, p->rank_xPyM, p->rank_xPyP);



  p->shiftx = (p->slicex)*(p->myposx);
  p->shifty = (p->slicey)*(p->myposy);
}



/* void halo_field(double *field, hydro_params p)
 *
 * Do halo communication of the variable in field to neighbours.
 */
void halo_field(double ***field, hydro_params p) {

  MPI_Request request;
  MPI_Status stat;


  clock_t start, end;

  MPI_Barrier(MPI_COMM_WORLD);

  start = clock();

  int mpi_counter = 1;
  int x, y, z;


  /* <<
   *
   * SEND LEFT
   */

  for(y = 1; y <= p.slicey; y++) {

    MPI_Sendrecv(&field[1][y][0],
		 p.Lz, MPI_DOUBLE, p.rank_xM, mpi_counter,
		 &field[p.slicex+1][y][0],
		 p.Lz, MPI_DOUBLE, p.rank_xP, mpi_counter,
		 MPI_COMM_WORLD, &stat);
    
    mpi_counter++;
  }



  /* >>
   *
   * SEND RIGHT
   */

  for(y = 1; y <= p.slicey; y++) {

    MPI_Sendrecv(&field[p.slicex][y][0],
		 p.Lz, MPI_DOUBLE, p.rank_xP, mpi_counter,
		 &field[0][y][0],
		 p.Lz, MPI_DOUBLE, p.rank_xM, mpi_counter,
		 MPI_COMM_WORLD, &stat);
    
    mpi_counter++;
    
  }




  /* /\ /\
   *
   * SEND UP
   */

  for(x = 1; x <= p.slicex; x++) {

      MPI_Sendrecv(&field[x][p.slicey][0],
		   p.Lz, MPI_DOUBLE, p.rank_yP, mpi_counter,
		   &field[x][0][0],
		   p.Lz, MPI_DOUBLE, p.rank_yM, mpi_counter,
		   MPI_COMM_WORLD, &stat);

      mpi_counter++;
  }


  /* \/ \/
   *
   * SEND DOWN
   */


  for(x = 1; x <= p.slicex; x++) {


      MPI_Sendrecv(&field[x][1][0],
		   p.Lz, MPI_DOUBLE, p.rank_yM, mpi_counter,
		   &field[x][p.slicey+1][0],
		   p.Lz, MPI_DOUBLE, p.rank_yP, mpi_counter,
		   MPI_COMM_WORLD, &stat);


      mpi_counter++;
  }


  /* [ ']
   *
   * SEND UP AND RIGHT
   */



  MPI_Sendrecv(&field[p.slicex][p.slicey][0],
	       p.Lz, MPI_DOUBLE, p.rank_xPyP, mpi_counter,
	       &field[0][0][0],
	       p.Lz, MPI_DOUBLE, p.rank_xMyM, mpi_counter,
	       MPI_COMM_WORLD, &stat);
  
  mpi_counter++;

  /* [, ]
   *
   * SEND DOWN AND LEFT
   */



  MPI_Sendrecv(&field[1][1][0],
	       p.Lz, MPI_DOUBLE, p.rank_xMyM, mpi_counter,
	       &field[p.slicex+1][p.slicey+1][0],
	       p.Lz, MPI_DOUBLE, p.rank_xPyP, mpi_counter,
	       MPI_COMM_WORLD, &stat);
  
  mpi_counter++;
  

  /* [' ]
   *
   * SEND UP AND LEFT
   */

      

  MPI_Sendrecv(&field[1][p.slicey][0],
	       p.Lz, MPI_DOUBLE, p.rank_xMyP, mpi_counter,
	       &field[p.slicex+1][0][0],
	       p.Lz, MPI_DOUBLE, p.rank_xPyM, mpi_counter,
	       MPI_COMM_WORLD, &stat);
  
  mpi_counter++;
  

  /* [ ,]
   *
   * SEND DOWN AND RIGHT
   */

  MPI_Sendrecv(&field[p.slicex][1][0],
	       p.Lz, MPI_DOUBLE, p.rank_xPyM, mpi_counter,
	       &field[0][p.slicey+1][0],
	       p.Lz, MPI_DOUBLE, p.rank_xMyP, mpi_counter,
	       MPI_COMM_WORLD, &stat);
  
  mpi_counter++;


  end = clock();

  comms_time += ((double) (end - start)) / CLOCKS_PER_SEC;
}

    

double reduce_sum(double result, hydro_params p) {

  double total = 0.0;
  MPI_Allreduce(&result, &total, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  return total;

}

int reduce_sum_int(int result, hydro_params p) {

  int total = 0;
  MPI_Allreduce(&result, &total, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  return total;

}


double reduce_max(double result, hydro_params p) {

  double total = 0.0;
  MPI_Allreduce(&result, &total, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  return total;

}

double reduce_min(double result, hydro_params p) {

  double total = 0.0;
  MPI_Allreduce(&result, &total, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
  return total;

}


int reduce_and(int result, hydro_params p) {

  int total = 0;
  MPI_Allreduce(&result, &total, 1, MPI_INT, MPI_LAND, MPI_COMM_WORLD);
  return total;

}


void init_comms_time(hydro_params *p) {
  comms_time = 0.0;
}

double get_comms_time(hydro_params *p) {
  return comms_time;
}



#else // MPI



void layout(hydro_params *p) {
  fprintf(stderr,"Layout: doing nothing: no parallelism\n");
  p->fieldN = (p->Lx+2)*(p->Ly+2)*(p->Lz);
  p->slicex = p->Lx;
  p->slicey = p->Ly;
  p->shiftx = 0;
  p->shifty = 0;
}


void halo_field(double ***field, hydro_params p) {


  clock_t start, end;

  start = clock();

  /* <<
   *
   * SEND LEFT
   */

  int mpi_counter = 1;
  int x, y, z;


  for(y = 1; y <= p.Ly; y++) {
    for(z = 0; z < p.Lz; z++) {

      field[p.Lx+1][y][z] = field[1][y][z];

      mpi_counter++;
    }
  }


  /* >>
   *
   * SEND RIGHT
   */

  for(y = 1; y <= p.Ly; y++) {
    for(z = 0; z < p.Lz; z++) {


      field[0][y][z] = field[p.Lx][y][z];

      mpi_counter++;
    }
  }




  /* /\ /\
   *
   * SEND UP
   */

  for(x = 1; x <= p.Lx; x++) {
    for(z = 0; z < p.Lz; z++) {


      field[x][0][z] = field[x][p.Ly][z];

      mpi_counter++;
    }
  }



  /* \/ \/
   *
   * SEND DOWN
   */


  for(x = 1; x <= p.Lx; x++) {
    for(z = 0; z < p.Lz; z++) {


      field[x][p.Lx+1][z] = field[x][1][z];

      mpi_counter++;
    }
  }



  /* [ ']
   *
   * SEND UP AND RIGHT
   */



  for(z = 0; z < p.Lz; z++) {
      
    field[0][0][z] = field[p.Lx][p.Ly][z];

    mpi_counter++;
  }

  /* [, ]
   *
   * SEND DOWN AND LEFT
   */


  for(z = 0; z < p.Lz; z++) {
      
    field[p.Lx+1][p.Ly+1][z] = field[1][1][z];

    mpi_counter++;
  }

  /* [' ]
   *
   * SEND UP AND LEFT
   */

  for(z = 0; z < p.Lz; z++) {
      
    field[p.Lx+1][0][z] = field[1][p.Ly][z];

    mpi_counter++;
  }


  /* [ ,]
   *
   * SEND DOWN AND RIGHT
   */

  for(z = 0; z < p.Lz; z++) {
      
    field[0][p.Ly+1][z] = field[p.Lx][1][z];

    mpi_counter++;
  }


  end = clock();

  comms_time += ((double) (end - start)) / CLOCKS_PER_SEC;
  
}

double reduce_sum(double result, hydro_params p) {
  // do nothing
  return result;
}

int reduce_sum_int(int result, hydro_params p) {
  // Do nothing
  return result;
}


double reduce_max(double result, hydro_params p) {
  // do nothing
  return result;
}

double reduce_min(double result, hydro_params p) {
  // do nothing
  return result;
}

int reduce_and(int result, hydro_params p) {

  // do nothing
  return result;

}


void init_comms_time(hydro_params *p) {
  comms_time = 0.0;
}

double get_comms_time(hydro_params *p) {
  return comms_time;
}


#endif // MPI
