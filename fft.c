#include "hydro.h"

#ifdef FFT

void fft(hydro_fields f, hydro_params p) {

  fftw_mpi_init();

  ptrdiff_t x_thickness, x_start, alloc_local;

  ptrdiff_t n0 = p.Lx;
  ptrdiff_t n1 = p.Ly;
  ptrdiff_t n2 = p.Lz;

  int slab;

  MPI_Status status;


  alloc_local = fftw_mpi_local_size_3d(n0, n1, n2, MPI_COMM_WORLD, &x_thickness, &x_start);

  fprintf(stderr,"FFTW told me to make a slab dx=%d thick starting at x=%d\n",
	  (int)x_thickness,
	  (int)x_start);
  fprintf(stderr,"And I will need %d data\n", (int)alloc_local);

  slab = (int)x_thickness;


  if(((int)x_thickness) != p.Lx/p.size) {
    fprintf(stderr,"Can\'t go on like this: FFTW told me to use a silly layout!\n");
    exit(-42);
  }

  if(((int)x_thickness) == 0) {
    fprintf(stderr,"Can\'t go on like this: dx=0!\n");
    exit(-43);
  }

  fftw_complex *in = fftw_alloc_complex(alloc_local);
  fftw_complex *out = fftw_alloc_complex(alloc_local);



  double *slice = (double *)malloc(slab*p.Ly*p.Lz*sizeof(double));

  double *ourdata = (double *)malloc(p.slicex*p.slicey*p.Lz*sizeof(double));

  int x,y, z;

  // Prepare our data for the FFT (stress-energy tensor?)

  // This is test data:
  for(x=0;x<p.slicex;x++) {
    for(y=0;y<p.slicey;y++) {
      for(z=0;z<p.Lz;z++) {
	ourdata[x*p.slicey*p.Lz + y*p.Lz + z] =
	  cos(((double)(x+p.shiftx))*2.0*M_PI/((double)p.Lx))
	  *cos(((double)(y+p.shifty))*2.0*M_PI/((double)p.Ly))
	  *cos(((double)z)*2.0*M_PI/((double)p.Lz));
	//	fprintf(stderr,"ourdata %lf\n", ourdata[x*p.slicey*p.Lz + y*p.Lz + z]);
	// (double)((x+p.shiftx)*10000 + (y+p.shifty)*100 + z);
      }
    }
  }


  int nx = p.Lx/p.slicex;
  int ny = p.Ly/p.slicey;

  int ry;

  for(x=0; x<p.Lx; x++) {

    if( (x >= x_start) && (x<(x_start + x_thickness)) ) {
      for(ry = 0; ry < ny; ry++) {
	fprintf(stderr,"Node %d: I own x=%d, so receiving from %d\n", p.rank, x, ry*nx + x/p.slicex);
	if(x/p.slicex == p.myposx && ry == p.myposy) {
	  fprintf(stderr,"Node %d: self receive\n", p.rank);
	  memcpy(&slice[(x-x_start)*p.Ly*p.Lz + ry*p.slicey*p.Lz],
		 &ourdata[(x-p.shiftx)*p.slicey*p.Lz],
		 p.slicey*p.Lz*sizeof(double));
	  continue;
	}
	fprintf(stderr,"Node %d: expecting tag %d\n", p.rank, x*ny+ry);
	MPI_Recv(&slice[(x-x_start)*p.Ly*p.Lz + ry*p.slicey*p.Lz],
		 p.slicey*p.Lz,
		 MPI_DOUBLE,
		 ry*nx + x/p.slicex,
		 x*ny + ry,
		 MPI_COMM_WORLD,
		 &status);
      }
      fprintf(stderr,"Node %d: Got it!\n", p.rank);
    } else if(x >= p.shiftx && x < (p.shiftx + p.slicex)) {
      //     fprintf(stderr,"slab is %d\n",slab);
      fprintf(stderr, "Node %d: I have some of x=%d, so sending to owner %d\n", p.rank, x, x/slab);
      fprintf(stderr, "Node %d: attaching tag %d\n", p.rank, x*ny + p.myposy);
      MPI_Send(&ourdata[(x-p.shiftx)*p.slicey*p.Lz],
	       p.slicey*p.Lz,
	       MPI_DOUBLE,
	       x/slab,
	       x*ny + p.myposy,
	       MPI_COMM_WORLD);
      fprintf(stderr, "Node %d: Sent it!\n", p.rank);

    } else {
      fprintf(stderr,"x=%d, rank is %d, nothing to do with me\n", x, p.rank);
    }

	MPI_Barrier(MPI_COMM_WORLD);

	if(!p.rank)
	  fprintf(stderr,"Everyone in sync after x=%d\n\n\n\n",x);


    // then wait until that slice of x is sorted
    
  }
  // Alloc slab array and do communication to get it into place

  


  for(x=0;x<slab;x++) {
    for(y=0;y<p.Ly;y++) {
	for(z=0;z<p.Lz;z++) {
	  fprintf(stderr,"SLAB %d %d %d : %g\n", x+x_start, y,  z, slice[x*p.Ly*p.Lz + y*p.Lz + z]);
	  in[x*p.Ly*p.Lz + y*p.Lz + z][0] =  slice[x*p.Ly*p.Lz + y*p.Lz + z];
	  in[x*p.Ly*p.Lz + y*p.Lz + z][1] = 0.0;
	}
    }
  }


  // Now planning
  
  fftw_plan plan = fftw_mpi_plan_dft_3d(p.Lx, p.Ly, p.Lz, in, out, MPI_COMM_WORLD,
					FFTW_FORWARD, FFTW_ESTIMATE);

  
  fftw_mpi_execute_dft(plan, in, out);


  for(x=0;x<slab;x++) {
    for(y=0;y<p.Ly;y++) {
      for(z=0;z<p.Lz;z++) {
	fprintf(stderr,"OUT %d %d %d: %g\n", x+x_start, y,  z,out[x*p.Ly*p.Lz + y*p.Lz + z][0]);
      }
    }
  }

  fftw_destroy_plan(plan);
  

  fftw_free(in);
  fftw_free(out);
  fftw_mpi_cleanup();
}

#endif // FFT
