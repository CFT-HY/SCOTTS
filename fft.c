/* fft.c
 *
 * Fourier transform of a field (not currently used)
 */
#include "hydro.h"


#ifdef FFT



void fft_field(hydro_fields f, hydro_params p) {

  ptrdiff_t x_thickness, x_start, alloc_local;

  ptrdiff_t n0 = p.Lx;
  ptrdiff_t n1 = p.Ly;
  ptrdiff_t n2 = p.Lz;

  int slab;

  MPI_Status status;

  double kmode, modsq;

  int x, y, z;

  double ***en = make_field(p);

  double *trim_en = (double *)malloc(p.slicex*p.slicey*p.Lz*sizeof(double));

  fftw_mpi_init();

  energy_density(f, p, en);

  double checken = 0.0;
  for(x=1; x<=p.slicex; x++) {
    for(y=1; y<=p.slicey; y++) {
      for(z=0; z<p.Lz; z++) {
	checken += en[x][y][z];

	trim_en[(x-1)*p.slicey*p.Lz + (y-1)*p.Lz + z] = en[x][y][z];
      }
    }
  }

  fprintf(stderr,"Check: energy density claimed %lf, total energy %lf\n",
	  reduce_sum(checken, p), reduce_sum(total_energy(f, p), p));


  alloc_local = fftw_mpi_local_size_3d(n0, n1, n2,
				       MPI_COMM_WORLD, &x_thickness, &x_start);


  
#ifdef FFT_DEBUG
  fprintf(stderr, "FFTW told me to make a slab dx=%d thick starting at x=%d\n",
	  (int)x_thickness,
	  (int)x_start);
  fprintf(stderr, "And I will need %d data\n", (int)alloc_local);
#endif // FFT_DEBUG

  slab = (int) x_thickness;

  if(((int)x_thickness) != p.Lx/p.size) {
    fprintf(stderr,
	    "Giving up in FFT: FFTW told me to use a silly layout!\n");
    exit(-42);
  }

  if(((int)x_thickness) == 0) {
    fprintf(stderr, "Giving up in FFT: dx=0!\n");
    exit(-43);
  }

  fftw_complex *in = fftw_alloc_complex(alloc_local);
  fftw_complex *out = fftw_alloc_complex(alloc_local);

  double *slice = (double *)malloc(slab*p.Ly*p.Lz*sizeof(double));


  /*
  // Prepare and make use of test data
  double *ourdata = (double *)malloc(p.slicex*p.slicey*p.Lz*sizeof(double));

  int x,y, z;

  // Prepare our data for the FFT (stress-energy tensor?)

  // This is test data:
  for(x = 0; x < p.slicex; x++) {
    for(y = 0; y < p.slicey; y++) {
      for(z = 0; z < p.Lz; z++) {
	ourdata[x*p.slicey*p.Lz + y*p.Lz + z] =
	  cos(((double)(x+p.shiftx))*2.0*M_PI/((double)p.Lx))
	  *cos(((double)(y+p.shifty))*2.0*M_PI/((double)p.Ly))
	  *cos(((double)z)*2.0*M_PI/((double)p.Lz));
	//	fprintf(stderr, "ourdata %lf\n",
	                        ourdata[x*p.slicey*p.Lz + y*p.Lz + z]);
	// (double)((x+p.shiftx)*10000 + (y+p.shifty)*100 + z);
      }
    }
  }
  */

  int nx = p.Lx/p.slicex;
  int ny = p.Ly/p.slicey;

  int ry;

  for(x=0; x<p.Lx; x++) {

    // Check whether we are the target node for this slab
    if((x >= x_start) && ( x < (x_start + x_thickness))) {

      for(ry = 0; ry < ny; ry++) {

#ifdef FFT_DEBUG
	fprintf(stderr,"Node %d: I own x=%d, so receiving from %d\n",
		p.rank, x, ry*nx + x/p.slicex);
#endif // FFT_DEBUG

	if((x/p.slicex == p.myposx) && (ry == p.myposy)) {

#ifdef FFT_DEBUG
	  fprintf(stderr,"Node %d: self receive\n", p.rank);
#endif // FFT_DEBUG

	  memcpy(&slice[(x-x_start)*p.Ly*p.Lz + ry*p.slicey*p.Lz],
		 &trim_en[(x-p.shiftx)*p.slicey*p.Lz],
		 p.slicey*p.Lz*sizeof(double));
	  continue;
	}

#ifdef FFT_DEBUG
	fprintf(stderr,"Node %d: expecting tag %d\n", p.rank, x*ny+ry);
#endif // FFT_DEBUG

	MPI_Recv(&slice[(x-x_start)*p.Ly*p.Lz + ry*p.slicey*p.Lz],
		 p.slicey*p.Lz,
		 MPI_DOUBLE,
		 ry*nx + x/p.slicex,
		 x*ny + ry,
		 MPI_COMM_WORLD,
		 &status);
      }

#ifdef FFT_DEBUG
      fprintf(stderr,"Node %d: Got it!\n", p.rank);
#endif // FFT_DEBUG

    } else if((x >= p.shiftx) && (x < (p.shiftx + p.slicex))) {

#ifdef FFT_DEBUG
      //     fprintf(stderr,"slab is %d\n",slab);
      fprintf(stderr,
	      "Node %d: I have some of x=%d, so sending to owner %d\n",
	      p.rank, x, x/slab);
      fprintf(stderr, "Node %d: attaching tag %d\n",
	      p.rank, x*ny + p.myposy);
#endif // FFT_DEBUG

      MPI_Send(&trim_en[(x-p.shiftx)*p.slicey*p.Lz],
	       p.slicey*p.Lz,
	       MPI_DOUBLE,
	       x/slab,
	       x*ny + p.myposy,
	       MPI_COMM_WORLD);

#ifdef FFT_DEBUG
      fprintf(stderr, "Node %d: Sent it!\n", p.rank);

    } else {
      fprintf(stderr,"x=%d, rank is %d, nothing to do with me\n", x, p.rank);
#endif // FFT_DEBUG

    }

#ifdef FFT_DEBUG
	MPI_Barrier(MPI_COMM_WORLD);

	if(!p.rank)
	  fprintf(stderr,"Everyone in sync after x=%d\n\n\n\n",x);
#endif // FFT_DEBUG

    // then wait until that slice of x is sorted
    
  }
  // Alloc slab array and do communication to get it into place

  
#ifdef FFT_DEBUG
  fprintf(stderr,"Done slabulation...\n");
#endif // FFT_DEBUG

  checken = 0.0;
  for(x=0; x<slab; x++) {
    for(y=0; y<p.Ly; y++) {
      for(z=0; z<p.Lz; z++) {
	checken += slice[x*p.Ly*p.Lz + y*p.Lz + z];
      }
    }
  }

  fprintf(stderr,"Second check: energy density claimed %lf\n",
	  reduce_sum(checken, p));




  for(x=0;x<slab;x++) {
    for(y=0;y<p.Ly;y++) {
      for(z=0;z<p.Lz;z++) {
	
#ifdef FFT_DEBUG
	fprintf(stderr, "SLAB %d %d %d : %g\n",
		x + x_start, y,  z, slice[x*p.Ly*p.Lz + y*p.Lz + z]);
#endif // FFT_DEBUG
	
	in[x*p.Ly*p.Lz + y*p.Lz + z][0] = slice[x*p.Ly*p.Lz + y*p.Lz + z];
	in[x*p.Ly*p.Lz + y*p.Lz + z][1] = 0.0;
      }
    }
  }


  // Now planning
  
  fftw_plan plan = fftw_mpi_plan_dft_3d(p.Lx, p.Ly, p.Lz,
					in, out, MPI_COMM_WORLD,
					FFTW_FORWARD, FFTW_ESTIMATE);

  
  fftw_mpi_execute_dft(plan, in, out);


  // Temporary until we get the binning sorted out
  char fftdest[200];

  sprintf(fftdest,"fft-%d.txt",p.rank);

  FILE *fp = fopen(fftdest,"w");
 

  for(x=0;x<slab;x++) {
    for(y=0;y<p.Ly;y++) {
      for(z=0;z<p.Lz;z++) {

	//	if((x+x_start)==0 && y==0 && z==0)
	//	  continue;

	kmode = sqrt(
		     ((double)((x+x_start)*(x+x_start)))/((double)(p.Lx*p.Lx))
		     + ((double)(y*y))/((double)(p.Ly*p.Ly))
		     + ((double)(z*z))/((double)(p.Lz*p.Lz))
		     )*2.0*M_PI;

	modsq = out[x*p.Ly*p.Lz + y*p.Lz + z][0]
	  *out[x*p.Ly*p.Lz + y*p.Lz + z][0]
	  + out[x*p.Ly*p.Lz + y*p.Lz + z][1]
	  *out[x*p.Ly*p.Lz + y*p.Lz + z][1];

	fprintf(fp, "%lf %d %d %d %g\n", kmode,
		x + ((int)x_start), y, z, modsq);

      }
    }
  }

  close(fp);

  fftw_destroy_plan(plan);
  
  fftw_free(in);
  fftw_free(out);

  fftw_mpi_cleanup();

  free_field(p, en);
}






#endif // FFT
