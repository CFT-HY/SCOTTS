/* fft.c
 *
 * Fourier transform (and power spectrum) of a field
 */
#include "hydro.h"


#ifdef FFT


/* fft_field(hydro_fields f, hydro_params p, double ***field, int step)
 *
 * Carries out the FFT of the supplied field, computes the power
 * spectrum, and stores it binned and labelled with the timestep
 * in a file.
 *
 * For the fluid velocity power spectrum, see velps.c
 * For the GW power spetrum, see gw.c
 *
 * Unlike velps.c and gw.c there is no normalisation done here.
 */
void fft_field(hydro_fields f, hydro_params p, double ***field, int step) {

  ptrdiff_t x_thickness, x_start, alloc_local;

  ptrdiff_t n0 = p.Lx;
  ptrdiff_t n1 = p.Ly;
  ptrdiff_t n2 = p.Lz;

  int slab;

  MPI_Status status;

  double kmode, modsq;

  int x, y, z;
  int i;

  

  double *trim_field = (double *)malloc(p.slicex*p.slicey*p.Lz*sizeof(double));

  double start = clock();

  fftw_mpi_init();

  double checken = 0.0;
  for(x=1; x<=p.slicex; x++) {
    for(y=1; y<=p.slicey; y++) {
      for(z=0; z<p.Lz; z++) {

	trim_field[(x-1)*p.slicey*p.Lz + (y-1)*p.Lz + z] = field[x][y][z];
      }
    }
  }


  alloc_local = fftw_mpi_local_size_3d(n0, n1, n2,
				       MPI_COMM_WORLD, &x_thickness, &x_start);

  // slab is the thickness the 'used' nodes will FFT
  // x_thickness is the thickness the current node will FFT
  slab = (int) x_thickness;

  if(((int)x_thickness) == 0) {
    /*
    fprintf(stderr, "Rank %d nothing to do for FFT: dx=0, x_start=%d!\n",
	    p.rank, (int)x_start);
    */
    slab = 1;
  }

  if(((int)x_thickness) < p.Lx/p.size) {
    fprintf(stderr,
	    "Rank %d giving up in FFT: FFTW told me to use a silly layout!\n",
	    p.rank);
    die(-42);
  }



  fftw_complex *in = fftw_alloc_complex(alloc_local);
  fftw_complex *out = fftw_alloc_complex(alloc_local);

  double *slice = (double *)malloc(slab*p.Ly*p.Lz*sizeof(double));

  int nx = p.Lx/p.slicex;
  int ny = p.Ly/p.slicey;

  int ry;



  // Now planning  
  fftw_plan plan = fftw_mpi_plan_dft_3d(p.Lx, p.Ly, p.Lz,
					in, out, MPI_COMM_WORLD,
					FFTW_FORWARD, FFTW_ESTIMATE);

  for(x=0; x<p.Lx; x++) {

    // Check whether we are the target node for this slab
    if((x >= x_start) && ( x < (x_start + x_thickness))) {

      //      fprintf(stderr, "%d: my slice, x=%d\n", p.rank, x);
      for(ry = 0; ry < ny; ry++) {

	if((x/p.slicex == p.myposx) && (ry == p.myposy)) {

	  memcpy(&slice[(x-x_start)*p.Ly*p.Lz + ry*p.slicey*p.Lz],
		 &trim_field[(x-p.shiftx)*p.slicey*p.Lz],
		 p.slicey*p.Lz*sizeof(double));
	  continue;
	}

	MPI_Recv(&slice[(x-x_start)*p.Ly*p.Lz + ry*p.slicey*p.Lz],
		 p.slicey*p.Lz,
		 MPI_DOUBLE,
		 ry*nx + x/p.slicex,
		 x*ny + ry,
		 MPI_COMM_WORLD,
		 &status);

	/*	
	fprintf(stderr, "%d: got pencil with x=%d\n", p.rank,
		x);
	*/
      }

    } else if((x >= p.shiftx) && (x < (p.shiftx + p.slicex))) {

      /*
      fprintf(stderr, "%d: my pencil, x=%d to dest %d/%d\n",
	      p.rank, x, x, slab);
      */

      MPI_Send(&trim_field[(x-p.shiftx)*p.slicey*p.Lz],
	       p.slicey*p.Lz,
	       MPI_DOUBLE,
	       x/slab,
	       x*ny + p.myposy,
	       MPI_COMM_WORLD);
    }

    // then wait until that slice of x is sorted
    
  }
  // Alloc slab array and do communication to get it into place

  printf0(p, "FFT: Done slicing\n");



  for(x=0;x<x_thickness;x++) {
    for(y=0;y<p.Ly;y++) {
      for(z=0;z<p.Lz;z++) {
	
	in[x*p.Ly*p.Lz + y*p.Lz + z][0] = slice[x*p.Ly*p.Lz + y*p.Lz + z];
	in[x*p.Ly*p.Lz + y*p.Lz + z][1] = 0.0;
      }
    }
  }



  
  fftw_mpi_execute_dft(plan, in, out);


  int nbins = minof3_int(p.Lx, p.Ly, p.Lz);
  double mink = 0.0;
  double maxk = 2.0*M_PI;

  double dk = (maxk-mink)/((double)nbins);

  double *bins = (double *)malloc(nbins*sizeof(double));
  int *counts = (int *)malloc(nbins*sizeof(int));

  for(i=0;i<nbins;i++) {
    bins[i] = 0.0;
    counts[i] = 0;
  }

  int whichbin;

  for(x=0;x<x_thickness;x++) {
    for(y=0;y<p.Ly;y++) {
      for(z=0;z<p.Lz;z++) {
        if(((x+((int)x_start))>p.Lx/2) || (y> p.Ly/2) || (z>p.Lz/2))
          continue;

	kmode = sqrt(
		     ((double)((x+x_start)*(x+x_start)))/((double)(p.Lx*p.Lx))
		     + ((double)(y*y))/((double)(p.Ly*p.Ly))
		     + ((double)(z*z))/((double)(p.Lz*p.Lz))
		     )*2.0*M_PI;

	modsq = out[x*p.Ly*p.Lz + y*p.Lz + z][0]
	  *out[x*p.Ly*p.Lz + y*p.Lz + z][0]
	  + out[x*p.Ly*p.Lz + y*p.Lz + z][1]
	  *out[x*p.Ly*p.Lz + y*p.Lz + z][1];

        whichbin = (int)floor(kmode/dk);
	bins[whichbin] += modsq;
	counts[whichbin]++;

      }
    }
  }


  double red_value;
  int red_count;

  for(i=0;i<nbins;i++) {

    red_value = reduce_sum(bins[i], p);
    red_count = reduce_sum_int(counts[i], p);

    bins[i] = red_value;
    counts[i] = red_count;
  }


  double thisk = dk/2.0;

  if(!p.rank) {

    char fftdest[200];

    sprintf(fftdest,"fft-%d.txt", step);
    
    FILE *fp = fopen(fftdest,"w");
      

    for(i=0;i<nbins;i++) {

    
      char fftdest[200];
      
 
      fprintf(fp, "%lf %g %d\n",
	      thisk, bins[i], counts[i]);

      thisk = thisk + dk;
      
    }

    fclose(fp);
  }
      
  free(bins);
  free(counts);

  free(slice);

  fftw_destroy_plan(plan);
  
  fftw_free(in);
  fftw_free(out);

  free(trim_field);

  fftw_mpi_cleanup();

}






#endif // FFT
