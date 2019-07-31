/** @file fft.c
 *
 * Fourier transform (and power spectrum) of a field.
 */
#include "hydro.h"


#ifdef FFT


/** Perform the FFT of a supplied field, and write the power spectrum
 * to a file.
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
void fft_field(float ***field, hydro_params p, int step, char *label) {

  ptrdiff_t x_thickness, x_start, alloc_local;

  ptrdiff_t n0 = p.Lx;
  ptrdiff_t n1 = p.Ly;
  ptrdiff_t n2 = p.Lz;

  //  int slab;

  MPI_Status status;

  float kmode, modsq;

  int x, y, z;
  int i;

  

  float *trim_field = (float *)malloc(p.slicex*p.slicey*p.Lz*sizeof(float));

  float start = clock();

  fftwf_mpi_init();

  float checken = 0.0;
  for(x=1; x<=p.slicex; x++) {
    for(y=1; y<=p.slicey; y++) {
      for(z=0; z<p.Lz; z++) {

	trim_field[(x-1)*p.slicey*p.Lz + (y-1)*p.Lz + z] = field[x][y][z];
      }
    }
  }


  alloc_local = fftwf_mpi_local_size_3d(n0, n1, n2,
				       MPI_COMM_WORLD, &x_thickness, &x_start);

  int *thicknesses = malloc(p.size*sizeof(int));
  int *starts = malloc(p.size*sizeof(int));
  int *map = malloc(p.Lx*sizeof(int));
			    

  if(p.rank) {
    MPI_Send(&x_thickness, 1, MPI_INTEGER, 0, p.rank, MPI_COMM_WORLD);
    MPI_Send(&x_start, 1, MPI_INTEGER, 0, p.rank, MPI_COMM_WORLD);

  } else {
    thicknesses[0] = x_thickness;
    starts[0] = x_start;

    for(i=1; i<p.size; i++) {
      MPI_Recv(&thicknesses[i], 1, MPI_INTEGER, i, i, MPI_COMM_WORLD, &status);
      MPI_Recv(&starts[i], 1, MPI_INTEGER, i, i, MPI_COMM_WORLD, &status);

      //      fprintf(stderr,"rank %d, thickness %d, start %d\n", i, thicknesses[i], starts[i]);
    }


    for(x=0; x<p.Lx; x++) {
      map[x] = -1;
      for(i=0; i<p.size; i++) {
	if(starts[i] <= x && (starts[i] + thicknesses[i]) > x) {
	  map[x] = i;
	  //	  fprintf(stderr, "therefore node %d is responsible for x=%d\n", i, x);
	  break;
	}
      }
      if(map[x] == -1) {
	fprintf(stderr,"cannot find a node responsible for x=%d!\n", x);
	die(-99);
      }

    }
  }
    
  printf0(p,"broadcasting our results\n");
  MPI_Bcast(map, p.Lx, MPI_INTEGER, 0, MPI_COMM_WORLD);
  printf0(p,"now to layout\n");

  /*
  // slab is the thickness the 'used' nodes will FFT
  // x_thickness is the thickness the current node will FFT
  slab = (int) x_thickness;

  int maxslab = reduce_max(slab, p);

  if(((int)x_thickness) == 0) {
    fprintf(stderr, "Rank %d nothing to do for FFT: dx=0, x_start=%d! maxslab=%d\n",
	    p.rank, (int)x_start, maxslab);
    slab = maxslab;
  }
  /* else if(((int)x_thickness) < p.Lx/p.size) {
    fprintf(stderr,
	    "Rank %d giving up in FFT: FFTW told me to use a silly layout!\n",
	    p.rank);
    die(-42);
  }
  */



  fftwf_complex *in = fftwf_alloc_complex(alloc_local);
  fftwf_complex *out = fftwf_alloc_complex(alloc_local);

  float *slice = (float *)malloc(x_thickness*p.Ly*p.Lz*sizeof(float));

  int nx = p.Lx/p.slicex;
  int ny = p.Ly/p.slicey;

  int ry;



  // Now planning  
  fftwf_plan plan = fftwf_mpi_plan_dft_3d(p.Lx, p.Ly, p.Lz,
					in, out, MPI_COMM_WORLD,
					FFTW_FORWARD, FFTW_ESTIMATE);

  for(x=0; x<p.Lx; x++) {

    // Check whether we are the target node for this slab
    if(map[x] == p.rank) {

      //      fprintf(stderr, "%d: my slice, x=%d\n", p.rank, x);
      for(ry = 0; ry < ny; ry++) {

	if((x/p.slicex == p.myposx) && (ry == p.myposy)) {
	  memcpy(&slice[(x-x_start)*p.Ly*p.Lz + ry*p.slicey*p.Lz],
		 &trim_field[(x-p.shiftx)*p.slicey*p.Lz],
		 p.slicey*p.Lz*sizeof(float));
	  continue;
	}

	MPI_Recv(&slice[(x-x_start)*p.Ly*p.Lz + ry*p.slicey*p.Lz],
		 p.slicey*p.Lz,
		 MPI_FLOAT,
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
	      p.rank, x, x, map[x]);
      */

      MPI_Send(&trim_field[(x-p.shiftx)*p.slicey*p.Lz],
	       p.slicey*p.Lz,
	       MPI_FLOAT,
	       map[x],
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



  
  fftwf_mpi_execute_dft(plan, in, out);


  int nbins = minof3_int(p.Lx, p.Ly, p.Lz);
  float mink = 0.0;
  float maxk = 2.0*M_PI;

  float dk = (maxk-mink)/((float)nbins);

  float *bins = (float *)malloc(nbins*sizeof(float));
  int *counts = (int *)malloc(nbins*sizeof(int));

  for(i=0;i<nbins;i++) {
    bins[i] = 0.0;
    counts[i] = 0;
  }

  int whichbin;

  int true_x, true_y, true_z;

  for(x=0;x<x_thickness;x++) {
    for(y=0;y<p.Ly;y++) {
      for(z=0;z<p.Lz;z++) {

	/*
        if(((x+((int)x_start))>p.Lx/2) || (y> p.Ly/2) || (z>p.Lz/2))
          continue;
	*/


        if(x+x_start > p.Lx/2)
          true_x = p.Lx - (x+x_start);
        else
          true_x = x+x_start;

        if(y > p.Ly/2)
          true_y = p.Ly - y;
        else
          true_y = y;

        if(z > p.Lz/2)
          true_z = p.Lz - z;
	else
          true_z = z;

        kmode = sqrt(
                     ((float)(true_x*true_x))/((float)(p.Lx*p.Lx))
                     + ((float)(true_y*true_y))/((float)(p.Ly*p.Ly))
                     + ((float)(true_z*true_z))/((float)(p.Lz*p.Lz))
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


  float red_value;
  int red_count;

  for(i=0;i<nbins;i++) {

    red_value = reduce_sum(bins[i], p);
    red_count = reduce_sum_int(counts[i], p);

    bins[i] = red_value;
    counts[i] = red_count;
  }


  float thisk = dk/2.0;

  if(!p.rank) {

    char fftdest[200];

    sprintf(fftdest,"fft-%s-%d.txt",label,step);
    
    FILE *fp = fopen(fftdest,"w");
      

    for(i=0;i<nbins;i++) {

    
      char fftdest[200];
      
 
      fprintf(fp, "%lf %g %d\n",
	      thisk, bins[i], counts[i]);

      thisk = thisk + dk;
      
    }

    fclose(fp);
  }

  printf0(p, "Basic FFT done\n");
      
  free(bins);
  free(counts);

  free(slice);

  fftwf_destroy_plan(plan);
  
  fftwf_free(in);
  fftwf_free(out);

  free(trim_field);

  fftwf_mpi_cleanup();

}

#ifndef SCALAR

/** Create internal energy field e = E/W and then perform fft of this field.
 *
 * Not a fan of this implementation, must be something neater/other 
 * place for this to go.
 */
void fft_e(hydro_fields f, hydro_params p, int step, char *label){
  int x, y, z;
  float ***e = make_field(p);
  for(x = 1; x <= p.slicex; x++) {
    for(y = 1; y <= p.slicey; y++) {
      for(z = 0; z < p.Lz; z++) {
	e[x][y][z] = f.E[x][y][z]/f.W[x][y][z];
      }
    }
  }
  halo_field(e, p);

  fft_field(e, p, step, label);

  free_field(p, e);
}

#endif //!SCALAR




#endif // FFT
