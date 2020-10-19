/** @file scalarps.c
 *
 * Calculates power spectrum of a scalar field.
 */
#include "hydro.h"


#ifdef FFT


/** Compute and bin the power spectrum of the supplied scalar field in
 * momentum space, and write the power spectrum to a file.
 *
 * Computes the power spectrum of supplied field in fourier space, and
 * stores it binned and labelled with the timestep in a file.
 *
 * See fft_func.c for the `fft_field` functions to transform fields
 * into fourier space. The momentum space fields are decomposed in pencils.
 *
 * For the calculation of vector field power spectra, see vectorps.c
 * For tensor power spetra, see tensorps.c
 *
 */
void scalarps(hydro_params p, fftwf_complex *field, int step, char *label) {

  ptrdiff_t x_thickness, x_start, alloc_local;

  ptrdiff_t n0 = p.Lx;
  ptrdiff_t n1 = p.Ly;
  ptrdiff_t n2 = p.Lz;

  //  int slab;

  MPI_Status status;

  float kmode, modsq;

  int x, y, z;
  int i;

  


  float start = clock();

  alloc_local = fftwf_mpi_local_size_3d(n0, n1, n2,
					MPI_COMM_WORLD, &x_thickness, &x_start);

  // Now we perform binning
  float fft_norm = (1.0/(((float)p.Lx)*((float)p.Ly)*((float)p.Lz)));

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

	modsq = field[x*p.Ly*p.Lz + y*p.Lz + z][0]
	  *field[x*p.Ly*p.Lz + y*p.Lz + z][0]
	  + field[x*p.Ly*p.Lz + y*p.Lz + z][1]
	  *field[x*p.Ly*p.Lz + y*p.Lz + z][1];

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

    // Do normalisation here because we don't do it during fft_scalar...
    bins[i] = red_value/(fft_norm*fft_norm);
    counts[i] = red_count;
  }


  float thisk = dk/2.0;

  if(!p.rank) {

    char fftdest[200];
    if(label != NULL){
      if(*label){
        sprintf(fftdest,"%s-ps-step%d.txt",label,step);
      }
    }
    else{
      sprintf(fftdest,"ps-step%d.txt",step);
    }
    
    FILE *fp = fopen(fftdest,"w");
      

    for(i=0;i<nbins;i++) {

    
      char fftdest[200];
      
      fprintf(fp, "%lf %g %d\n",
	      thisk/(p.a*p.dx), bins[i], counts[i]);

      thisk = thisk + dk;
      
    }

    fclose(fp);
  }

  free(bins);
  free(counts);
  
  float end = clock();
  if(label != NULL){
    if(*label){
      printf0(p, "%s scalar power spectrum calculation took %lf\n", label,
        ((float) (end - start)) / CLOCKS_PER_SEC);
    }
  }
  else{
    printf0(p, "scalar power spectrum calculation took %lf\n",
      ((float) (end - start)) / CLOCKS_PER_SEC);
  }
    
}



#endif // FFT
