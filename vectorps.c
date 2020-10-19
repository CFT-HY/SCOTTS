/** @file vectorps.c
 *
 * File for binning vector power spectras. Computes total ETC,
 * longitudinal ETC, and rotational ETC.
 *
 * FFT functions in fft_func.c
 * Projectors are in projectors.c
 * Tensor power spectra are in tensorps.c,
 * Scalar power spectrum calculation is in scalarps.c.
 * Unequal time correlators are in uetc.c
 */
#include "hydro.h"


#ifdef FFT

/** Combine the data from each slice (supplied separately by each
 * rank to this function), bin, compute the histogram and output.
 */
void histogram(hydro_params p, float *slice, char *filename,
	       int slab, int x_start) {
  // The rest proceeds as in gw.c

  int i, x, y, z;

  float kmode;

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


  for(x=0;x<slab;x++) {
    for(y=0;y<p.Ly;y++) {
      for(z=0;z<p.Lz;z++) {

	/*
        if(((x+x_start)>p.Lx/2) || (y> p.Ly/2) || (z>p.Lz/2))
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

	whichbin = (int)floor(kmode/dk);
	bins[whichbin] += slice[x*p.Ly*p.Lz + y*p.Lz + z];
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

  // spokesman does the final analysis
  if(!p.rank) {


  FILE *fp = fopen(filename,"w");

    for(i=0;i<nbins;i++) {


      fprintf(fp, "%lf %g %d\n",
	      thisk/(p.a*p.dx), ((float)(i+1))*bins[i], counts[i]);

      thisk = thisk + dk;
    }

  fclose(fp);

  }

  free(bins);
  free(counts);
}





/** Take the supplied momentum space vector field, apply projectors to
 * split it into rotational, and longitudinal components, compute the
 * power spectrum of said components and write these and the total
 * power spectrum to file.
 *
 * See projectors.c for the projectors used here.
 *
 * See fft_func.c for the `fft_field` functions to transform fields
 * into fourier space. The momentum space fields are decomposed in pencils.
 *
 * For the calculation of scalar field power spectra, see scalarps.c
 * For tensor power spetra, see tensorps.c
 */

void vectorps(hydro_params p, fftwf_complex **outcpts, int step, char *label) {

  ptrdiff_t x_thickness, x_start, alloc_local;

  ptrdiff_t n0 = p.Lx;
  ptrdiff_t n1 = p.Ly;
  ptrdiff_t n2 = p.Lz;

  //  int slab;

  MPI_Status status;

  float modsq;

  int x, y, z;
  int i;

  if(label != NULL){
    if(*label){
      printf0(p, "Starting %s vector power spectrum calculation.\n",label);
    }
  }
  else{
   printf0(p, "Starting vector power spectrum calculation.\n");
  }
  float start = clock();

  alloc_local = fftwf_mpi_local_size_3d(n0, n1, n2,
				       MPI_COMM_WORLD, &x_thickness, &x_start);

  /*
   * At this point, we should have a normed FFT
   * of each component in outcpts[i][k]
   */

  float *slice_rot = (float *)malloc(x_thickness*p.Ly*p.Lz*sizeof(float));
  float *slice_div = (float *)malloc(x_thickness*p.Ly*p.Lz*sizeof(float));
  float *slice_tot = (float *)malloc(x_thickness*p.Ly*p.Lz*sizeof(float));

  // The brains of the operation:
  // Turn the FFT'd vector into projected power spectrum
  split_vector(p, x_start, x_thickness, slice_rot, slice_div, slice_tot, outcpts);



  char fftdest[200];
  if(label != NULL){
    if(*label){
      sprintf(fftdest,"%s-rot-ps-step%d.txt", label,step);
    }
  }
  else{
    sprintf(fftdest,"rot-ps-step%d.txt",step);
  }
  histogram(p, slice_rot, fftdest, x_thickness, x_start);
  if(label != NULL){
    if(*label){
      sprintf(fftdest,"%s-div-ps-step%d.txt", label,step);
    }
  }
  else{
    sprintf(fftdest,"div-ps-step%d.txt",step);
  }
  histogram(p, slice_div, fftdest, x_thickness, x_start);
  if(label != NULL){
    if(*label){
      sprintf(fftdest,"%s-tot-ps-step%d.txt", label,step);
    }
  }
  else{
    sprintf(fftdest,"tot-ps-step%d.txt",step);
  }
  histogram(p, slice_tot, fftdest, x_thickness, x_start);


  // Tidy up

  free(slice_rot);
  free(slice_div);
  free(slice_tot);


  float end = clock();
  if(label != NULL){
    if(*label){
      printf0(p, "%s vector power spectrum calculation took %lf\n", label,
        ((float) (end - start)) / CLOCKS_PER_SEC);
    }
  }
  else{
    printf0(p, "vector power spectrum calculation took %lf\n",
      ((float) (end - start)) / CLOCKS_PER_SEC);
  }
}

#endif // FFT
