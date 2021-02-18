/** @file tensorps.c
 *
 * File for binning tensor power spectras. Computes transverse
 * traceless ETCs.
 *
 * FFT functions in fft_func.c
 * Projectors are in projectors.c
 * Vector power spectra are in vectorps.c,
 * Scalar power spectrum calculation is in scalarps.c.
 * Unequal time correlators are in uetc.c
 */
#include "hydro.h"


#ifdef FFT




/** Calculates tensor power spectrum. Takes a tensor in momentum
 * space, and then projects it to the transverse traceless gauge and
 * computes the power spectrum.
 *
 */
float tensorps(hydro_params p, fftwf_complex **outcpts, int step, char *label) {

  ptrdiff_t x_thickness, x_start, alloc_local;

  ptrdiff_t n0 = p.Lx;
  ptrdiff_t n1 = p.Ly;
  ptrdiff_t n2 = p.Lz;

  //  int slab;

  MPI_Status status;

  float kmode, modsq;

  int x, y, z;
  int i;

  // *p.a*p.a*p.a*p.dx*p.dx*p.dx

  if(label != NULL){
      if(*label){
	printf0(p, "Starting %s TT tensor PS calculation.\n", label);
      }
  }
  else{
    printf0(p, "Starting TT tensor PS calculation.\n", label);
  }
  float start = clock();


  alloc_local = fftwf_mpi_local_size_3d(n0, n1, n2,
				       MPI_COMM_WORLD,
				       &x_thickness, &x_start);

  
  float *slice = (float *)malloc(x_thickness*p.Ly*p.Lz*sizeof(float));

  printf0(p, "Transverse traceless projection...\n");
  // At this point, we should have a normed FFT

  // Turn the FFT'd tensor into transverse traceless gauge.
  tens_proj(p, x_start, x_thickness, slice, outcpts);

  if(label != NULL){
    if(*label){
      printf0(p, "Calculating total mod square of %s TT tensor.\n", label);
    }
  }
  else{
    printf0(p, "Calculating total mod square of TT tensor.\n");
  }
  float local_TT = 0.0;

  for(x=0; x<x_thickness; x++) {
    for(y=0; y<p.Ly; y++) {
      for(z=0; z<p.Lz; z++) {
	local_TT += slice[x*p.Ly*p.Lz + y*p.Lz + z];
	
      }
    }
  }



  float total_TT = reduce_sum(local_TT, p);

  if(label != NULL){
    if(*label){
      printf0(p,
	      "Unnormalised mod square of %s TT tensor [density]"
	      "claimed %6.10lf\n", label, total_TT);
    }
  }
  else{
    printf0(p,
	    "Unnormalised mod square of TT tensor [density]"
	    "claimed %6.10lf\n", total_TT);
  }
  

  // Bin the power spectrum on the fly

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

	// For binning we use momentum space index
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
  //  float comovingk, thisf, thisdiff, thisomega;

  // not sure (includes h^2)
  //  float omegarad = 3.5e-5;

  //  float doffrac = pow(1.0/100.0,1.0/3.0);

  // spokesman does the final analysis
  if(!p.rank) {

  char fftdest[200];

  if(label != NULL){
      if(*label){
	sprintf(fftdest,"%s-TT-ps-step%d.txt",label, step);
      }
  }
  else{
    sprintf(fftdest,"TT-ps-step%d.txt", step);
  }
  
  FILE *fp = fopen(fftdest,"w");

    for(i=0;i<nbins;i++) {
      
      fprintf(fp, "%lf %g %d\n",
	      thisk/(p.dx*p.a), ((float)(i+1))*bins[i], counts[i]);

      thisk = thisk + dk;
    }

  fclose(fp);

  }

  free(bins);
  free(counts);


  // Tidy up

  free(slice);


  float end = clock();

  if(label != NULL){
    if(*label){
      printf0(p, "%s TT tensor PS calculation took %lf\n", label,
	      ((float) (end - start)) / CLOCKS_PER_SEC);
    }
  }
  else{
  printf0(p, "TT tensor PS calculation took %lf\n",
	  ((float) (end - start)) / CLOCKS_PER_SEC);
  }

  return total_TT;
}




#endif // FFT
