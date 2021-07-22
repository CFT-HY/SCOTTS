/** @file uetc.c
 *
 * File containing unequal time correlator calculations.
 *
 * FFT functions in fft_func.c,
 * Projectors are in projectors.c,
 * Tensor power spectra are in tensorps.c,
 * Vector power spectra are in vectorps.c,
 * Scalar power spectrum calculation is in scalarps.c.
 */
#include "hydro.h"


#ifdef FFT




/**  Initialise the UETC calculation for all fields for which we
 * calculate UETCs with the momentum space
 */
void init_uetc(hydro_fields f, hydro_params p, fft_fields fft_f) {


  int i, x, y, z;

  if(p.fft_shst){
  // Initialise UETC for energy momentum tensor.
  float ****Tij = make_tensor(p);
  stress_energy(f, p, Tij);

  fft_tensor(p, fft_f, Tij, fft_f.initial_Tij, 1.0);
  
  
  free_tensor(p, Tij);
  }
#ifndef SCALAR
  if(p.fft_vel){
    fft_vector(p, fft_f, f.V, fft_f.initial_V);
  }
#endif
}


/**  Actually calculate the transverse traceless uetc of a tensor between two timesteps and store in
 * an appropriately named file.
 */
void uetc_tensor(hydro_params p, fftwf_complex **tensor_then,
		 fftwf_complex **tensor_now, int step_then, int step_now,
		 char *label){

  MPI_Status status;

  long ksitesq;
  
  int x, y, z;
  int i;

  if(label != NULL){
    if(*label){
      printf0(p, "Starting %s tensor (TT) UETC calculation.\n",label);
    }
  }
  else{
   printf0(p, "Starting tensor (TT) UETC calculation.\n");
  }
  float start = clock();

  ptrdiff_t alloc_local, x_thickness, x_start;

  ptrdiff_t n0 = p.Lx;
  ptrdiff_t n1 = p.Ly;
  ptrdiff_t n2 = p.Lz;
  MPI_Comm fftw_comm;
  int stride = (p.size > p.Lx) ? ((int)(p.size/p.Lx)) : 1;
  int color = p.rank%stride;
  MPI_Comm_split(MPI_COMM_WORLD, color, p.rank,
		 &fftw_comm);

  if(color == 0) {
    alloc_local = fftwf_mpi_local_size_3d(n0, n1, n2,
					  fftw_comm,
					  &x_thickness,
					  &x_start);
  } else{
    alloc_local = 0;
    x_thickness = 0;
    x_start = 0;
  }


  /////////////// Project to get transverse traceless tensor uetc ///////////////////

  float *slice_re = (float *)malloc(x_thickness*p.Ly*p.Lz*sizeof(float));
  float *slice_im = (float *)malloc(x_thickness*p.Ly*p.Lz*sizeof(float));


  // could make strictly > because equal time is just shear stress
  if(step_now >= p.uetcstart) {

    printf0(p,"projecting\n");

    uetc_tens_proj(p, x_start, x_thickness, slice_re, slice_im,
		 tensor_then, tensor_now);


    printf0(p,"binning\n");

    int nbins = minof3_int(p.Lx, p.Ly, p.Lz);
    float mink = 0.0;
    float maxk = 2.0*M_PI;

    float dk = (maxk-mink)/((float)nbins);

    float *bins_re = (float *)malloc(nbins*sizeof(float));
    float *bins_im = (float *)malloc(nbins*sizeof(float));
    int *counts = (int *)malloc(nbins*sizeof(int));

    for(i=0;i<nbins;i++) {
      bins_re[i] = 0.0;
      bins_im[i] = 0.0;
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
	  ksitesq = true_x*true_x + true_y*true_y + true_z*true_z;
	
	  whichbin = (int)sqrt(ksitesq);
	  bins_re[whichbin] += slice_re[x*p.Ly*p.Lz + y*p.Lz + z];
	  bins_im[whichbin] += slice_im[x*p.Ly*p.Lz + y*p.Lz + z];
	  counts[whichbin]++;

	}
      }
    }


    float red_value;
    int red_count;

    for(i=0;i<nbins;i++) {

      red_value = reduce_sum(bins_re[i], p);
      bins_re[i] = red_value;

      red_value = reduce_sum(bins_im[i], p);
      bins_im[i] = red_value;

      red_count = reduce_sum_int(counts[i], p);
      counts[i] = red_count;
    }


    float thisk = dk/2.0;

    // spokesman does the final analysis
    if(!p.rank) {

      char fftdest[200];

      if(label != NULL){
	if(*label){
	  sprintf(fftdest,"%s-TT-uetc-ref%d-step%d.txt", label, step_then, step_now);
	}
      }
      else{
	sprintf(fftdest,"TT-uetc-ref%d-step%d.txt", step_then, step_now);
      }
      FILE *fp = fopen(fftdest,"w");

      for(i=0;i<nbins;i++) {

	fprintf(fp, "%lf %g %g %d %d %d\n",
		thisk, bins_re[i], bins_im[i], counts[i], step_then, step_now);

	thisk = thisk + dk;
      }

      fclose(fp);

    }

    free(bins_re);
    free(bins_im);
    free(counts);
  }

  // Tidy up

  free(slice_re);
  free(slice_im);


  MPI_Comm_free(&fftw_comm);
  float end = clock();

  if(label != NULL){
    if(*label){
      printf0(p, "%s tensor (TT) UETC calculation took %lf\n",label,
	      ((float) (end - start)) / CLOCKS_PER_SEC);
    }
  }
  else{
    printf0(p, "tensor (TT) UETC calculation took %lf\n",
	    ((float) (end - start)) / CLOCKS_PER_SEC);
  }
}


/**  Calculate the rotatioanl and divergent components of the uetc of
 * a vector field between two timesteps and store them in an
 * appropriately named file.
 */
void uetc_vector(hydro_params p, fftwf_complex **vector_then,
		 fftwf_complex **vector_now, int step_then, int step_now,
		 char *label){

  MPI_Status status;

  long ksitesq;

  int x, y, z;
  int i;

  if(label != NULL){
    if(*label){
      printf0(p, "Starting %s vector UETC calculation.\n",label);
    }
  }
  else{
   printf0(p, "Starting vector UETC calculation.\n");
  }
  float start = clock();

  ptrdiff_t alloc_local, x_thickness, x_start;

  ptrdiff_t n0 = p.Lx;
  ptrdiff_t n1 = p.Ly;
  ptrdiff_t n2 = p.Lz;
  MPI_Comm fftw_comm;
  int stride = (p.size > p.Lx) ? ((int)(p.size/p.Lx)) : 1;
  int color = p.rank%stride;
  MPI_Comm_split(MPI_COMM_WORLD, color, p.rank,
		 &fftw_comm);

  if(color == 0) {
    alloc_local = fftwf_mpi_local_size_3d(n0, n1, n2,
					  fftw_comm,
					  &x_thickness,
					  &x_start);
  } else{
    alloc_local = 0;
    x_thickness = 0;
    x_start = 0;
  }


  /////////////// Project to get rot and div components of vector uetc ///////////////////

  float *slice_rot_re = (float *)malloc(x_thickness*p.Ly*p.Lz*sizeof(float));
  float *slice_rot_im = (float *)malloc(x_thickness*p.Ly*p.Lz*sizeof(float));
  float *slice_div_re = (float *)malloc(x_thickness*p.Ly*p.Lz*sizeof(float));
  float *slice_div_im = (float *)malloc(x_thickness*p.Ly*p.Lz*sizeof(float));


  // could make strictly > because equal time is just shear stress
  if(step_now >= p.uetcstart) {

    printf0(p,"projecting\n");

    uetc_split_vector(p, x_start, x_thickness, slice_rot_re, slice_rot_im,
		      slice_div_re, slice_div_im, vector_then, vector_now);


    printf0(p,"binning\n");

    int nbins = minof3_int(p.Lx, p.Ly, p.Lz);
    float mink = 0.0;
    float maxk = 2.0*M_PI;

    float dk = (maxk-mink)/((float)nbins);
    
    float *bins_rot_re = (float *)malloc(nbins*sizeof(float));
    float *bins_rot_im = (float *)malloc(nbins*sizeof(float));
    float *bins_div_re = (float *)malloc(nbins*sizeof(float));
    float *bins_div_im = (float *)malloc(nbins*sizeof(float));

    int *counts = (int *)malloc(nbins*sizeof(int));

    for(i=0;i<nbins;i++) {
      bins_rot_re[i] = 0.0;
      bins_rot_im[i] = 0.0;
      bins_div_re[i] = 0.0;
      bins_div_im[i] = 0.0;

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
	  ksitesq = true_x*true_x + true_y*true_y + true_z*true_z;
	
	  whichbin = (int)sqrt(ksitesq);

	  bins_rot_re[whichbin] += slice_rot_re[x*p.Ly*p.Lz + y*p.Lz + z];
	  bins_rot_im[whichbin] += slice_rot_im[x*p.Ly*p.Lz + y*p.Lz + z];
	  bins_div_re[whichbin] += slice_div_re[x*p.Ly*p.Lz + y*p.Lz + z];
	  bins_div_im[whichbin] += slice_div_im[x*p.Ly*p.Lz + y*p.Lz + z];

	  counts[whichbin]++;

	}
      }
    }

    MPI_Barrier(MPI_COMM_WORLD);

    float red_value;
    int red_count;

    for(i=0;i<nbins;i++) {

      red_value = reduce_sum(bins_rot_re[i], p);
      bins_rot_re[i] = red_value;

      red_value = reduce_sum(bins_rot_im[i], p);
      bins_rot_im[i] = red_value;

      red_value = reduce_sum(bins_div_re[i], p);
      bins_div_re[i] = red_value;

      red_value = reduce_sum(bins_div_im[i], p);
      bins_div_im[i] = red_value;


      red_count = reduce_sum_int(counts[i], p);
      counts[i] = red_count;
    }


    float thisk = dk/2.0;

    // spokesman does the final analysis

    // Do rot first
    if(!p.rank) {

      char fftdest[200];

      if(label != NULL){
	if(*label){
	  sprintf(fftdest,"%s-rot-uetc-ref%d-step%d.txt", label, step_then, step_now);
	}
      }
      else{
	sprintf(fftdest,"rot-uetc-ref%d-step%d.txt", step_then, step_now);
      }
      FILE *fp1 = fopen(fftdest,"w");

      for(i=0;i<nbins;i++) {

	fprintf(fp1, "%lf %g %g %d %d %d\n",
		thisk/(p.a*p.dx), (thisk/dk)*bins_rot_re[i], (thisk/dk)*bins_rot_im[i], counts[i], step_then, step_now);

	thisk = thisk + dk;
      }

      // Now div

      thisk = dk/2.0;

      if(label != NULL){
	if(*label){
	  sprintf(fftdest,"%s-div-uetc-ref%d-step%d.txt", label, step_then, step_now);
	}
      }
      else{
	sprintf(fftdest,"div-uetc-ref%d-step%d.txt", step_then, step_now);
      }
      FILE *fp2 = fopen(fftdest,"w");

      for(i=0;i<nbins;i++) {

	fprintf(fp2, "%lf %g %g %d %d %d\n",
		thisk/(p.a*p.dx), (thisk/dk)*bins_div_re[i], (thisk/dk)*bins_div_im[i], counts[i], step_then, step_now);

	thisk = thisk + dk;
      }


      fclose(fp1);
      fclose(fp2);

    }

    free(bins_rot_re);
    free(bins_rot_im);
    free(bins_div_re);
    free(bins_div_im);

    free(counts);
  }

  // Tidy up

  free(slice_rot_re);
  free(slice_rot_im);
  free(slice_div_re);
  free(slice_div_im);


  MPI_Comm_free(&fftw_comm);

  float end = clock();

  if(label != NULL){
    if(*label){
      printf0(p, "%s vector UETC calculation took %lf\n",label,
	      ((float) (end - start)) / CLOCKS_PER_SEC);
    }
  }
  else{
    printf0(p, "vector UETC calculation took %lf\n",
	    ((float) (end - start)) / CLOCKS_PER_SEC);
  }
}


#endif // FFT
