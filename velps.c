/* velps.c
 *
 * Fourier transform and velocity/vorticity power spectrum.
 *
 * The gravitational wave power spectrum is in gw.c,
 * and a generic field power spectrum calculation is in fft.c.
 */
#include "hydro.h"


#ifdef FFT


/* vel_proj(int T, float kx, float ky, float kz);
 *
 * Calculate the projector for the transverse component of the
 * velocity.
 */
float vel_proj(int T, float kx, float ky, float kz) {

  float mag = sqrt(kx*kx + ky*ky + kz*kz);

  // Avoid overflow
  if(fabs(mag) < 0.000001)
    return 0.0;

  kx = kx/mag;
  ky = ky/mag;
  kz = kz/mag;

  float total = 0.0;

  switch(T) {

  case 11:
    return 1.0 - 1.0*kx*kx;
    break;
  case 21:
    return -1.0*kx*ky;
    break;
  case 31:
    return -1.0*kx*kz;
    break;
  case 12:
    return -1.0*kx*ky;
    break;
  case 13:
    return -1.0*kx*kz;
    break;
  case 22:
    return 1.0 - 1.0*ky*ky;
    break;
  case 32:
    return -1.0*ky*kz;
    break;
  case 23:
    return -1.0*ky*kz;
    break;
  case 33:
    return 1.0 - 1.0*kz*kz;
    break;

  default:
    fprintf(stderr,
	    "Unknown velocity projector element! Nonsense will ensue...\n");

  }

  return 0.0;
}


/* void split_and_power(hydro_params p, int x_start, int slab,
 *                    float *product, float *product_div,
 *                    fftwf_complex **vk)
 *
 * Split the FFT output vk into the transverse and longitudinal ('div')
 * components.
 */
void split_and_power(hydro_params p, int x_start, int slab,
		     float *product, float *product_div, float *product_tot,
		     fftwf_complex **vk) {

  int i, j;
  int x, y, z;
  float kx, ky, kz;

  float res_r, res_i;
  float resid_r, resid_i;
  float tot_r, tot_i;

  int true_x, true_y, true_z;

  for(x=0; x<slab; x++) {
    for(y=0; y<p.Ly; y++) {
      for(z=0; z<p.Lz; z++) {

        if(x+x_start > p.Lx/2)
          true_x = -(p.Lx - (x+x_start));
        else
          true_x = x+x_start;

        if(y > p.Ly/2)
          true_y = -(p.Ly - y);
        else
          true_y = y;

        if(z > p.Lz/2)
          true_z = -(p.Lz - z);
	else
          true_z = z;


        kx = sqrt((2.0 - 2.0*cos(((float)(true_x))*2.0*M_PI/(((float)p.Lx))))/(p.dx*p.dx));
	ky = sqrt((2.0 - 2.0*cos(((float)(true_y))*2.0*M_PI/(((float)p.Ly))))/(p.dx*p.dx));
        kz = sqrt((2.0 - 2.0*cos(((float)(true_z))*2.0*M_PI/(((float)p.Lz))))/(p.dx*p.dx));


	//        kx = ((float)(x+x_start))*2.0*M_PI/((float)p.Lx);
	//        ky = ((float)y)*2.0*M_PI/((float)p.Ly);
	//        kz = ((float)z)*2.0*M_PI/((float)p.Lz);

        product[x*p.Ly*p.Lz + y*p.Lz + z] = 0.0;
        product_div[x*p.Ly*p.Lz + y*p.Lz + z] = 0.0;
        product_tot[x*p.Ly*p.Lz + y*p.Lz + z] = 0.0;


        for(i=1; i<=3; i++) {

	  res_r = 0.0;
	  res_i = 0.0;

	  resid_r = 0.0;
	  resid_i = 0.0;

	  tot_r = vk[i-1][x*p.Ly*p.Lz + y*p.Lz + z][0];
	  tot_i = vk[i-1][x*p.Ly*p.Lz + y*p.Lz + z][1];

          for(j=1; j<=3; j++) {

	    // Transverse components
	    res_r += vel_proj(i*10 + j, kx, ky, kz)
	      *vk[j-1][x*p.Ly*p.Lz + y*p.Lz + z][0];
	    res_i += vel_proj(i*10 + j, kx, ky, kz)
	      *vk[j-1][x*p.Ly*p.Lz + y*p.Lz + z][1];

	    // And longitudinal...
	    if(i == j) {
	    resid_r += (1.0 - vel_proj(i*10 + j, kx, ky, kz))
	      *vk[j-1][x*p.Ly*p.Lz + y*p.Lz + z][0];
	    resid_i += (1.0 - vel_proj(i*10 + j, kx, ky, kz))
	      *vk[j-1][x*p.Ly*p.Lz + y*p.Lz + z][1];
	    } else {
	    resid_r += (-1.0*vel_proj(i*10 + j, kx, ky, kz))
	      *vk[j-1][x*p.Ly*p.Lz + y*p.Lz + z][0];
	    resid_i += (-1.0*vel_proj(i*10 + j, kx, ky, kz))
	      *vk[j-1][x*p.Ly*p.Lz + y*p.Lz + z][1];
	    }
	  }

	  product[x*p.Ly*p.Lz + y*p.Lz + z] +=
	    res_r*res_r + res_i*res_i;

	  product_div[x*p.Ly*p.Lz + y*p.Lz + z] +=
	    resid_r*resid_r + resid_i*resid_i;

	  product_tot[x*p.Ly*p.Lz + y*p.Lz + z] +=
	    tot_r*tot_r + tot_i*tot_i;

	}

      }
    }
  }

}



/* void histogram(hydro_params p, float *slice, char *filename,
 *               int slab, int x_start)
 *
 * Combine the data from each slice (supplied separately by each
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





/* fft_vel(hydro_fields f, hydro_params p)
 *
 * Main method. Calculate the velocity power spectrum.
 */
void fft_vel(hydro_fields f, hydro_params p, int step, float ****vectorfield) {

#ifndef SCALAR
  ptrdiff_t x_thickness, x_start, alloc_local;

  ptrdiff_t n0 = p.Lx;
  ptrdiff_t n1 = p.Ly;
  ptrdiff_t n2 = p.Lz;

  int slab;

  MPI_Status status;

  float modsq;

  int x, y, z;
  int i;

  float fft_norm = (1.0/(((float)p.Lx)*((float)p.Ly)*((float)p.Lz)));
  // *p.a*p.a*p.a*p.dx*p.dx*p.dx

  float *trim = (float *)malloc(p.slicex*p.slicey*p.Lz*sizeof(float));

  printf0(p, "Starting velocity power spectrum FFT calculation.\n");

  float start = clock();

  fftwf_mpi_init();

  alloc_local = fftwf_mpi_local_size_3d(n0, n1, n2,
				       MPI_COMM_WORLD, &x_thickness, &x_start);


  // slab is the thickness the 'used' nodes wil FFT
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



  fftwf_complex *in = fftwf_alloc_complex(alloc_local);
  fftwf_complex *out = fftwf_alloc_complex(alloc_local);

  fftwf_complex **outcpts = (fftwf_complex **)malloc(3*sizeof(fftwf_complex *));
  
  for(i=0;i<3;i++) {
    outcpts[i] = fftwf_alloc_complex(alloc_local);
  }

  float *slice = (float *)malloc(slab*p.Ly*p.Lz*sizeof(float));
  float *slice_div = (float *)malloc(slab*p.Ly*p.Lz*sizeof(float));
  float *slice_tot = (float *)malloc(slab*p.Ly*p.Lz*sizeof(float));


  int nx = p.Lx/p.slicex;
  int ny = p.Ly/p.slicey;

  int ry;

  // Now planning  
  fftwf_plan plan = fftwf_mpi_plan_dft_3d(p.Lx, p.Ly, p.Lz,
					in, out, MPI_COMM_WORLD,
					FFTW_BACKWARD, FFTW_ESTIMATE);

  for(i = 0; i < 3; i++) {

    float check = 0.0;    
    
    for(x=1; x<=p.slicex; x++) {
      for(y=1; y<=p.slicey; y++) {
	for(z=0; z<p.Lz; z++) {
	  check += vectorfield[i][x][y][z];
	  trim[(x-1)*p.slicey*p.Lz + (y-1)*p.Lz + z] = vectorfield[i][x][y][z];
	}
      }
    }

    //    fprintf(stderr,"Check: uij %lf\n", reduce_sum(check, p));
    
    for(x=0; x<p.Lx; x++) {
      
      // Check whether we are the target node for this slab
      if((x >= x_start) && ( x < (x_start + x_thickness))) {
	
	for(ry = 0; ry < ny; ry++) {

	  //	  fprintf(stderr, "%d: my slice, x=%d\n", p.rank, x);
	  
	  if((x/p.slicex == p.myposx) && (ry == p.myposy)) {
	    
	    memcpy(&slice[(x-x_start)*p.Ly*p.Lz + ry*p.slicey*p.Lz],
		   &trim[(x-p.shiftx)*p.slicey*p.Lz],
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
	}
	
	//	fprintf(stderr, "%d: got pencil with x=%d\n", p.rank, x);

	
      } else if((x >= p.shiftx) && (x < (p.shiftx + p.slicex))) {
	
	
	/*
	fprintf(stderr, "%d: my pencil, x=%d to dest %d/%d\n",
		p.rank, x, x, slab);
	*/
	
      MPI_Send(&trim[(x-p.shiftx)*p.slicey*p.Lz],
	       p.slicey*p.Lz,
	       MPI_FLOAT,
	       x/slab,
	       x*ny + p.myposy,
	       MPI_COMM_WORLD);
      }
      
    }

    float total = 0.0;

    printf0(p, "Vel FFT: Done slicing for cpt %d\n", i);

    for(x=0;x<x_thickness;x++) {
      for(y=0;y<p.Ly;y++) {
	for(z=0;z<p.Lz;z++) {
	  in[x*p.Ly*p.Lz + y*p.Lz + z][0] = slice[x*p.Ly*p.Lz + y*p.Lz + z];
	  in[x*p.Ly*p.Lz + y*p.Lz + z][1] = 0.0;
	  total += in[x*p.Ly*p.Lz + y*p.Lz + z][0];
	}
      }
    }
  
    fftwf_mpi_execute_dft(plan, in, out);
    
    total = 0.0;
    
    for(x=0;x<x_thickness;x++) {
      for(y=0;y<p.Ly;y++) {
	for(z=0;z<p.Lz;z++) {
	  outcpts[i][x*p.Ly*p.Lz + y*p.Lz + z][0] 
	    = fft_norm*out[x*p.Ly*p.Lz + y*p.Lz + z][0];
	  outcpts[i][x*p.Ly*p.Lz + y*p.Lz + z][1]
	    = fft_norm*out[x*p.Ly*p.Lz + y*p.Lz + z][1];

	  //	  total += outcpts[i][x*p.Ly*p.Lz + y*p.Lz + z][0];
	}
      }
    }


  }

  /*
   * At this point, we should have a normed FFT
   * of each component in outcpts[i][k]
   */


  // The brains of the operation:
  // Turn the FFT'd tensor into projected power spectrum
  split_and_power(p, x_start, x_thickness, slice, slice_div, slice_tot, outcpts);



  char fftdest[200];

  sprintf(fftdest,"rot-ps-step%d.txt", step);
  histogram(p, slice, fftdest, x_thickness, x_start);
  sprintf(fftdest,"div-ps-step%d.txt", step);
  histogram(p, slice_div, fftdest, x_thickness, x_start);
  sprintf(fftdest,"tot-ps-step%d.txt", step);
  histogram(p, slice_tot, fftdest, x_thickness, x_start);


  // Tidy up

  free(slice);
  free(slice_div);
  free(slice_tot);
  free(trim);

  fftwf_destroy_plan(plan);
  
  fftwf_free(in);
  fftwf_free(out);

  for(i=0;i<3;i++)
    fftwf_free(outcpts[i]);
  free(outcpts);

  fftwf_mpi_cleanup();

  float end = clock();

  printf0(p, "Velocity power spectrum FFT calculation took %lf\n",
	  ((float) (end - start)) / CLOCKS_PER_SEC);

#endif // SCALAR
}




#endif // FFT
