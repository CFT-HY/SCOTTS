/* velps.c
 *
 * Fourier transform and velocity/vorticity power spectrum.
 */
#include "hydro.h"


#ifdef FFT


/* vel_proj(int T, double kx, double ky, double kz);
 */
double vel_proj(int T, double kx, double ky, double kz) {

  double mag = sqrt(kx*kx + ky*ky + kz*kz);

  // Avoid overflow
  if(fabs(mag) < 0.000001)
    return 0.0;

  kx = kx/mag;
  ky = ky/mag;
  kz = kz/mag;

  double total = 0.0;

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
    fprintf(stderr,"Unknown velocity projector element! Nonsense will ensue...\n");

  }

  return 0.0;
}


void split_and_power(hydro_params p, int x_start, int slab,
		     double *product, double *product_div, fftw_complex **vk) {

  int i, j;
  int x, y, z;
  double kx, ky, kz;

  double res_r, res_i;
  double resid_r, resid_i;

  for(x=0; x<slab; x++) {
    for(y=0; y<p.Ly; y++) {
      for(z=0; z<p.Lz; z++) {

        kx = ((double)(x+x_start))*2.0*M_PI/((double)p.Lx);
        ky = ((double)y)*2.0*M_PI/((double)p.Ly);
        kz = ((double)z)*2.0*M_PI/((double)p.Lz);

        product[x*p.Ly*p.Lz + y*p.Lz + z] = 0.0;
        product_div[x*p.Ly*p.Lz + y*p.Lz + z] = 0.0;


        for(i=1; i<=3; i++) {

	  res_r = 0.0;
	  res_i = 0.0;

	  resid_r = 0.0;
	  resid_i = 0.0;
	  
          for(j=1; j<=3; j++) {
	    res_r += vel_proj(i*10 + j, kx, ky, kz)*vk[j-1][x*p.Ly*p.Lz + y*p.Lz + z][0];
	    res_i += vel_proj(i*10 + j, kx, ky, kz)*vk[j-1][x*p.Ly*p.Lz + y*p.Lz + z][1];

	    resid_r += (1.0 - vel_proj(i*10 + j, kx, ky, kz))*vk[j-1][x*p.Ly*p.Lz + y*p.Lz + z][0];
	    resid_i += (1.0 - vel_proj(i*10 + j, kx, ky, kz))*vk[j-1][x*p.Ly*p.Lz + y*p.Lz + z][1];
	  }

	  product[x*p.Ly*p.Lz + y*p.Lz + z] +=
	    res_r*res_r + res_i*res_i;

	  product_div[x*p.Ly*p.Lz + y*p.Lz + z] +=
	    resid_r*resid_r + resid_i*resid_i;

	}

      }
    }
  }

}


void histogram(hydro_params p, double *slice, char *filename, int slab, int x_start) {
  // The rest proceeds as in gw.c

  int i, x, y, z;

  double kmode;

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



  for(x=0;x<slab;x++) {
    for(y=0;y<p.Ly;y++) {
      for(z=0;z<p.Lz;z++) {

        if(((x+x_start)>p.Lx/2) || (y> p.Ly/2) || (z>p.Lz/2))
          continue;
	
	kmode = sqrt(
		     ((double)((x+x_start)*(x+x_start)))/((double)(p.Lx*p.Lx))
		     + ((double)(y*y))/((double)(p.Ly*p.Ly))
		     + ((double)(z*z))/((double)(p.Lz*p.Lz))
		     )*2.0*M_PI;

	whichbin = (int)floor(kmode/dk);
	bins[whichbin] += slice[x*p.Ly*p.Lz + y*p.Lz + z];
	counts[whichbin]++;

	//	if(step)
	// fprintf(stderr,"%d %lf\n",whichbin,slice[x*p.Ly*p.Lz + y*p.Lz + z]);

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

  // spokesman does the final analysis
  if(!p.rank) {


  FILE *fp = fopen(filename,"w");

    for(i=0;i<nbins;i++) {


      fprintf(fp, "%lf %g %d\n",
	      thisk, bins[i], counts[i]);

      thisk = thisk + dk;
    }

  fclose(fp);

  }

  free(bins);
  free(counts);
}





/* fft_vel(hydro_fields f, hydro_params p)
 *
 * Calculate the velocity power spectrum.
 */
void fft_vel(hydro_fields f, hydro_params p, int step) {

  ptrdiff_t x_thickness, x_start, alloc_local;

  ptrdiff_t n0 = p.Lx;
  ptrdiff_t n1 = p.Ly;
  ptrdiff_t n2 = p.Lz;

  int slab;

  MPI_Status status;

  double modsq;

  int x, y, z;
  int i;

  double fft_norm = p.a*p.a*p.a*p.dx*p.dx*p.dx/pow(2.0*M_PI,1.5);

  double *trim = (double *)malloc(p.slicex*p.slicey*p.Lz*sizeof(double));

  printf0(p, "Starting FFT stuff.\n");

  double start = clock();

  fftw_mpi_init();

  alloc_local = fftw_mpi_local_size_3d(n0, n1, n2,
				       MPI_COMM_WORLD, &x_thickness, &x_start);


  slab = (int) x_thickness;

  if(((int)x_thickness) != p.Lx/p.size) {
    fprintf(stderr,
	    "Giving up in FFT: FFTW told me to use a silly layout!\n");
    die(-42);
  }

  if(((int)x_thickness) == 0) {
    fprintf(stderr, "Giving up in FFT: dx=0!\n");
    die(-43);
  }


  fftw_complex *in = fftw_alloc_complex(alloc_local);
  fftw_complex *out = fftw_alloc_complex(alloc_local);

  fftw_complex **outcpts = (fftw_complex **)malloc(3*sizeof(fftw_complex *));
  
  for(i=0;i<3;i++) {
    outcpts[i] = fftw_alloc_complex(alloc_local);
  }

  double *slice = (double *)malloc(slab*p.Ly*p.Lz*sizeof(double));
  double *slice_div = (double *)malloc(slab*p.Ly*p.Lz*sizeof(double));


  int nx = p.Lx/p.slicex;
  int ny = p.Ly/p.slicey;

  int ry;

  // Now planning  
  fftw_plan plan = fftw_mpi_plan_dft_3d(p.Lx, p.Ly, p.Lz,
					in, out, MPI_COMM_WORLD,
					FFTW_BACKWARD, FFTW_ESTIMATE);

  for(i = 0; i < 3; i++) {

    double check = 0.0;    
    
    for(x=1; x<=p.slicex; x++) {
      for(y=1; y<=p.slicey; y++) {
	for(z=0; z<p.Lz; z++) {
	  check += f.V[i][x][y][z];
	  trim[(x-1)*p.slicey*p.Lz + (y-1)*p.Lz + z] = f.V[i][x][y][z];
	}
      }
    }

    //    fprintf(stderr,"Check: uij %lf\n", reduce_sum(check, p));
    
    for(x=0; x<p.Lx; x++) {
      
      // Check whether we are the target node for this slab
      if((x >= x_start) && ( x < (x_start + x_thickness))) {
	
	for(ry = 0; ry < ny; ry++) {
	  
	  if((x/p.slicex == p.myposx) && (ry == p.myposy)) {
	    
	    memcpy(&slice[(x-x_start)*p.Ly*p.Lz + ry*p.slicey*p.Lz],
		   &trim[(x-p.shiftx)*p.slicey*p.Lz],
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
	}
	
	
      } else if((x >= p.shiftx) && (x < (p.shiftx + p.slicex))) {
	
	
      MPI_Send(&trim[(x-p.shiftx)*p.slicey*p.Lz],
	       p.slicey*p.Lz,
	       MPI_DOUBLE,
	       x/slab,
	       x*ny + p.myposy,
	       MPI_COMM_WORLD);
      }
      
    }

    double total = 0.0;

    for(x=0;x<slab;x++) {
      for(y=0;y<p.Ly;y++) {
	for(z=0;z<p.Lz;z++) {
	  in[x*p.Ly*p.Lz + y*p.Lz + z][0] = slice[x*p.Ly*p.Lz + y*p.Lz + z];
	  in[x*p.Ly*p.Lz + y*p.Lz + z][1] = 0.0;
	  total += in[x*p.Ly*p.Lz + y*p.Lz + z][0];
	}
      }
    }


    /*
    if(!i) {
      fprintf(stderr,"total incpts %d is %6.10lf\n", p.rank,
	      reduce_sum(total,p));
    }

    if(!p.rank && !i)
      fprintf(stderr,"inslice 114: %.12lf\n", in[p.Ly*p.Lz + p.Lz + 4][0]);
    */

  
    fftw_mpi_execute_dft(plan, in, out);
    
    total = 0.0;
    
    for(x=0;x<slab;x++) {
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

  // At this point, we should have a normed FFT of each component in outcpts[i][k]



  // The brains of the operation:
  // Turn the FFT'd tensor into projected power spectrum
  split_and_power(p, x_start, slab, slice, slice_div, outcpts);



  char fftdest[200];

  sprintf(fftdest,"rot-ps-step%d.txt", step);
  histogram(p, slice, fftdest, slab, x_start);
  sprintf(fftdest,"div-ps-step%d.txt", step);
  histogram(p, slice_div, fftdest, slab, x_start);


  // Tidy up

  free(slice);
  free(slice_div);
  free(trim);

  fftw_destroy_plan(plan);
  
  fftw_free(in);
  fftw_free(out);

  for(i=0;i<3;i++)
    fftw_free(outcpts[i]);


  fftw_mpi_cleanup();

  double end = clock();

  printf0(p, "FFT stuff took %lf\n",
	  ((double) (end - start)) / CLOCKS_PER_SEC);

}




#endif // FFT
