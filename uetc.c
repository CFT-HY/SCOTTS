/* uetc.c
 */
#include "hydro.h"


#ifdef FFT




/* void gwproject(hydro_params p, int x_start, int slab,
 *              double *product, fftw_complex **udot); 
 *
 * Project from the six unique degrees of freedom in udot,
 * (obtained here in momentum space from a distributed shared memory
 * FFT) to obtain |hij|^2 in momentum space, cf eq. 30 in the
 * paper cited at the start.
 *
 * slab is the thickness in the x-direction of the local domain
 * x_start is the physical start of this region
 * y and z have the full spatial extent p.Ly and p.Lz
 */ 
void init_uetc(hydro_fields f, hydro_params p) {


  int i, x, y, z;

  double ****Tij = make_tensor(p);
  stress_energy(f, p, Tij);


  for(i = 0; i < TENSOR_CPTS; i++) {
    for(x = 1; x <= p.slicex; x++) {
      for(y = 1; y <= p.slicey; y++) {
        for(z = 0; z < p.Lz; z++) {

          f.initial_Tij[i][x][y][z] = Tij[i][x][y][z];
        }
      }
    }
  }


  free_tensor(p, Tij);
}



void uetc_project(hydro_params p, int x_start, int slab,
		  double *product, fftw_complex **Tij_then, fftw_complex **Tij_now) {

  int i, j, l, m;
  int x, y, z;
  double kx, ky, kz;


  for(x=0; x<slab; x++) {
    for(y=0; y<p.Ly; y++) {
      for(z=0; z<p.Lz; z++) {

	kx = ((double)(x+x_start))*2.0*M_PI/((double)p.Lx);
	ky = ((double)y)*2.0*M_PI/((double)p.Ly);
	kz = ((double)z)*2.0*M_PI/((double)p.Lz);
	
	product[x*p.Ly*p.Lz + y*p.Lz + z] = 0.0;

	// Sum over indices on lambda_{ij,lm} and udot_{ij}(k)
	for(i=1; i<=3; i++) {
	  for(j=1; j<=3; j++) {
	    for(l=1; l<=3; l++) {
	      for(m=1; m<=3; m++) {

		product[x*p.Ly*p.Lz + y*p.Lz + z] +=
		  lambda(i, j, l, m, kx, ky, kz)
		  *(Tij_then[indexof(i,j)][x*p.Ly*p.Lz + y*p.Lz + z][0]
		   *Tij_now[indexof(l,m)][x*p.Ly*p.Lz + y*p.Lz + z][0]
		    + Tij_then[indexof(i,j)][x*p.Ly*p.Lz + y*p.Lz + z][1]
		    *Tij_now[indexof(l,m)][x*p.Ly*p.Lz + y*p.Lz + z][1]);
	      }
	    }
	  }
	}

      }
    }
  }
}




/* fft_uetc(hydro_fields f, hydro_params p)
 *
 */
void fft_uetc(hydro_fields f, hydro_params p, int step) {

  ptrdiff_t x_thickness, x_start, alloc_local;

  ptrdiff_t n0 = p.Lx;
  ptrdiff_t n1 = p.Ly;
  ptrdiff_t n2 = p.Lz;

  int slab;

  MPI_Status status;

  double kmode, modsq;

  int x, y, z;
  int i;

  double fft_norm = p.a*p.a*p.a*p.dx*p.dx*p.dx/pow(2.0*M_PI,1.5);

  double *trim_then = (double *)malloc(p.slicex*p.slicey*p.Lz*sizeof(double));
  double *trim_now = (double *)malloc(p.slicex*p.slicey*p.Lz*sizeof(double));

  printf0(p, "Starting UETC calculation.\n");

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

  /////////////// here compute other thing to be squished ////////////////

  double ****Tij = make_tensor(p);
  stress_energy(f, p, Tij);

  /////////////// done /////////////////


  fftw_complex *in = fftw_alloc_complex(alloc_local);
  fftw_complex *out = fftw_alloc_complex(alloc_local);

  fftw_complex **outcpts_then = (fftw_complex **)malloc(6*sizeof(fftw_complex *));
  fftw_complex **outcpts_now = (fftw_complex **)malloc(6*sizeof(fftw_complex *));
  
  for(i=0;i<TENSOR_CPTS;i++) {
    outcpts_then[i] = fftw_alloc_complex(alloc_local);
    outcpts_now[i] = fftw_alloc_complex(alloc_local);
  }

  double *slice_then = (double *)malloc(slab*p.Ly*p.Lz*sizeof(double));
  double *slice_now = (double *)malloc(slab*p.Ly*p.Lz*sizeof(double));

  double *slice = (double *)malloc(slab*p.Ly*p.Lz*sizeof(double));


  int nx = p.Lx/p.slicex;
  int ny = p.Ly/p.slicey;

  int ry;

  // Now planning  
  fftw_plan plan = fftw_mpi_plan_dft_3d(p.Lx, p.Ly, p.Lz,
					in, out, MPI_COMM_WORLD,
					FFTW_BACKWARD, FFTW_ESTIMATE);

  for(i = 0; i < TENSOR_CPTS; i++) {

    double check = 0.0;    
    
    for(x=1; x<=p.slicex; x++) {
      for(y=1; y<=p.slicey; y++) {
	for(z=0; z<p.Lz; z++) {
	  trim_then[(x-1)*p.slicey*p.Lz + (y-1)*p.Lz + z] = f.initial_Tij[i][x][y][z];
	  trim_now[(x-1)*p.slicey*p.Lz + (y-1)*p.Lz + z] = Tij[i][x][y][z];
	}
      }
    }

    //    fprintf(stderr,"Check: uij %lf\n", reduce_sum(check, p));
    
    for(x=0; x<p.Lx; x++) {
      
      // Check whether we are the target node for this slab
      if((x >= x_start) && ( x < (x_start + x_thickness))) {
	
	for(ry = 0; ry < ny; ry++) {
	  
	  if((x/p.slicex == p.myposx) && (ry == p.myposy)) {
	    
	    memcpy(&slice_then[(x-x_start)*p.Ly*p.Lz + ry*p.slicey*p.Lz],
		   &trim_then[(x-p.shiftx)*p.slicey*p.Lz],
		   p.slicey*p.Lz*sizeof(double));

	    memcpy(&slice_now[(x-x_start)*p.Ly*p.Lz + ry*p.slicey*p.Lz],
		   &trim_now[(x-p.shiftx)*p.slicey*p.Lz],
		   p.slicey*p.Lz*sizeof(double));
	    continue;
	  }
	  
	  
	  MPI_Recv(&slice_then[(x-x_start)*p.Ly*p.Lz + ry*p.slicey*p.Lz],
		   p.slicey*p.Lz,
		   MPI_DOUBLE,
		   ry*nx + x/p.slicex,
		   x*ny + ry,
		   MPI_COMM_WORLD,
		   &status);

	  MPI_Recv(&slice_now[(x-x_start)*p.Ly*p.Lz + ry*p.slicey*p.Lz],
		   p.slicey*p.Lz,
		   MPI_DOUBLE,
		   ry*nx + x/p.slicex,
		   p.Lx*ny + x*ny + ry,
		   MPI_COMM_WORLD,
		   &status);
	}
	
	
      } else if((x >= p.shiftx) && (x < (p.shiftx + p.slicex))) {
	
	
      MPI_Send(&trim_then[(x-p.shiftx)*p.slicey*p.Lz],
	       p.slicey*p.Lz,
	       MPI_DOUBLE,
	       x/slab,
	       x*ny + p.myposy,
	       MPI_COMM_WORLD);

      MPI_Send(&trim_now[(x-p.shiftx)*p.slicey*p.Lz],
	       p.slicey*p.Lz,
	       MPI_DOUBLE,
	       x/slab,
	       p.Lx*ny + x*ny + p.myposy,
	       MPI_COMM_WORLD);
      }
      
    }

    double total = 0.0;

    for(x=0;x<slab;x++) {
      for(y=0;y<p.Ly;y++) {
	for(z=0;z<p.Lz;z++) {
	  in[x*p.Ly*p.Lz + y*p.Lz + z][0] = slice_then[x*p.Ly*p.Lz + y*p.Lz + z];
	  in[x*p.Ly*p.Lz + y*p.Lz + z][1] = 0.0;
	}
      }
    }

    fftw_mpi_execute_dft(plan, in, out);

    for(x=0;x<slab;x++) {
      for(y=0;y<p.Ly;y++) {
	for(z=0;z<p.Lz;z++) {
	  outcpts_then[i][x*p.Ly*p.Lz + y*p.Lz + z][0] 
	    = fft_norm*out[x*p.Ly*p.Lz + y*p.Lz + z][0];
	  outcpts_then[i][x*p.Ly*p.Lz + y*p.Lz + z][1]
	    = fft_norm*out[x*p.Ly*p.Lz + y*p.Lz + z][1];

	}
      }
    }



    for(x=0;x<slab;x++) {
      for(y=0;y<p.Ly;y++) {
	for(z=0;z<p.Lz;z++) {
	  in[x*p.Ly*p.Lz + y*p.Lz + z][0] = slice_now[x*p.Ly*p.Lz + y*p.Lz + z];
	  in[x*p.Ly*p.Lz + y*p.Lz + z][1] = 0.0;
	}
      }
    }

    fftw_mpi_execute_dft(plan, in, out);

    for(x=0;x<slab;x++) {
      for(y=0;y<p.Ly;y++) {
	for(z=0;z<p.Lz;z++) {
	  outcpts_now[i][x*p.Ly*p.Lz + y*p.Lz + z][0] 
	    = fft_norm*out[x*p.Ly*p.Lz + y*p.Lz + z][0];
	  outcpts_now[i][x*p.Ly*p.Lz + y*p.Lz + z][1]
	    = fft_norm*out[x*p.Ly*p.Lz + y*p.Lz + z][1];

	}
      }
    }

  }


  uetc_project(p, x_start, slab, slice, outcpts_then, outcpts_now);


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

  char fftdest[200];

  sprintf(fftdest,"uetc-step%d.txt", step);

  FILE *fp = fopen(fftdest,"w");

    for(i=0;i<nbins;i++) {

      fprintf(fp, "%lf %g %d %d\n",
	      thisk, bins[i], counts[i], step);

      thisk = thisk + dk;
    }

  fclose(fp);

  }

  free(bins);
  free(counts);

  // Tidy up

  free(slice);
  free(slice_then);
  free(slice_now);
  free(trim_then);
  free(trim_now);

  fftw_destroy_plan(plan);
  
  fftw_free(in);
  fftw_free(out);

  for(i=0;i<TENSOR_CPTS;i++) {
    fftw_free(outcpts_then[i]);
    fftw_free(outcpts_now[i]);
  }


  free_tensor(p, Tij);

  fftw_mpi_cleanup();

  double end = clock();

  printf0(p, "UETC FFT calculation took %lf\n",
	  ((double) (end - start)) / CLOCKS_PER_SEC);

}




#endif // FFT
