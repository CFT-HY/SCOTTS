/* uetc.c
 *
 * Unequal time correlator calculation
 */
#include "hydro.h"


#ifdef FFT




/* void init_uetc(hydro_fields f, hydro_params p)
 *
 * Initialise the UETC calculation with the Tij's from the
 * current timestep.
 */ 
void init_uetc(hydro_fields f, hydro_params p) {


  int i, x, y, z;

  float ****Tij = make_tensor(p);
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


/* void uetc_project(hydro_params p, int x_start, int slab,
 *                 float *product_re, float *product_im,
 *                 fftwf_complex **Tij_then, fftwf_complex **Tij_now)
 *
 * Like gwproject (see gw.c), but for two separate timesteps, and takes two
 * different tensors.
 */
void uetc_project(hydro_params p, int x_start, int slab,
		  float *product_re, float *product_im,
		  fftwf_complex **Tij_then, fftwf_complex **Tij_now) {

  int i, j, l, m;
  int x, y, z;
  float kx, ky, kz;

  int true_x, true_y, true_z;

  float s_x, s_y, s_z;

  for(x=0; x<slab; x++) {
    for(y=0; y<p.Ly; y++) {
      for(z=0; z<p.Lz; z++) {

        s_x = 1.0;
        s_y = 1.0;
        s_z = 1.0;

        if(x+x_start > p.Lx/2) {
          true_x = -(p.Lx - (x+x_start));
          s_x = -1.0;
        } else {
          true_x = x+x_start;
          s_x = 1.0;
        }

        if(y > p.Ly/2) {
          true_y = -(p.Ly - y);
          s_y = -1.0;
        } else {
          true_y = y;
          s_y = 1.0;
        }

        if(z > p.Lz/2) {
          true_z = -(p.Lz - z);
          s_z = -1.0;
        } else {
          true_z = z;
          s_z = 1.0;
        }

	/*
	kx = ((float)(x+x_start))*2.0*M_PI/((float)p.Lx);
	ky = ((float)y)*2.0*M_PI/((float)p.Ly);
	kz = ((float)z)*2.0*M_PI/((float)p.Lz);
	*/
	

        kx = s_x*sqrt((2.0 - 2.0*cos(((float)(true_x))*2.0*M_PI/(((float)p.Lx))))/(p.dx*p.dx));
	ky = s_y*sqrt((2.0 - 2.0*cos(((float)(true_y))*2.0*M_PI/(((float)p.Ly))))/(p.dx*p.dx));
	kz = s_z*sqrt((2.0 - 2.0*cos(((float)(true_z))*2.0*M_PI/(((float)p.Lz))))/(p.dx*p.dx));

	product_re[x*p.Ly*p.Lz + y*p.Lz + z] = 0.0;
	product_im[x*p.Ly*p.Lz + y*p.Lz + z] = 0.0;

	// Sum over indices on lambda_{ij,lm} and udot_{ij}(k)
	for(i=1; i<=3; i++) {
	  for(j=1; j<=3; j++) {
	    for(l=1; l<=3; l++) {
	      for(m=1; m<=3; m++) {

		product_re[x*p.Ly*p.Lz + y*p.Lz + z] +=
		  lambda(i, j, l, m, kx, ky, kz)
		  *(Tij_then[indexof(i,j)][x*p.Ly*p.Lz + y*p.Lz + z][0]
		   *Tij_now[indexof(l,m)][x*p.Ly*p.Lz + y*p.Lz + z][0]
		    + Tij_then[indexof(i,j)][x*p.Ly*p.Lz + y*p.Lz + z][1]
		    *Tij_now[indexof(l,m)][x*p.Ly*p.Lz + y*p.Lz + z][1]);

		product_im[x*p.Ly*p.Lz + y*p.Lz + z] +=
		  lambda(i, j, l, m, kx, ky, kz)
		  *(-Tij_then[indexof(i,j)][x*p.Ly*p.Lz + y*p.Lz + z][0]
		   *Tij_now[indexof(l,m)][x*p.Ly*p.Lz + y*p.Lz + z][1]
		    + Tij_then[indexof(i,j)][x*p.Ly*p.Lz + y*p.Lz + z][1]
		    *Tij_now[indexof(l,m)][x*p.Ly*p.Lz + y*p.Lz + z][0]);
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
 * Actually calculate the uetc at the current timestep and
 * store in an appropriately named file.
 */
void fft_uetc(hydro_fields f, hydro_params p, int step) {

  ptrdiff_t x_thickness, x_start, alloc_local;

  ptrdiff_t n0 = p.Lx;
  ptrdiff_t n1 = p.Ly;
  ptrdiff_t n2 = p.Lz;

  //  int slab;

  MPI_Status status;

  float kmode, modsq;

  int x, y, z;
  int i;

  float fft_norm = (1.0/(((float)p.Lx)*((float)p.Ly)*((float)p.Lz)))
    *(1.0/sqrt(32.0*M_PI));
  // *p.a*p.a*p.a*p.dx*p.dx*p.dx

  float *trim_then = (float *)malloc(p.slicex*p.slicey*p.Lz*sizeof(float));
  float *trim_now = (float *)malloc(p.slicex*p.slicey*p.Lz*sizeof(float));

  printf0(p, "Starting UETC calculation.\n");

  float start = clock();

  fftwf_mpi_init();

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
          //      fprintf(stderr, "therefore node %d is responsible for x=%d\n", i, x);
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
  */

  /////////////// here compute other thing to be squished ////////////////

  float ****Tij = make_tensor(p);
  stress_energy(f, p, Tij);

  printf0(p,"generated stress-energy\n");

  /////////////// done /////////////////


  fftwf_complex *in = fftwf_alloc_complex(alloc_local);
  fftwf_complex *out = fftwf_alloc_complex(alloc_local);

  fftwf_complex **outcpts_then =
    (fftwf_complex **)malloc(6*sizeof(fftwf_complex *));
  fftwf_complex **outcpts_now =
    (fftwf_complex **)malloc(6*sizeof(fftwf_complex *));
  
  for(i=0;i<TENSOR_CPTS;i++) {
    outcpts_then[i] = fftwf_alloc_complex(alloc_local);
    outcpts_now[i] = fftwf_alloc_complex(alloc_local);
  }

  float *slice_then = (float *)malloc(x_thickness*p.Ly*p.Lz*sizeof(float));
  float *slice_now = (float *)malloc(x_thickness*p.Ly*p.Lz*sizeof(float));

  float *slice_re = (float *)malloc(x_thickness*p.Ly*p.Lz*sizeof(float));
  float *slice_im = (float *)malloc(x_thickness*p.Ly*p.Lz*sizeof(float));


  int nx = p.Lx/p.slicex;
  int ny = p.Ly/p.slicey;

  int ry;

  // Now planning  
  fftwf_plan plan = fftwf_mpi_plan_dft_3d(p.Lx, p.Ly, p.Lz,
					in, out, MPI_COMM_WORLD,
					FFTW_BACKWARD, FFTW_ESTIMATE);

  for(i = 0; i < TENSOR_CPTS; i++) {

    float check = 0.0;    
    
    for(x=1; x<=p.slicex; x++) {
      for(y=1; y<=p.slicey; y++) {
	for(z=0; z<p.Lz; z++) {
	  trim_then[(x-1)*p.slicey*p.Lz + (y-1)*p.Lz + z] 
	    = f.initial_Tij[i][x][y][z];
	  trim_now[(x-1)*p.slicey*p.Lz + (y-1)*p.Lz + z]
	    = Tij[i][x][y][z];
	}
      }
    }

    //    fprintf(stderr,"Check: uij %lf\n", reduce_sum(check, p));
    
    for(x=0; x<p.Lx; x++) {
      
      // Check whether we are the target node for this slab
      if(map[x] == p.rank) {
	
	for(ry = 0; ry < ny; ry++) {
	  
	  if((x/p.slicex == p.myposx) && (ry == p.myposy)) {
	    
	    memcpy(&slice_then[(x-x_start)*p.Ly*p.Lz + ry*p.slicey*p.Lz],
		   &trim_then[(x-p.shiftx)*p.slicey*p.Lz],
		   p.slicey*p.Lz*sizeof(float));

	    memcpy(&slice_now[(x-x_start)*p.Ly*p.Lz + ry*p.slicey*p.Lz],
		   &trim_now[(x-p.shiftx)*p.slicey*p.Lz],
		   p.slicey*p.Lz*sizeof(float));
	    continue;
	  }
	  
	  
	  MPI_Recv(&slice_then[(x-x_start)*p.Ly*p.Lz + ry*p.slicey*p.Lz],
		   p.slicey*p.Lz,
		   MPI_FLOAT,
		   ry*nx + x/p.slicex,
		   x*ny + ry,
		   MPI_COMM_WORLD,
		   &status);

	  MPI_Recv(&slice_now[(x-x_start)*p.Ly*p.Lz + ry*p.slicey*p.Lz],
		   p.slicey*p.Lz,
		   MPI_FLOAT,
		   ry*nx + x/p.slicex,
		   p.Lx*ny + x*ny + ry,
		   MPI_COMM_WORLD,
		   &status);
	}
	
	
      } else if((x >= p.shiftx) && (x < (p.shiftx + p.slicex))) {
	
	
      MPI_Send(&trim_then[(x-p.shiftx)*p.slicey*p.Lz],
	       p.slicey*p.Lz,
	       MPI_FLOAT,
	       map[x],
	       x*ny + p.myposy,
	       MPI_COMM_WORLD);

      MPI_Send(&trim_now[(x-p.shiftx)*p.slicey*p.Lz],
	       p.slicey*p.Lz,
	       MPI_FLOAT,
	       map[x],
	       p.Lx*ny + x*ny + p.myposy,
	       MPI_COMM_WORLD);
      }
      
    }

    float total = 0.0;

    for(x=0;x<x_thickness;x++) {
      for(y=0;y<p.Ly;y++) {
	for(z=0;z<p.Lz;z++) {
	  in[x*p.Ly*p.Lz + y*p.Lz + z][0] = slice_then[x*p.Ly*p.Lz + y*p.Lz + z];
	  in[x*p.Ly*p.Lz + y*p.Lz + z][1] = 0.0;
	}
      }
    }

    printf0(p, "executing FFT for early tensor, cpt %d/%d\n", i, TENSOR_CPTS);

    fftwf_mpi_execute_dft(plan, in, out);

    for(x=0;x<x_thickness;x++) {
      for(y=0;y<p.Ly;y++) {
	for(z=0;z<p.Lz;z++) {
	  outcpts_then[i][x*p.Ly*p.Lz + y*p.Lz + z][0] 
	    = fft_norm*out[x*p.Ly*p.Lz + y*p.Lz + z][0];
	  outcpts_then[i][x*p.Ly*p.Lz + y*p.Lz + z][1]
	    = fft_norm*out[x*p.Ly*p.Lz + y*p.Lz + z][1];

	}
      }
    }



    for(x=0;x<x_thickness;x++) {
      for(y=0;y<p.Ly;y++) {
	for(z=0;z<p.Lz;z++) {
	  in[x*p.Ly*p.Lz + y*p.Lz + z][0]
	    = slice_now[x*p.Ly*p.Lz + y*p.Lz + z];
	  in[x*p.Ly*p.Lz + y*p.Lz + z][1]
	    = 0.0;
	}
      }
    }

    printf0(p, "executing FFT for late tensor, cpt %d/%d\n", i, TENSOR_CPTS);

    fftwf_mpi_execute_dft(plan, in, out);

    for(x=0;x<x_thickness;x++) {
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

  // could make strictly > because equal time is just shear stress
  if(step >= p.uetcstart) {

    printf0(p,"projecting\n");

    uetc_project(p, x_start, x_thickness, slice_re, slice_im,
		 outcpts_then, outcpts_now);


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
	

	  kmode = sqrt(
		       ((float)(true_x*true_x))/((float)(p.Lx*p.Lx))
		       + ((float)(true_y*true_y))/((float)(p.Ly*p.Ly))
		       + ((float)(true_z*true_z))/((float)(p.Lz*p.Lz))
		       )*2.0*M_PI;
	  
	  whichbin = (int)floor(kmode/dk);
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
  
      sprintf(fftdest,"uetc-step%d.txt", step);
      
      FILE *fp = fopen(fftdest,"w");
  
      for(i=0;i<nbins;i++) {

	fprintf(fp, "%lf %g %g %d %d\n",
		thisk, bins_re[i], bins_im[i], counts[i], step);

	thisk = thisk + dk;
      }

      fclose(fp);

    }

    free(bins_re);
    free(bins_im);
    free(counts);
  }

  if(step >= p.uetcstart) {


    printf0(p,"projecting for equal time\n");

    uetc_project(p, x_start, x_thickness, slice_re, slice_im,
		 outcpts_now, outcpts_now);


    printf0(p,"binning for equal time\n");

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
	

	  kmode = sqrt(
		       ((float)(true_x*true_x))/((float)(p.Lx*p.Lx))
		       + ((float)(true_y*true_y))/((float)(p.Ly*p.Ly))
		       + ((float)(true_z*true_z))/((float)(p.Lz*p.Lz))
		       )*2.0*M_PI;
	  
	  whichbin = (int)floor(kmode/dk);
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
  
      sprintf(fftdest,"shearstress-step%d.txt", step);
      
      FILE *fp = fopen(fftdest,"w");
  
      for(i=0;i<nbins;i++) {

	fprintf(fp, "%lf %g %g %d %d\n",
		thisk, bins_re[i], bins_im[i], counts[i], step);

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
  free(slice_then);
  free(slice_now);
  free(trim_then);
  free(trim_now);

  fftwf_destroy_plan(plan);
  
  fftwf_free(in);
  fftwf_free(out);

  for(i=0;i<TENSOR_CPTS;i++) {
    fftwf_free(outcpts_then[i]);
    fftwf_free(outcpts_now[i]);
  }


  free_tensor(p, Tij);

  fftwf_mpi_cleanup();

  float end = clock();

  printf0(p, "UETC FFT calculation took %lf\n",
	  ((float) (end - start)) / CLOCKS_PER_SEC);

}




#endif // FFT
