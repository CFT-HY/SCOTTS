/* fft.c
 *
 * Fourier transform (and power spectrum) of a field
 */
#include "hydro.h"


#ifdef INITPS

/* get_normal
 *
 * Lazy and inefficient way of getting gaussian
 * random number
 */
double get_normal(double mean, double dev) {
  double u1, u2;
  
  u1 = drand48();
  u2 = drand48();

  return dev*(sqrt(-2.0*log(u1))*sin(2.0*M_PI*u2)-mean);
}

/* spectrum
 * 
 */
void spectrum(double ksq, hydro_params p, fftw_complex *res) 
{ 

  double phase;

  double L = p.Lx;


  if(fabs(ksq) < 0.000001) {
    (*res)[0] = 0.0;
    (*res)[1] = 0.0;
  } else {

    (*res)[0] = get_normal(0.0, 0.1*exp(-0.25*sqrt(ksq)*((double)L)));
    //    (*res)[1] = get_normal(0.0, sqrt(0.5*hbar/sqrt(ksq + m*m)));
    // Anders said to use separate Gaussian values for each cpt
    /*
    if(fabs(ksq - 0.585786) < 0.0001) {
      printf("%g\n", (*res)[0]);
    }
    */

    (*res)[1] = (*res)[0];

    phase = drand48();

    // Random phase
    (*res)[0] = (*res)[0]*cos(2.0*M_PI*phase);
    (*res)[1] = (*res)[1]*sin(2.0*M_PI*phase);

  } 
}



/* init_ps(hydro_fields f, hydro_params p, double ***field)
 *
 * Initialises the field with a power spectrum.
 * Used to convince reluctant referees.
 */
void init_ps(hydro_fields f, hydro_params p, double ****field) {

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

  //  fftw_mpi_init();

  double checken = 0.0;




  //  alloc_local = fftw_mpi_local_size_3d(n0, n1, n2,
  //				       MPI_COMM_WORLD, &x_thickness, &x_start);

  /*
  slab = (int) x_thickness;

  if(((int)x_thickness) != p.Lx/p.size) {
    fprintf(stderr,
	    "Giving up in FFT: FFTW told me to use a silly layout!\n");
    die(-42);
  }
  */

  /*
  if(((int)x_thickness) == 0) {
    fprintf(stderr, "Giving up in FFT: dx=0!\n");
    die(-43);
  }
  */

  //  fftw_complex *in = fftw_alloc_complex(alloc_local);
  //  fftw_complex *swap_in = fftw_alloc_complex(alloc_local);
  //  fftw_complex *out = fftw_alloc_complex(alloc_local);

  fftw_complex **in = (fftw_complex **)malloc(3*sizeof(fftw_complex *));
  fftw_complex **out = (fftw_complex **)malloc(3*sizeof(fftw_complex *));

  in[0] = fftw_alloc_complex(p.Lx*p.Ly*p.Lz);
  out[0] = fftw_alloc_complex(p.Lx*p.Ly*p.Lz);

  in[1] = fftw_alloc_complex(p.Lx*p.Ly*p.Lz);
  out[1] = fftw_alloc_complex(p.Lx*p.Ly*p.Lz);

  in[2] = fftw_alloc_complex(p.Lx*p.Ly*p.Lz);
  out[2] = fftw_alloc_complex(p.Lx*p.Ly*p.Lz);


  // double *slice = (double *)malloc(slab*p.Ly*p.Lz*sizeof(double));

  int nx = p.Lx/p.slicex;
  int ny = p.Ly/p.slicey;

  int ry;


  /*
   * Set up modes with desired power spectrum...
   */
  int true_x;

  double ksq;

  //  int i;


  for(x=0;x<p.Lx;x++) {
    for(y=0;y<p.Ly;y++) {
      for(z=0;z<p.Lz;z++) {
	//	true_x = x + x_start;

	ksq = (2.0 - 2.0*cos(((double)x)*2.0*M_PI/(((double)p.Lx))))/(p.dx*p.dx);
	ksq += (2.0 - 2.0*cos(((double)y)*2.0*M_PI/(((double)p.Ly))))/(p.dx*p.dx);
	ksq += (2.0 - 2.0*cos(((double)z)*2.0*M_PI/(((double)p.Lz))))/(p.dx*p.dx);




	for(i=0; i<3; i++) {
	  spectrum(ksq, p, &(in[i][x*p.Ly*p.Lz + y*p.Lz + z]));
	  //	  	  in[i][x*p.Ly*p.Lz + y*p.Lz + z][0] = 0.0;
	  //	  	  in[i][x*p.Ly*p.Lz + y*p.Lz + z][1] = 0.0;
	}

	//		in[1][x*p.Ly*p.Lz + y*p.Lz + z][0] = 1.0;
      }
    }
  }

  // and swap them around
  /*
  if(p.size - p.rank -1 == p.rank) {
    memcpy(swap_in, in, x_thickness*p.Ly*p.Lz*sizeof(fftw_complex));
  } else {

    MPI_Send(in,
	     2*x_thickness*p.Ly*p.Lz,
	     MPI_DOUBLE,
	     p.size - p.rank -1,
	     x_start,
	     MPI_COMM_WORLD);

    MPI_Recv(swap_in,
	     2*x_thickness*p.Ly*p.Lz,
	     MPI_DOUBLE,
	     p.size-p.rank-1,
	     x_start,
	     MPI_COMM_WORLD,
	     &status);
  }

//  fprintf(stderr,"x_thickness is %d\n", x_thickness);
  */

  // moduli will not work in parallel, also not swapped
  for(x=0;x<p.Lx;x++) {
    for(y=0;y<p.Ly;y++) {
      for(z=0;z<p.Lz;z++) {
	for(i=0; i<3; i++) {
	  in[i][x*p.Ly*p.Lz + y*p.Lz + z][0] = in[i][((p.Lx-x)%p.Lx)*p.Ly*p.Lz
					       + ((p.Ly-y)%p.Ly)*p.Lz + ((p.Lz-z)%p.Lz)][0];
	  in[i][x*p.Ly*p.Lz + y*p.Lz + z][1] = -1.0*in[i][((p.Lx-x)%p.Lx)*p.Ly*p.Lz
						    + ((p.Ly-y)%p.Ly)*p.Lz + ((p.Lz-z)%p.Lz)][1];

	  //	true_x = x;
	  
	  if((x == 0 || x == p.Lx/2) && (y == 0 || y == p.Ly/2) && (z == 0 || z == p.Lz/2)) {
	    in[i][x*p.Ly*p.Lz + y*p.Lz + z][0] = sqrt(in[i][((p.Lx-x)%p.Lx)*p.Ly*p.Lz
						      + ((p.Ly-y)%p.Ly)*p.Lz + ((p.Lz-z)%p.Lz)][0]
						   *in[i][((p.Lx-x)%p.Lx)*p.Ly*p.Lz
						       + ((p.Ly-y)%p.Ly)*p.Lz + ((p.Lz-z)%p.Lz)][0]
						   + in[i][((p.Lx-x)%p.Lx)*p.Ly*p.Lz
							+ ((p.Ly-y)%p.Ly)*p.Lz + ((p.Lz-z)%p.Lz)][1]
						   *in[i][((p.Lx-x)%p.Lx)*p.Ly*p.Lz
						       + ((p.Ly-y)%p.Ly)*p.Lz + ((p.Lz-z)%p.Lz)][1]);
	    in[i][x*p.Ly*p.Lz + y*p.Lz + z][1] = 0.0;
	  }
	}
      }
    }
  }



  double kx, ky, kz;
  int j;

  double in_proj_re[3];
  double in_proj_im[3];


  int reruns = 0;

  double resid, stuff;


  for(reruns = 0; reruns < 5; reruns++) {
    resid = 0.0;
    stuff = 0.0;
    
    // Project out divergenceless bit!
    for(x=0;x<p.Lx;x++) {
      for(y=0;y<p.Ly;y++) {
	for(z=0;z<p.Lz;z++) {



	  
	  kx = sqrt((2.0 - 2.0*cos(((double)x)*2.0*M_PI/(((double)p.Lx))))/(p.dx*p.dx));
	  // ((double)x)*2.0*M_PI/((double)p.Lx);
	  ky = sqrt((2.0 - 2.0*cos(((double)y)*2.0*M_PI/(((double)p.Lx))))/(p.dx*p.dx));
	  // ((double)y)*2.0*M_PI/((double)p.Ly);
	  kz = sqrt((2.0 - 2.0*cos(((double)z)*2.0*M_PI/(((double)p.Lx))))/(p.dx*p.dx));
	  // ((double)z)*2.0*M_PI/((double)p.Lz);
	  
	  in_proj_re[0] = 0.0;
	  in_proj_re[1] = 0.0;
	  in_proj_re[2] = 0.0;
	  
	  in_proj_im[0] = 0.0;
	  in_proj_im[1] = 0.0;
	  in_proj_im[2] = 0.0;

	  //	  if(x == 1 && y == 1 && z == 1 && !p.rank)
	  //	     fprintf(stderr,"site before (%g,%g,%g)\n", in[0][x*p.Ly*p.Lz + y*p.Lz + z][0] , in[1][x*p.Ly*p.Lz + y*p.Lz + z][0], in[2][x*p.Ly*p.Lz + y*p.Lz + z][0]);
	  
	  //	  if(x == 1 && y == 1 && z == 1 && !p.rank)
	  //	    fprintf(stderr,"vel proj:\n");
	  
	  for(i=1; i<=3; i++) {

	    double res_re = 0.0;
	    double res_im = 0.0;

	    for(j=1; j<=3; j++) {
	      
	      in_proj_re[i-1] += vel_proj(i*10 + j, kx, ky, kz)*in[j-1][x*p.Ly*p.Lz + y*p.Lz + z][0];
	      in_proj_im[i-1] += vel_proj(i*10 + j, kx, ky, kz)*in[j-1][x*p.Ly*p.Lz + y*p.Lz + z][1];

	      // if(x == 1 && y == 1 && z == 1 && !p.rank)
	      // fprintf(stderr,"%g ",vel_proj(i*10 + j, kx, ky, kz));
	      // if(x == 1 && y == 1 && z == 1 && !p.rank)
	      // fprintf(stderr,"entry %d: %g += %g*%g\n", i, in_proj_re[i-1],
	      // vel_proj(i*10 + j, kx, ky, kz), in[j-1][x*p.Ly*p.Lz + y*p.Lz + z][0]);



	      if(i == j) {
		res_re += (1.0-vel_proj(i*10 + j, kx, ky, kz))
		  *in[j-1][x*p.Ly*p.Lz + y*p.Lz + z][0];
		res_im += (1.0-vel_proj(i*10 + j, kx, ky, kz))
		  *in[j-1][x*p.Ly*p.Lz + y*p.Lz + z][1];
	      } else {
		res_re += (-1.0*vel_proj(i*10 + j, kx, ky, kz))
		  *in[j-1][x*p.Ly*p.Lz + y*p.Lz + z][0];
		res_im += (-1.0*vel_proj(i*10 + j, kx, ky, kz))
		  *in[j-1][x*p.Ly*p.Lz + y*p.Lz + z][1];
	      }

	      
			
	    }
	      stuff += in_proj_re[i-1]*in_proj_re[i-1] + in_proj_im[i-1]*in_proj_im[i-1];
	      resid += res_re*res_re + res_im*res_im;


	    //	      if(x == 1 && y == 1 && z == 1 && !p.rank)
	    //		fprintf(stderr,"\n");
	  }


	  in[0][x*p.Ly*p.Lz + y*p.Lz + z][0] = in_proj_re[0];
	  in[1][x*p.Ly*p.Lz + y*p.Lz + z][0] = in_proj_re[1];
	  in[2][x*p.Ly*p.Lz + y*p.Lz + z][0] = in_proj_re[2];

	  in[0][x*p.Ly*p.Lz + y*p.Lz + z][1] = in_proj_im[0];
	  in[1][x*p.Ly*p.Lz + y*p.Lz + z][1] = in_proj_im[1];
	  in[2][x*p.Ly*p.Lz + y*p.Lz + z][1] = in_proj_im[2];

	  //	  if(x == 1 && y == 1 && z == 1 && !p.rank)
	  //	     fprintf(stderr,"site now (%g,%g,%g)\n", in[0][x*p.Ly*p.Lz + y*p.Lz + z][0] , in[1][x*p.Ly*p.Lz + y*p.Lz + z][0], in[2][x*p.Ly*p.Lz + y*p.Lz + z][0]);


	}
      }
    }

    if(!p.rank)
      fprintf(stderr,"run %d, stuff is %g, resid is now %g\n", reruns, stuff, resid);
  }

  // Now planning  
  fftw_plan plan = fftw_plan_dft_3d(p.Lx, p.Ly, p.Lz,
				    in[0], out[0],
				    FFTW_FORWARD, FFTW_ESTIMATE);


  
  fftw_execute_dft(plan, in[0], out[0]);
  fftw_execute_dft(plan, in[1], out[1]);
  fftw_execute_dft(plan, in[2], out[2]);


  /*
   * Now copy out stuff from slab arrays into pencils. Reverse of what we usually do.
   */



  for(x=1;x<=p.slicex;x++) {
    for(y=1;y<=p.slicey;y++) {
      for(z=0;z<p.Lz;z++) {	
	for(i=0;i<3;i++) {
	  field[i][x][y][z] = out[i][(x+p.shiftx-1)*p.Ly*p.Lz + (y+p.shifty-1)*p.Lz + z][0];


	/*	
	printf("%d %d %g %g\n", p.rank,
	       (x+p.shiftx-1)*p.Ly*p.Lz + (y+p.shifty-1)*p.Lz + z,
	       out[(x+p.shiftx-1)*p.Ly*p.Lz + (y+p.shifty-1)*p.Lz + z][0],
	       out[(x+p.shiftx-1)*p.Ly*p.Lz + (y+p.shifty-1)*p.Lz + z][1]);
	*/
	//	slice[x*p.Ly*p.Lz + y*p.Lz + z] = out[x*p.Ly*p.Lz + y*p.Lz + z][0];

	}
      }
    }
  }

  halo_field(field[0], p);
  halo_field(field[1], p);
  halo_field(field[2], p);

  //  fprintf(stderr,"111 = %g\n", field[1][1][1]);

  //  exit(0);

  /*
  for(x=0; x<p.Lx; x++) {

    // Check whether we are the target node for this slab
    if((x >= x_start) && ( x < (x_start + x_thickness))) {

      for(ry = 0; ry < ny; ry++) {

	if((x/p.slicex == p.myposx) && (ry == p.myposy)) {

	  memcpy(&trim_field[(x-p.shiftx)*p.slicey*p.Lz],
		 &slice[(x-x_start)*p.Ly*p.Lz + ry*p.slicey*p.Lz],
		 p.slicey*p.Lz*sizeof(double));
	  continue;
	}

	MPI_Send(&slice[(x-x_start)*p.Ly*p.Lz + ry*p.slicey*p.Lz],
		 p.slicey*p.Lz,
		 MPI_DOUBLE,
		 ry*nx + x/p.slicex,
		 x*ny + ry,
		 MPI_COMM_WORLD);
      }

    } else if((x >= p.shiftx) && (x < (p.shiftx + p.slicex))) {

      MPI_Recv(&trim_field[(x-p.shiftx)*p.slicey*p.Lz],
	       p.slicey*p.Lz,
	       MPI_DOUBLE,
	       x/slab,
	       x*ny + p.myposy,
	       MPI_COMM_WORLD,
	       &status);
    }

    // then wait until that slice of x is sorted
    
  }
  // Alloc slab array and do communication to get it into place

  for(x=1; x<=p.slicex; x++) {
    for(y=1; y<=p.slicey; y++) {
      for(z=0; z<p.Lz; z++) {

	field[x][y][z] = trim_field[(x-1)*p.slicey*p.Lz + (y-1)*p.Lz + z];
      }
    }
  }  


  */





  //  free(slice);

  fftw_destroy_plan(plan);
  
  fftw_free(in[0]);
  fftw_free(out[0]);
  fftw_free(in[1]);
  fftw_free(out[1]);
  fftw_free(in[2]);
  fftw_free(out[2]);

  free(in);
  free(out);

  free(trim_field);

  //  fftw_mpi_cleanup();

}






#endif // INITPS
