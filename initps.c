/* fft.c
 *
 * Fourier transform (and power spectrum) of a field
 */
#include "hydro.h"


#ifdef INITPS


float get_momtot(hydro_fields f, hydro_params p) {

  int x, y, z, xmax;


  float momtot = 0.0;

  // Just search for maxmimum
  for(x = 1; x <= p.slicex; x++) {
    for(y = 1; y <= p.slicey; y++) {
      for(z = 0; z < p.Lz; z++) {
        momtot += sqrt(f.Z[0][x][y][z]*f.Z[0][x][y][z]
          + f.Z[1][x][y][z]*f.Z[1][x][y][z]
                       + f.Z[2][x][y][z]*f.Z[2][x][y][z]);

      }
    }
  }

  return momtot;
}

void norm_power(hydro_fields f, hydro_params p, float ****field) {

  int i, x, y, z;
  // should be normalised to volume, this makes it okay
  float momtot = reduce_sum(get_momtot(f, p),p)/((float)p.Lx*p.Ly*p.Lz);

  printf0(p, "momtot per site is %g\n", momtot);
  for(i=0; i<3; i++) {
    for(x=1; x<=p.slicex; x++) {
      for(y=1; y<=p.slicey; y++) {
        for(z=0; z<p.Lz; z++) {
          field[i][x][y][z] = p.initps*field[i][x][y][z]/momtot;
        }
      }
    }
    halo_field(field[i], p);
  }
}

/* project_down
 *
 * Project out rotational or divergent bits
 */
void project_down(hydro_params p, fftwf_complex **in, int shift_x, int x_thickness, int times) {

  int x, y, z;
  float kx, ky, kz;

  float in_proj_re[3];
  float in_proj_im[3];

  float res_re, res_im;

  int i, j;

  float resid, stuff;

  float tdx = 0.5;

  int reruns;

  int true_x, true_y, true_z;

  float s_x, s_y, s_z;

  for(reruns = 0; reruns < times; reruns++) {
    resid = 0.0;
    stuff = 0.0;
    
    // Project out divergenceless bit!
    for(x=0;x<x_thickness;x++) {
      for(y=0;y<p.Ly;y++) {
	for(z=0;z<p.Lz;z++) {

	  s_x = 1.0;
	  s_y = 1.0;
	  s_z = 1.0;

	  if(x+shift_x > p.Lx/2) {
	    true_x = -(p.Lx - (x+shift_x));
	    s_x = -1.0;
	  } else {
	    true_x = x+shift_x;
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
	  kx = sqrt((2.0 - 2.0*cos(((float)(x+shift_x))*2.0*M_PI/(((float)p.Lx))))/(tdx*tdx));
	  ky = sqrt((2.0 - 2.0*cos(((float)y)*2.0*M_PI/(((float)p.Ly))))/(tdx*tdx));
	  kz = sqrt((2.0 - 2.0*cos(((float)z)*2.0*M_PI/(((float)p.Lz))))/(tdx*tdx));
	  */

	  kx = s_x*sqrt((2.0 - 2.0*cos(((float)(true_x))*2.0*M_PI/(((float)p.Lx))))/(p.dx*p.dx));
	  ky = s_y*sqrt((2.0 - 2.0*cos(((float)(true_y))*2.0*M_PI/(((float)p.Ly))))/(p.dx*p.dx));
	  kz = s_z*sqrt((2.0 - 2.0*cos(((float)(true_z))*2.0*M_PI/(((float)p.Lz))))/(p.dx*p.dx));
	  
	  in_proj_re[0] = 0.0;
	  in_proj_re[1] = 0.0;
	  in_proj_re[2] = 0.0;
	  
	  in_proj_im[0] = 0.0;
	  in_proj_im[1] = 0.0;
	  in_proj_im[2] = 0.0;

	  
	  for(i=1; i<=3; i++) {

	    res_re = 0.0;
	    res_im = 0.0;

	    for(j=1; j<=3; j++) {

	      // #ifndef DIVPS	      
	      /* Rot?
	      in_proj_re[i-1] += vel_proj(i*10 + j, kx, ky, kz)*in[j-1][x*p.Ly*p.Lz + y*p.Lz + z][0];
	      in_proj_im[i-1] += vel_proj(i*10 + j, kx, ky, kz)*in[j-1][x*p.Ly*p.Lz + y*p.Lz + z][1];
	      */
	      // #else
	      // #warning DIVPS enabled
	      if(i == j) {
		in_proj_re[i-1] += (1.0-vel_proj(i*10 + j, kx, ky, kz))*in[j-1][x*p.Ly*p.Lz + y*p.Lz + z][0];
		in_proj_im[i-1] += (1.0-vel_proj(i*10 + j, kx, ky, kz))*in[j-1][x*p.Ly*p.Lz + y*p.Lz + z][1];
	      } else {
		in_proj_re[i-1] += (-1.0*vel_proj(i*10 + j, kx, ky, kz))*in[j-1][x*p.Ly*p.Lz + y*p.Lz + z][0];
		in_proj_im[i-1] += (-1.0*vel_proj(i*10 + j, kx, ky, kz))*in[j-1][x*p.Ly*p.Lz + y*p.Lz + z][1];
	      }
	      // #endif




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



	  }


	  in[0][x*p.Ly*p.Lz + y*p.Lz + z][0] = in_proj_re[0];
	  in[1][x*p.Ly*p.Lz + y*p.Lz + z][0] = in_proj_re[1];
	  in[2][x*p.Ly*p.Lz + y*p.Lz + z][0] = in_proj_re[2];

	  in[0][x*p.Ly*p.Lz + y*p.Lz + z][1] = in_proj_im[0];
	  in[1][x*p.Ly*p.Lz + y*p.Lz + z][1] = in_proj_im[1];
	  in[2][x*p.Ly*p.Lz + y*p.Lz + z][1] = in_proj_im[2];

	}
      }
    }

    if(!p.rank)
      fprintf(stderr,"run %d, stuff is %g, resid is now %g\n", reruns, stuff, resid);
  }

}


/* get_normal
 *
 * Lazy and inefficient way of getting gaussian
 * random number
 */
float get_normal(float mean, float dev) {
  float u1, u2;
  
  u1 = drand48();
  u2 = drand48();

  return dev*(sqrt(-2.0*log(u1))*sin(2.0*M_PI*u2)-mean);
}

/* spectrum
 * 
 */
void spectrum(float ksq, hydro_params p, fftwf_complex *res) 
{ 

  float phase;

  float L = p.Lx;


  if(fabs(ksq) < 0.000001 || fabs(ksq) < p.initcutoff) {
    (*res)[0] = 0.0;
    (*res)[1] = 0.0;
  } else {


    //  (*res)[0] = p.initps*exp(-0.25*sqrt(ksq)*p.dx*((float)L));    
    //    (*res)[0] = get_normal(0.0, p.initps*exp(-0.25*sqrt(ksq)*p.dx*((float)L)));
    (*res)[0] = get_normal(0.0, exp(-1.0*sqrt(ksq)*p.dx*p.initlength));


    (*res)[1] = (*res)[0];

    phase = drand48();

    // Random phase
    (*res)[0] = (*res)[0]*cos(2.0*M_PI*phase);
    (*res)[1] = (*res)[1]*sin(2.0*M_PI*phase);

  } 
}








/* init_ps(hydro_fields f, hydro_params p, float ***field)
 *
 * Initialises the field with a power spectrum.
 * Used to convince reluctant referees.
 */
void init_ps(hydro_fields f, hydro_params p, float ****field) {


  ptrdiff_t x_thickness, x_start, alloc_local;


  ptrdiff_t n0 = p.Lx;
  ptrdiff_t n1 = p.Ly;
  ptrdiff_t n2 = p.Lz;

  int slab;

  MPI_Status status;

  float kmode, modsq;

  int x, y, z;
  int i, j;

  float *trim = (float *)malloc(p.slicex*p.slicey*p.Lz*sizeof(float));


  fftwf_mpi_init();

  alloc_local = fftwf_mpi_local_size_3d(n0, n1, n2,
                                       MPI_COMM_WORLD,
                                       &x_thickness, &x_start);



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


  fftwf_complex **in = (fftwf_complex **)malloc(3*sizeof(fftwf_complex *));
  fftwf_complex **swap_in = (fftwf_complex **)malloc(3*sizeof(fftwf_complex *));

  fftwf_complex **out = (fftwf_complex **)malloc(3*sizeof(fftwf_complex *));

  in[0] = fftwf_alloc_complex(x_thickness*p.Ly*p.Lz);
  swap_in[0] = fftwf_alloc_complex(x_thickness*p.Ly*p.Lz);

  out[0] = fftwf_alloc_complex(x_thickness*p.Ly*p.Lz);

  in[1] = fftwf_alloc_complex(x_thickness*p.Ly*p.Lz);
  swap_in[1] = fftwf_alloc_complex(x_thickness*p.Ly*p.Lz);

  out[1] = fftwf_alloc_complex(x_thickness*p.Ly*p.Lz);

  in[2] = fftwf_alloc_complex(x_thickness*p.Ly*p.Lz);
  swap_in[2] = fftwf_alloc_complex(x_thickness*p.Ly*p.Lz);

  out[2] = fftwf_alloc_complex(x_thickness*p.Ly*p.Lz);



  int nx = p.Lx/p.slicex;
  int ny = p.Ly/p.slicey;

  int ry;


  /*
   * Set up modes with desired power spectrum...
   */

  float ksq;

  //  int i;

  int true_x, true_y, true_z;


  for(x=x_start;x<x_start+x_thickness;x++) {
    for(y=0;y<p.Ly;y++) {
      for(z=0;z<p.Lz;z++) {

	if(x > p.Lx/2) {
	  true_x = -(p.Lx - x);
	} else {
	  true_x = x;
	}

	if(y > p.Ly/2) {
	  true_y = -(p.Ly - y);
	} else {
	  true_y = y;
	}

	if(z > p.Lz/2) {
	  true_z = -(p.Lz - z);
	} else {
	  true_z = z;
	}

	ksq = (2.0 - 2.0*cos(((float)true_x)*2.0*M_PI/(((float)p.Lx))))/(p.dx*p.dx);
	ksq += (2.0 - 2.0*cos(((float)true_y)*2.0*M_PI/(((float)p.Ly))))/(p.dx*p.dx);
	ksq += (2.0 - 2.0*cos(((float)true_z)*2.0*M_PI/(((float)p.Lz))))/(p.dx*p.dx);

	for(i=0; i<3; i++) {
	  spectrum(ksq, p, &(in[i][(x-x_start)*p.Ly*p.Lz + y*p.Lz + z]));
	}

      }
    }
  }

  project_down(p, in, x_start, x_thickness, 5);

  /*
  for(x=x_start; x<(x_start+x_thickness); x++) {
    fprintf(stderr, "** node %d, x %d:\t%10.3g %10.3g\t **\n", p.rank, x, in[0][(x-x_start)*p.Ly*p.Lz][0], in[0][(x-x_start)*p.Ly*p.Lz][1]);
  } 
  */

  printf0(p, "Local spectrum initialisation done\n");

  // At this point each node knows what it ought to have, but we need to 'hermitianise' it
  // Need to do swapping


  // These are the minimum and maximum (INCLUSIVE) x-slice coordinates needed by this node
  int min_needed = p.Lx-(x_start+x_thickness)+1;
  int max_needed = p.Lx-x_start;

  for(i=0; i<3; i++) {

      
    //	fprintf(stderr, "Node %d  for %d-%d needs %d-%d\n", p.rank, x_start, x_start+x_thickness-1,  min_needed, max_needed);


	for(j=0; j<=p.Lx; j++) {

	  if(j>= x_start && j < x_start+x_thickness) {


	    //	    fprintf(stderr, "Rank %d: SEND row %02d to node %d (rank %d, size %d)", p.rank, j, ((p.Lx-j))/x_thickness % p.size, p.rank, p.size);

	    if(((p.Lx-j))/x_thickness % p.size == p.rank) {
	      // copy - done on recv side
	      //	      fprintf(stderr, " (c)\n");
	    } else {
	      //	      fprintf(stderr, " (m)\n");
	      MPI_Send(&(((float *)in[i])[((j-x_start)%p.Lx)*p.Ly*p.Lz*2]), 2*p.Ly*p.Lz, MPI_FLOAT, ((p.Lx-j))/x_thickness % p.size, j, MPI_COMM_WORLD);
	    }
	  }


	  if(j>=min_needed && j<=max_needed) {
	    //	    fprintf(stderr, "Rank %d: RECV row %02d from node %d and becomes swap row %d", p.rank, j % p.Lx, (j % p.Lx)/x_thickness, j-min_needed);

	    if((j % p.Lx)/x_thickness == p.rank) {
	      //	      fprintf(stderr, " (c)\n");
	      // copy - need to check
	      //	      fprintf(stderr,"rank %d copying row %d, (local row %d) %g %g to self\n", p.rank, j, (j-x_start)%p.Lx, in[i][((j-x_start)%p.Lx)*p.Ly*p.Lz][0], in[i][((j-x_start)%p.Lx)*p.Ly*p.Lz][1]);
	      memcpy(&(((fftwf_complex *)swap_in[i])[(j-min_needed)*p.Ly*p.Lz]), &(((fftwf_complex *)in[i])[((j-x_start)%p.Lx)*p.Ly*p.Lz]), p.Ly*p.Lz*sizeof(fftwf_complex));

	      //	      fprintf(stderr, "rank %d copy is %g %g\n", p.rank, swap_in[i][(j-min_needed)*p.Ly*p.Lz][0], swap_in[i][(j-min_needed)*p.Ly*p.Lz][1]);
	    } else {	     
	      //	      fprintf(stderr, " (m)\n");
	      MPI_Recv(&(((float *)swap_in[i])[(j-min_needed)*p.Ly*p.Lz*2]), 2*p.Ly*p.Lz, MPI_FLOAT, (j % p.Lx)/x_thickness, j, MPI_COMM_WORLD, &status);
	    }
	  }



    }
      
  }

  // at this point swap_in has the stuff needed to do the conjugate properly

  //  MPI_Barrier(MPI_COMM_WORLD);

  // Debug messages
  //  printf0(p, "Done swap_in initialisation\n", in[0][4][0], in[0][4][1], swap_in[0][4][0], swap_in[0][4][1]);
  //  int swap_x_start = p.Lx - x_start - x_thickness;
  //  printf0(p, "Swap_x_start is %d, x_start is %d, x_thickness is %d\n", swap_x_start, x_start, x_thickness);


  for(x=x_start; x < x_start+x_thickness; x++) {

    /*
    // Debug message
    fprintf(stderr, "Rank %d: USE Row %d has dual row %d which is swap row %d\n",
	    p.rank, x, (p.Lx-x)%p.Lx,
	    ((((p.Lx-x)  )% x_thickness) + x_thickness - 1) % x_thickness );
    */

    // This is the row in the swap_in where one can find (p.Lx-x)
    int swap_row = (((p.Lx-x) % x_thickness) + x_thickness - 1) % x_thickness;

    for(y=0;y<p.Ly;y++) {
      for(z=0;z<p.Lz;z++) {
	for(i=0; i<3; i++) {

	  // Conjugate the system. Treat special cases first.

	  if((x == 0 || x == p.Lx/2) && (y == 0 || y == p.Ly/2) && (z == 0 || z == p.Lz/2)) {
	    

	    /* These sites need to be pure real.
	     */	    
	    in[i][(x-x_start)*p.Ly*p.Lz + y*p.Lz + z][0] = sqrt(swap_in[i][swap_row*p.Ly*p.Lz
						      + ((p.Ly-y)%p.Ly)*p.Lz + ((p.Lz-z)%p.Lz)][0]
						      *swap_in[i][swap_row*p.Ly*p.Lz
						       + ((p.Ly-y)%p.Ly)*p.Lz + ((p.Lz-z)%p.Lz)][0]
						      + swap_in[i][swap_row*p.Ly*p.Lz
							+ ((p.Ly-y)%p.Ly)*p.Lz + ((p.Lz-z)%p.Lz)][1]
						      *swap_in[i][swap_row*p.Ly*p.Lz
						      + ((p.Ly-y)%p.Ly)*p.Lz + ((p.Lz-z)%p.Lz)][1]);
	  
	    in[i][(x-x_start)*p.Ly*p.Lz + y*p.Lz + z][1] = 0.0;
	    

	     } 

	  else if((x == 0 || x == p.Lx/2)) {



	    /* need local mirror image on these slices, so ignore swap rows
	     * if we used swap_in we'd not end up with something that is hermitian because otherwise
	     * [<>          []] initially would be swapped (and conjugated) to become
	     * [[]          <>]
	     */
	    in[i][(x-x_start)*p.Ly*p.Lz + y*p.Lz + z][0] = in[i][(x-x_start)*p.Ly*p.Lz
							 + ((p.Ly-y)%p.Ly)*p.Lz + ((p.Lz-z)%p.Lz)][0];
	    in[i][(x-x_start)*p.Ly*p.Lz + y*p.Lz + z][1] = -1.0*in[i][(x-x_start)*p.Ly*p.Lz
							 + ((p.Ly-y)%p.Ly)*p.Lz + ((p.Lz-z)%p.Lz)][1];


	  } else if(x >= p.Lx/2) {
	    /* in this case we conjugate the stuff in the lower
	       slices, using swap_row to work out where it was put.
	    */

	    in[i][(x-x_start)*p.Ly*p.Lz + y*p.Lz + z][0] = 1.0*swap_in[i][swap_row*p.Ly*p.Lz
								      + ((p.Ly-y)%p.Ly)*p.Lz + ((p.Lz-z)%p.Lz)][0];
	    in[i][(x-x_start)*p.Ly*p.Lz + y*p.Lz + z][1] = -1.0*swap_in[i][swap_row*p.Ly*p.Lz
									   + ((p.Ly-y)%p.Ly)*p.Lz + ((p.Lz-z)%p.Lz)][1];

	  } 
	  // otherwise do nothing


	}
      }
    }
  }



  //  MPI_Barrier(MPI_COMM_WORLD);
  //  printf0(p, "Done swapping\n");


  // Now planning  
  fftwf_plan plan = fftwf_mpi_plan_dft_3d(p.Lx, p.Ly, p.Lz,
				    in[0], out[0], MPI_COMM_WORLD,
				    FFTW_FORWARD, FFTW_ESTIMATE);







  //  MPI_Barrier(MPI_COMM_WORLD);
  //  printf0(p, "Done planning\n");

  /*
  for(x=x_start; x<(x_start+x_thickness); x++) {
    fprintf(stderr, "node %d, x %d:\t%10.3g %10.3g\t ->\n", p.rank, x, in[0][(x-x_start)*p.Ly*p.Lz][0], in[0][(x-x_start)*p.Ly*p.Lz][1]);
  } 
  */

  fftwf_mpi_execute_dft(plan, in[0], out[0]);
  fftwf_mpi_execute_dft(plan, in[1], out[1]);
  fftwf_mpi_execute_dft(plan, in[2], out[2]);

 
  MPI_Barrier(MPI_COMM_WORLD);
  //  printf0(p, "Done FFTing, example value %g %g\n", out[0][0][0], out[0][0][1]);

  /*
  for(x=x_start; x<(x_start+x_thickness); x++) {
    fprintf(stderr, "node %d, x %d:\t->\t%10.3g %10.3g\t\n", p.rank, x, out[0][(x-x_start)*p.Ly*p.Lz][0], out[0][(x-x_start)*p.Ly*p.Lz][1]);
  } 

  */

  float *slice = (float *)malloc(slab*p.Ly*p.Lz*sizeof(float));


  float maximag = 0.0;

  for(i=0; i<3; i++) {

    // prepare slice

    for(x=0; x<slab; x++) {
      for(y=0; y<p.Ly; y++) {
	for(z=0; z<p.Lz; z++) {
	  if(fabs(out[i][x*p.Ly*p.Lz + y*p.Lz + z][1]) > fabs(maximag)) {
	    maximag = out[i][x*p.Ly*p.Lz + y*p.Lz + z][1];
	  }
	  slice[x*p.Ly*p.Lz + y*p.Lz + z] = out[i][x*p.Ly*p.Lz + y*p.Lz + z][0];
	}
      }
    }

    //    fprintf(stderr,"maximag is %g\n", maximag);



    // segfault error below here

    for(x=0; x<p.Lx; x++) {

      // Check whether we are the source node for this slab                                                                                                                                                                                  
      if((x >= x_start) && ( x < (x_start + x_thickness))) {

	for(ry = 0; ry < ny; ry++) {

	  if((x/p.slicex == p.myposx) && (ry == p.myposy)) {


	    memcpy(&trim[(x-p.shiftx)*p.slicey*p.Lz],
		   &slice[(x-x_start)*p.Ly*p.Lz + ry*p.slicey*p.Lz],
		   p.slicey*p.Lz*sizeof(float));

	    continue;
	  }


	  MPI_Send(&slice[(x-x_start)*p.Ly*p.Lz + ry*p.slicey*p.Lz],
		   p.slicey*p.Lz,
		   MPI_FLOAT,
		   ry*nx + x/p.slicex,
		   x*ny + ry,
		   MPI_COMM_WORLD);

	}


      } else if((x >= p.shiftx) && (x < (p.shiftx + p.slicex))) {


	MPI_Recv(&trim[(x-p.shiftx)*p.slicey*p.Lz],
		 p.slicey*p.Lz,
		 MPI_FLOAT,
		 x/slab,
		 x*ny + p.myposy,
		 MPI_COMM_WORLD,
		 &status);

      }
    }


    // untrim
    
    /*
    for(x=1;x<=p.slicex;x++) {
      for(y=1;y<=p.slicey;y++) {
	for(z=0;z<p.Lz;z++) {	
	  for(i=0;i<3;i++) {

	    if(p.Lx > p.Lx) {
	      field[i][x][y][z] = out[i][2*(x+p.shiftx-1)*p.Ly*p.Lz + 2*(y+p.shifty-1)*p.Lz + 2*z][0];
	    } else {
	      field[i][x][y][z] = out[i][(x+p.shiftx-1)*p.Ly*p.Lz + (y+p.shifty-1)*p.Lz + z][0];
	    }

	  }
	}
      }

    */
    
    for(x=1; x<=p.slicex; x++) {
      for(y=1; y<=p.slicey; y++) {
        for(z=0; z<p.Lz; z++) {
          field[i][x][y][z] = trim[(x-1)*p.slicey*p.Lz + (y-1)*p.Lz + z];
        }
      }
    }
    


    halo_field(field[i], p);
  }


  free(slice);
  free(trim);

  fftwf_destroy_plan(plan);
  
  fftwf_free(in[0]);
  fftwf_free(swap_in[0]);
  fftwf_free(out[0]);
  fftwf_free(in[1]);
  fftwf_free(swap_in[1]);
  fftwf_free(out[1]);
  fftwf_free(in[2]);
  fftwf_free(swap_in[2]);
  fftwf_free(out[2]);

  free(in);
  free(out);



}






#endif // INITPS
