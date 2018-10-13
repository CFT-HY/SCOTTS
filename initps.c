/* fft.c
 *
 * Fourier transform (and power spectrum) of a field
 */
#include "hydro.h"

#if defined(FFT) && ! defined(SCALAR)

// http://people.math.sc.edu/kellerlv/Quadratic_Interpolation.pdf 
// Lagrange Interp formula
float quadratic(float x1, float x2, float x3,
		float y1, float y2, float y3, float x) {

  return y1*((x - x2)*(x - x3))/((x1 - x2)*(x1 - x3)) 
	     + y2*((x - x1)*(x - x3))/((x2 - x1)*(x2 - x3))
	     + y3*((x - x1)*(x - x2))/((x3 - x1)*(x3 - x2));

  

}

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
          field[i][x][y][z] = p.initnorm*field[i][x][y][z]/momtot;
        }
      }
    }
    halo_field(field[i], p);
  }
}


void VtoZ(hydro_fields f, hydro_params p) {

  int i, x, y, z;

  // This assumes V has been initialised (but not necessarily haloed)

  // First, then, halo V
  halo_field(f.V[0], p);
  halo_field(f.V[1], p);
  halo_field(f.V[2], p);
 
  // Work out U from V (approximately)
  // (this can be made more precise by adopting what is in evolve.c)
  for(x=1; x<=p.slicex; x++) {
    for(y=1; y<=p.slicey; y++) {
      for(z=0; z<p.Lz; z++) {
	  f.W[x][y][z] = 1.0/sqrt(1 - (f.V[0][x][y][z]*f.V[0][x][y][z] +
				       f.V[1][x][y][z]*f.V[1][x][y][z] +
				       f.V[2][x][y][z]*f.V[2][x][y][z]));
	  
	  f.U[0][x][y][z] = f.W[x][y][z]*f.V[0][x][y][z];
	  f.U[1][x][y][z] = f.W[x][y][z]*f.V[1][x][y][z];
	  f.U[2][x][y][z] = f.W[x][y][z]*f.V[2][x][y][z];
      }
    }
  }

  // Halo W and U
  halo_field(f.W, p);
  halo_field(f.U[0], p);
  halo_field(f.U[1], p);
  halo_field(f.U[2], p);

  float eps = 1.0;
  float kappa = 1.0;


  // Work out Z and E from U and W (approximately)
  // (this can be made more precise by adopting what is in evolve.c)
  for(x=1; x<=p.slicex; x++) {
    for(y=1; y<=p.slicey; y++) {
      for(z=0; z<p.Lz; z++) {
	  
	f.E[x][y][z] = f.W[x][y][z]*eps;

	/*
	f.kappa[x][y][z] = kappa;

	f.p[x][y][z] = (f.kappa[x][y][z] - 1.0)*f.E[x][y][z]/f.W[x][y][z];
	*/

	f.p[x][y][z] = 1.0;
	f.kappa[x][y][z] = 1.0 + f.W[x][y][z]*f.p[x][y][z]/f.E[x][y][z];
      }
    }
  }
  

  // halo E and kappa
  halo_field(f.E, p);
  halo_field(f.kappa, p);




  float sigmabar;

  for(x=1; x<=p.slicex; x++) {
    for(y=1; y<=p.slicey; y++) {
      for(z=0; z<p.Lz; z++) {

	sigmabar = (f.kappa[x][y][z]
                    *f.E[x][y][z]
                    + f.kappa[x-1][y][z]
                    *f.E[x-1][y][z]
                    + f.kappa[x-1][y-1][z]
                    *f.E[x-1][y-1][z]
                    + f.kappa[x-1][y-1][((z+p.Lz-1)%p.Lz)]
                    *f.E[x-1][y-1][((z+p.Lz-1)%p.Lz)]
                    + f.kappa[x][y-1][((z+p.Lz-1)%p.Lz)]
                    *f.E[x][y-1][((z+p.Lz-1)%p.Lz)]
                    + f.kappa[x-1][y][((z+p.Lz-1)%p.Lz)]
                    *f.E[x-1][y][((z+p.Lz-1)%p.Lz)]
                    + f.kappa[x][y-1][z]
                    *f.E[x][y-1][z]
                    + f.kappa[x][y][((z-1+p.Lz)%p.Lz)]
                    *f.E[x][y][((z-1+p.Lz)%p.Lz)]
		    )/8.0;

	f.Z[0][x][y][z] = f.U[0][x][y][z]*sigmabar;
	f.Z[1][x][y][z] = f.U[1][x][y][z]*sigmabar;
	f.Z[2][x][y][z] = f.U[2][x][y][z]*sigmabar;
      }
    }
  }
  
  halo_field(f.Z[0], p);
  halo_field(f.Z[1], p);
  halo_field(f.Z[2], p);



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
	      //  Rot?
	      // v_i^\perp = P_{ij} v_j
	      // where P_{ij} = \delta_{ij} - \hat{k}_i \hat{k}_j
	      // and $\hat{k}$ is a unit vector in the $k$ direction.
	      /*
	      in_proj_re[i-1] += vel_proj(i*10 + j, kx, ky, kz)*in[j-1][x*p.Ly*p.Lz + y*p.Lz + z][0];
	      in_proj_im[i-1] += vel_proj(i*10 + j, kx, ky, kz)*in[j-1][x*p.Ly*p.Lz + y*p.Lz + z][1];
	      */


	      // #else
	      // #warning DIVPS enabled
	      // Div?
	      // v_i^\parallel = (\delta_{ij} - P_{ij}) v_j
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


/** Do a linear interpolation to compute the initial power spectrum.
 *
 * A key assumption is that the k-modes are _evenly distributed_ throughout
 * each bin, which isn't true! And increasing the order of the interpolation
 * will not really deal with this systematic error (I have tried quadratic).
 *
 * If this causes undue concern, there are two options as to how to fix it:
 * 1) When saving power spectrums, store the mean $k$ in each bin
 *    as well as the central $k$ for each bin, and use the mean $k$ (at least)
 *    when using power spectra as inputs.
 * 2) Recompute the mean $k$ here when building out the power spectrum again.
 *    That would make it harder to change lattice volumes between runs, as
 *    the relevant mode counts per bin for making the interpolation
 *    of the spectral density are for the original simulation
 *    that generated the input PS.
 */
void spectrum_interp(float ksq, hydro_params p, fftwf_complex *res,
		     float *k_bins, float *pow_bins, int n_bins) {

  float phase;

  float L = p.Lx;


  float amp, amp_last, amp_this, amp_next, grad;
  int i;


  amp_last = 0.0;
  amp_this = 0.0;
  amp_next = 0.0;
  amp = 0.0;

  // bin coordinates are midpoints
  float dk = 2.0*k_bins[0];

  float kmode = sqrt(fabs(ksq));
  int whichbin = (int)floor(kmode/dk);

  /*
  if(whichbin < (n_bins-1)) {
    amp = sqrt(pow_bins[whichbin]);
  }
  */


  if(whichbin > 0 && whichbin < (n_bins-1)) {
    // amp_last = sqrt(pow_bins[whichbin-1]);
    amp_this = sqrt(pow_bins[whichbin]);
    amp_next = sqrt(pow_bins[whichbin+1]);
    grad = (amp_next - amp_this)/dk;
    amp = amp_this + (kmode - ((float)(whichbin))*dk)*grad;

    // quadratic interpolation isn't any better!
    // amp = quadratic(((float)(whichbin-1))*dk, ((float)(whichbin))*dk,
    //		    ((float)(whichbin+1))*dk, amp_last, amp_this, amp_next,
    //		    kmode);


    /*

    if(isnan(amp)) {
      fprintf(stderr, "Our AMP IS A NAN!\n");
      fprintf(stderr, "kmode: %g, this: %g, bin %d, next: %g, grad: %g, amp: %g\n",
	      kmode, amp_this, whichbin, amp_next, grad, amp);
    }
    */

  }

  
  /*
  for(i = 0; i < n_bins; i++) {
    if((k_bins[i] - dk/2.0) <= sqrt(fabs(ksq)) && (k_bins[i] + dk/2.0) > sqrt(fabs(ksq))) {
      amp = sqrt(pow_bins[i]);
      break;
    }
  }
  */


  (*res)[0] = get_normal(0.0, amp);
  (*res)[1] = (*res)[0];

  phase = drand48();

  // Random phase
  (*res)[0] = (*res)[0]*cos(2.0*M_PI*phase);
  (*res)[1] = (*res)[1]*sin(2.0*M_PI*phase);


}




/* spectrum
 * 
 */
void spectrum(float ksq, hydro_params p, fftwf_complex *res) 
{ 

  float phase;

  float L = p.Lx;


  if(fabs(ksq) < 0.000001 || fabs(ksq) > p.initcutoff) {
    (*res)[0] = 0.0;
    (*res)[1] = 0.0;
  } else {


    //  (*res)[0] = p.initnorm*exp(-0.25*sqrt(ksq)*p.dx*((float)L));    
    //    (*res)[0] = get_normal(0.0, p.initnorm*exp(-0.25*sqrt(ksq)*p.dx*((float)L)));
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

  printf0(p, "Loading initial power spectrum from %s\n", p.initpsfile);


  ptrdiff_t x_thickness, x_start, alloc_local;


  ptrdiff_t n0 = p.Lx;
  ptrdiff_t n1 = p.Ly;
  ptrdiff_t n2 = p.Lz;

  //  int slab;

  MPI_Status status;

  float kmode, modsq;

  int x, y, z;
  int i, j;

  float *trim = (float *)malloc(p.slicex*p.slicey*p.Lz*sizeof(float));


  fftwf_mpi_init();

  alloc_local = fftwf_mpi_local_size_3d(n0, n1, n2,
                                       MPI_COMM_WORLD,
                                       &x_thickness, &x_start);


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

  printf0(p, "broadcasting our results\n");
  MPI_Bcast(map, p.Lx, MPI_INTEGER, 0, MPI_COMM_WORLD);
  printf0(p, "now to layout\n");



  fftwf_complex **in = (fftwf_complex **)malloc(3*sizeof(fftwf_complex *));
  in[0] = fftwf_alloc_complex(alloc_local);
  in[1] = fftwf_alloc_complex(alloc_local);
  in[2] = fftwf_alloc_complex(alloc_local);

  fftwf_complex **swap_in = (fftwf_complex **)malloc(3*sizeof(fftwf_complex *));
  swap_in[0] = fftwf_alloc_complex(x_thickness*p.Ly*p.Lz);
  swap_in[1] = fftwf_alloc_complex(x_thickness*p.Ly*p.Lz);
  swap_in[2] = fftwf_alloc_complex(x_thickness*p.Ly*p.Lz);

  fftwf_complex **out = (fftwf_complex **)malloc(3*sizeof(fftwf_complex *));
  out[0] = fftwf_alloc_complex(alloc_local);
  out[1] = fftwf_alloc_complex(alloc_local);
  out[2] = fftwf_alloc_complex(alloc_local);



  int nx = p.Lx/p.slicex;
  int ny = p.Ly/p.slicey;

  int ry;


  /*
   * Set up modes with desired power spectrum...
   */

  float ksq;

  //  int i;

  int true_x, true_y, true_z;

  //  int n_bins = p.Lx;
  float *k_bins = (float *)malloc(p.initpsbins*sizeof(float));
  float *pow_bins = (float *)malloc(p.initpsbins*sizeof(float));
  int in_bin;
  int items_read;
  float fudge;

  FILE *fp = fopen(p.initpsfile, "r");

  for(i = 0; i < p.initpsbins; i++) {

    if(feof(fp)) {
      printf0(p, "Fewer bins in file (%d) than initpsbins suggested (%d)\n",
	      i, p.initpsbins);
      die(100);
    }

    items_read = fscanf(fp, "%f%g%d", &k_bins[i], &pow_bins[i], &in_bin);

    if(items_read != 3) {
      printf0(p,
	      "Initpsfile %s: row %d: did not read a full line (3 items), "
	      "just %d\n",
	      p.initpsfile, i, items_read);
      die(100);
    }

    if(in_bin == 0) {
      pow_bins[i] = 0.0;
    } else {
      pow_bins[i] /= ((float)(((long int)in_bin)*((long int)(i+1))));
    }

    if(isnan(pow_bins[i])) {
      printf0(p, "Bin %d is a NaN\n", i);
    }
  }
  fclose(fp);

  fudge = (2.0*M_PI/((float)p.Lx))/((k_bins[1]-k_bins[0])*p.dx);
  if(fabs(fudge - 1.0) > 1e-6) {
    printf0(p, "Applying fudge factor %g^3 to power spectrum: "
	    "seems like dk or volume or dx is different\n",
	    fudge);

    for(i = 0; i < p.initpsbins; i++) {
      pow_bins[i] *= fudge*fudge*fudge;
    }
  }

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

	/*
	ksq = (2.0 - 2.0*cos(((float)true_x)*2.0*M_PI/(((float)p.Lx))))/(p.dx*p.dx);
	ksq += (2.0 - 2.0*cos(((float)true_y)*2.0*M_PI/(((float)p.Ly))))/(p.dx*p.dx);
	ksq += (2.0 - 2.0*cos(((float)true_z)*2.0*M_PI/(((float)p.Lz))))/(p.dx*p.dx);
	*/

	ksq = (
	       ((float)(true_x*true_x))/((float)(p.Lx*p.Lx))
	       + ((float)(true_y*true_y))/((float)(p.Ly*p.Ly))
	       + ((float)(true_z*true_z))/((float)(p.Lz*p.Lz))
	       )*(2.0*M_PI)*(2.0*M_PI)/(p.dx*p.dx);

	for(i=0; i<3; i++) {
	  // spectrum(ksq, p, &(in[i][(x-x_start)*p.Ly*p.Lz + y*p.Lz + z]));
	  spectrum_interp(ksq, p,
			  &(in[i][(x-x_start)*p.Ly*p.Lz + y*p.Lz + z]),
			  k_bins, pow_bins, p.initpsbins);
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

    printf0(p, "Direction %d being handled\n", i);

    // At this point each node knows what it ought to have, but we need to 'hermitianise' it
    // Need to do swapping






    for(j=0; j<=p.Lx; j++) {

      if(j>= x_start && j < x_start+x_thickness) {




	if(map[((p.Lx-j) % p.Lx)] == p.rank) {
	} else {
	  MPI_Send(&(((float *)in[i])[((j-x_start)%p.Lx)*p.Ly*p.Lz*2]),
		   2*p.Ly*p.Lz, MPI_FLOAT,
		   map[((p.Lx-j) % p.Lx)], j, MPI_COMM_WORLD);
	}
      }


      if(j>=min_needed && j<=max_needed) {

	if(map[(j % p.Lx)] == p.rank) {


	  memcpy(&(((fftwf_complex *)swap_in[i])[(j-min_needed)*p.Ly*p.Lz]),
		 &(((fftwf_complex *)in[i])[((j-x_start)%p.Lx)*p.Ly*p.Lz]),
		 p.Ly*p.Lz*sizeof(fftwf_complex));


	} else {
	  MPI_Recv(&(((float *)swap_in[i])[(j-min_needed)*p.Ly*p.Lz*2]),
		   2*p.Ly*p.Lz, MPI_FLOAT,
		   map[(j % p.Lx)], j, MPI_COMM_WORLD, &status);

	}
      }



    }
      


    // at this point swap_in has the stuff needed to do the conjugate properly


    printf0(p, "Starting swap\n");


    for(x=x_start; x < x_start+x_thickness; x++) {


      // This is the row in the swap_in where one can find (p.Lx-x)
      //      int swap_row = (((p.Lx-x) % x_thickness) + x_thickness - 1) % x_thickness;
      int swap_row = (((p.Lx-x) - min_needed));

      for(y=0;y<p.Ly;y++) {
	for(z=0;z<p.Lz;z++) {

	  // Conjugate the system. Treat special cases first.

	  if((x == 0 || x == p.Lx/2)
	     && (y == 0 || y == p.Ly/2)
	     && (z == 0 || z == p.Lz/2)) {
	    

	    /* These sites need to be pure real.
	     */	    
	    in[i][(x-x_start)*p.Ly*p.Lz + y*p.Lz + z][0] \
	      = sqrt(swap_in[i][swap_row*p.Ly*p.Lz
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
	     * if we used swap_in we'd not end up with something that is
	     * hermitian because otherwise:
	     * [<>          []]
	     * initially would be swapped (and conjugated)
	     * to become:
	     * [[]          <>]
	     */
	    in[i][(x-x_start)*p.Ly*p.Lz + y*p.Lz + z][0] = \
	      in[i][(x-x_start)*p.Ly*p.Lz
		      + ((p.Ly-y)%p.Ly)*p.Lz + ((p.Lz-z)%p.Lz)][0];
	    in[i][(x-x_start)*p.Ly*p.Lz + y*p.Lz + z][1] = \
	      -1.0*in[i][(x-x_start)*p.Ly*p.Lz
			   + ((p.Ly-y)%p.Ly)*p.Lz + ((p.Lz-z)%p.Lz)][1];


	  } else if(x >= p.Lx/2) {
	    /* in this case we conjugate the stuff in the lower
	       slices, using swap_row to work out where it was put.
	    */

	    in[i][(x-x_start)*p.Ly*p.Lz + y*p.Lz + z][0] = \
	      1.0*swap_in[i][swap_row*p.Ly*p.Lz
			       + ((p.Ly-y)%p.Ly)*p.Lz + ((p.Lz-z)%p.Lz)][0];
	    in[i][(x-x_start)*p.Ly*p.Lz + y*p.Lz + z][1] = \
	      -1.0*swap_in[i][swap_row*p.Ly*p.Lz
				+ ((p.Ly-y)%p.Ly)*p.Lz + ((p.Lz-z)%p.Lz)][1];

	  } 
	  // otherwise do nothing


	}
      }
    }

    MPI_Barrier(MPI_COMM_WORLD);

  }





  MPI_Barrier(MPI_COMM_WORLD);
  printf0(p, "Done swapping\n");


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

  float *slice = (float *)malloc(x_thickness*p.Ly*p.Lz*sizeof(float));


  float maximag = 0.0;

  for(i=0; i<3; i++) {

    // prepare slice

    for(x=0; x<x_thickness; x++) {
      for(y=0; y<p.Ly; y++) {
	for(z=0; z<p.Lz; z++) {
	  if(fabs(out[i][x*p.Ly*p.Lz + y*p.Lz + z][1]) > fabs(maximag)) {
	    maximag = out[i][x*p.Ly*p.Lz + y*p.Lz + z][1];
	  }
	  slice[x*p.Ly*p.Lz + y*p.Lz + z] = out[i][x*p.Ly*p.Lz + y*p.Lz + z][0];
	}
      }
    }

    fprintf(stderr,"maximag is %g\n", maximag);
    MPI_Barrier(MPI_COMM_WORLD);
    

    // segfault error below here

    for(x=0; x<p.Lx; x++) {

      // Check whether we are the source node for this slab                                                                                                                                                                                  
      if(map[x] == p.rank) {
	
        for(ry = 0; ry < ny; ry++) {
	  
          if((x/p.slicex == p.myposx) && (ry == p.myposy)) {
            //      fprintf(stderr, "rank %d: doing memcpy and local thickness is %d\n",
            //              p.rank, x_thickness);

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


	//        fprintf(stderr, "rank %d: receiving into location %p\n", p.rank,
	//		&trim[(x-p.shiftx)*p.slicey*p.Lz]);

	MPI_Recv(&trim[(x-p.shiftx)*p.slicey*p.Lz],
                 p.slicey*p.Lz,
                 MPI_FLOAT,
                 map[x],
                 x*ny + p.myposy,
                 MPI_COMM_WORLD,
                 &status);

      }
    }


    // untrim
    
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


#endif // FFT && !SCALAR


