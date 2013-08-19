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

    // (*res)[0] = p.initps*exp(-0.25*sqrt(ksq)*p.dx*((double)L));
    (*res)[0] = get_normal(0.0, p.initps*exp(-0.25*sqrt(ksq)*p.dx*((double)L)));

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

  double start = clock();

  //  fftw_mpi_init();

  double checken = 0.0;




  int sizex = 512;
  int sizey = 512;
  int sizez = 512;
  double tdx = 0.5;

  fftw_complex **in = (fftw_complex **)malloc(3*sizeof(fftw_complex *));
  fftw_complex **out = (fftw_complex **)malloc(3*sizeof(fftw_complex *));

  in[0] = fftw_alloc_complex(sizex*sizey*sizez);
  out[0] = fftw_alloc_complex(sizex*sizey*sizez);

  in[1] = fftw_alloc_complex(sizex*sizey*sizez);
  out[1] = fftw_alloc_complex(sizex*sizey*sizez);

  in[2] = fftw_alloc_complex(sizex*sizey*sizez);
  out[2] = fftw_alloc_complex(sizex*sizey*sizez);



  int nx = p.Lx/p.slicex;
  int ny = p.Ly/p.slicey;

  int ry;


  /*
   * Set up modes with desired power spectrum...
   */
  int true_x;

  double ksq;

  //  int i;


  for(x=0;x<sizex;x++) {
    for(y=0;y<sizey;y++) {
      for(z=0;z<sizez;z++) {

	ksq = (2.0 - 2.0*cos(((double)x)*2.0*M_PI/(((double)sizex))))/(tdx*tdx);
	ksq += (2.0 - 2.0*cos(((double)y)*2.0*M_PI/(((double)sizey))))/(tdx*tdx);
	ksq += (2.0 - 2.0*cos(((double)z)*2.0*M_PI/(((double)sizez))))/(tdx*tdx);

	for(i=0; i<3; i++) {
	  spectrum(ksq, p, &(in[i][x*sizey*sizez + y*sizez + z]));
	}

      }
    }
  }

  // moduli will not work in parallel, also not swapped
  for(x=0;x<sizex;x++) {
    for(y=0;y<sizey;y++) {
      for(z=0;z<sizez;z++) {
	for(i=0; i<3; i++) {
	  in[i][x*sizey*sizez + y*sizez + z][0] = in[i][((sizex-x)%sizex)*sizey*sizez
					       + ((sizey-y)%sizey)*sizez + ((sizez-z)%sizez)][0];
	  in[i][x*sizey*sizez + y*sizez + z][1] = -1.0*in[i][((sizex-x)%sizex)*sizey*sizez
						    + ((sizey-y)%sizey)*sizez + ((sizez-z)%sizez)][1];

	  if((x == 0 || x == sizex/2) && (y == 0 || y == sizey/2) && (z == 0 || z == sizez/2)) {
	    in[i][x*sizey*sizez + y*sizez + z][0] = sqrt(in[i][((sizex-x)%sizex)*sizey*sizez
						      + ((sizey-y)%sizey)*sizez + ((sizez-z)%sizez)][0]
						   *in[i][((sizex-x)%sizex)*sizey*sizez
						       + ((sizey-y)%sizey)*sizez + ((sizez-z)%sizez)][0]
						   + in[i][((sizex-x)%sizex)*sizey*sizez
							+ ((sizey-y)%sizey)*sizez + ((sizez-z)%sizez)][1]
						   *in[i][((sizex-x)%sizex)*sizey*sizez
						       + ((sizey-y)%sizey)*sizez + ((sizez-z)%sizez)][1]);
	    in[i][x*sizey*sizez + y*sizez + z][1] = 0.0;
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
    for(x=0;x<sizex;x++) {
      for(y=0;y<sizey;y++) {
	for(z=0;z<sizez;z++) {



	  
	  kx = sqrt((2.0 - 2.0*cos(((double)x)*2.0*M_PI/(((double)sizex))))/(tdx*tdx));
	  ky = sqrt((2.0 - 2.0*cos(((double)y)*2.0*M_PI/(((double)sizex))))/(tdx*tdx));
	  kz = sqrt((2.0 - 2.0*cos(((double)z)*2.0*M_PI/(((double)sizex))))/(tdx*tdx));
	  
	  in_proj_re[0] = 0.0;
	  in_proj_re[1] = 0.0;
	  in_proj_re[2] = 0.0;
	  
	  in_proj_im[0] = 0.0;
	  in_proj_im[1] = 0.0;
	  in_proj_im[2] = 0.0;

	  
	  for(i=1; i<=3; i++) {

	    double res_re = 0.0;
	    double res_im = 0.0;

	    for(j=1; j<=3; j++) {

#ifndef DIVPS	      
	      in_proj_re[i-1] += vel_proj(i*10 + j, kx, ky, kz)*in[j-1][x*sizey*sizez + y*sizez + z][0];
	      in_proj_im[i-1] += vel_proj(i*10 + j, kx, ky, kz)*in[j-1][x*sizey*sizez + y*sizez + z][1];
#else
#warning DIVPS enabled
	      if(i == j) {
		in_proj_re[i-1] += (1.0-vel_proj(i*10 + j, kx, ky, kz))*in[j-1][x*sizey*sizez + y*sizez + z][0];
		in_proj_im[i-1] += (1.0-vel_proj(i*10 + j, kx, ky, kz))*in[j-1][x*sizey*sizez + y*sizez + z][1];
	      } else {
		in_proj_re[i-1] += (-1.0*vel_proj(i*10 + j, kx, ky, kz))*in[j-1][x*sizey*sizez + y*sizez + z][0];
		in_proj_im[i-1] += (-1.0*vel_proj(i*10 + j, kx, ky, kz))*in[j-1][x*sizey*sizez + y*sizez + z][1];
	      }
#endif




	      if(i == j) {
		res_re += (1.0-vel_proj(i*10 + j, kx, ky, kz))
		  *in[j-1][x*sizey*sizez + y*sizez + z][0];
		res_im += (1.0-vel_proj(i*10 + j, kx, ky, kz))
		  *in[j-1][x*sizey*sizez + y*sizez + z][1];
	      } else {
		res_re += (-1.0*vel_proj(i*10 + j, kx, ky, kz))
		  *in[j-1][x*sizey*sizez + y*sizez + z][0];
		res_im += (-1.0*vel_proj(i*10 + j, kx, ky, kz))
		  *in[j-1][x*sizey*sizez + y*sizez + z][1];
	      }

	      
			
	    }
	      stuff += in_proj_re[i-1]*in_proj_re[i-1] + in_proj_im[i-1]*in_proj_im[i-1];
	      resid += res_re*res_re + res_im*res_im;



	  }


	  in[0][x*sizey*sizez + y*sizez + z][0] = in_proj_re[0];
	  in[1][x*sizey*sizez + y*sizez + z][0] = in_proj_re[1];
	  in[2][x*sizey*sizez + y*sizez + z][0] = in_proj_re[2];

	  in[0][x*sizey*sizez + y*sizez + z][1] = in_proj_im[0];
	  in[1][x*sizey*sizez + y*sizez + z][1] = in_proj_im[1];
	  in[2][x*sizey*sizez + y*sizez + z][1] = in_proj_im[2];

	}
      }
    }

    if(!p.rank)
      fprintf(stderr,"run %d, stuff is %g, resid is now %g\n", reruns, stuff, resid);
  }

  // Now planning  
  fftw_plan plan = fftw_plan_dft_3d(sizex, sizey, sizez,
				    in[0], out[0],
				    FFTW_FORWARD, FFTW_ESTIMATE);


  
  fftw_execute_dft(plan, in[0], out[0]);
  fftw_execute_dft(plan, in[1], out[1]);
  fftw_execute_dft(plan, in[2], out[2]);




  for(x=1;x<=p.slicex;x++) {
    for(y=1;y<=p.slicey;y++) {
      for(z=0;z<p.Lz;z++) {	
	for(i=0;i<3;i++) {

	  if(sizex > p.Lx) {
	    field[i][x][y][z] = out[i][2*(x+p.shiftx-1)*sizey*sizez + 2*(y+p.shifty-1)*sizez + 2*z][0];
	  } else {
	    field[i][x][y][z] = out[i][(x+p.shiftx-1)*p.Ly*p.Lz + (y+p.shifty-1)*p.Lz + z][0];
	  }

	}
      }
    }
  }

  halo_field(field[0], p);
  halo_field(field[1], p);
  halo_field(field[2], p);


  fftw_destroy_plan(plan);
  
  fftw_free(in[0]);
  fftw_free(out[0]);
  fftw_free(in[1]);
  fftw_free(out[1]);
  fftw_free(in[2]);
  fftw_free(out[2]);

  free(in);
  free(out);


}






#endif // INITPS
