/* gw.c
 *
 * Fourier transform and GW spectrum calculation.
 *
 * Projectors and expressions are principally taken from
 * "Gravitational wave background from reheating after hybrid inflation"
 * Garcia-Bellido, Figueroa and Sastre
 * PRD 77, 043517 (2008)
 */
#include "hydro.h"


#ifdef FFT


/* proj(int T, double kx, double ky, double kz);
 *
 * Projector P_{ij}, where ij is parameterised by T
 * Given above eq. 20 in the paper.
 */
double proj(int T, double kx, double ky, double kz) {

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
  case 22:
    return 1.0 - 1.0*ky*ky;
    break;
  case 32:
    return -1.0*ky*kz;
    break;
  case 33:
    return 1.0 - 1.0*kz*kz;
    break;

  default:
    fprintf(stderr,"Unknown projector element! Nonsense will ensue...\n");

  }

  return 0.0;
}


/* lambda(int i, int j, int l, int m, double kx, double ky, double kz)
 *
 * the lambda projector object
 * i, j, l, m are the indices (NB: one-based)
 * and kx, ky, kz are momentum components
 */
double lambda(int i, int j, int l, int m, double kx, double ky, double kz) {
  
  double total = 0.0;

  int cpt1, cpt2;



  if(i > j)
    cpt1 = i*10 + j;
  else
    cpt1 = j*10 + i;

  if(l > m)
    cpt2 = l*10 + m;
  else
    cpt2 = m*10 + l;

  total -= 0.5*proj(cpt1, kx, ky, kz)*proj(cpt2, kx, ky, kz);


  if(i > l)
    cpt1 = i*10 + l;
  else
    cpt1 = l*10 + i;

  if(j > m)
    cpt2 = j*10 + m;
  else
    cpt2 = m*10 + j;

  total += proj(cpt1, kx, ky, kz)*proj(cpt2, kx, ky, kz);

  return total;

}


/* indexof(int i, int j)
 *
 * For individual indices i and j associated with a single tensor,
 * map back to the six unique degrees of freedom encoded in our tensor object.
 */
int indexof(int i, int j)
{

  if(i == 1) {
    if(j == 1) {
      return CPT_11;
    } else if(j == 2) {
      return CPT_21;
    } else if(j == 3) {
      return CPT_31;
    }
  } else if(i == 2) {
    if(j == 1) {
      return CPT_21;
    } else if(j == 2) {
      return CPT_22;
    } else if(j == 3) {
      return CPT_32;
    }
  } else if(i == 3) {
    if(j == 1) {
      return CPT_31;
    } else if(j == 2) {
      return CPT_32;
    } else if(j == 3) {
      return CPT_33;
    }
  }

  fprintf(stderr,"Unrecognised indices: %d, %d\n", i, j);

  return 99;
}



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
void gwproject(hydro_params p, int x_start, int slab,
	       double *product, fftw_complex **udot) {

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
		  *(udot[indexof(i,j)][x*p.Ly*p.Lz + y*p.Lz + z][0]
		   *udot[indexof(l,m)][x*p.Ly*p.Lz + y*p.Lz + z][0]
		    + udot[indexof(i,j)][x*p.Ly*p.Lz + y*p.Lz + z][1]
		    *udot[indexof(l,m)][x*p.Ly*p.Lz + y*p.Lz + z][1]);
	      }
	    }
	  }
	}

      }
    }
  }
}





/* fft_tensor(hydro_fields f, hydro_params p)
 *
 * Take the udotij's, which are the un-projected linearised
 * metric perturbations, FFT them, project them in momentum space
 * and then calculate a power spectrum for them.
 *
 * The projector calculations (above) are slightly fiddly, so if
 * you decide to calculate the power spectrum several times during the
 * simulation be prepared to pre-calculate the above.
 */
void fft_tensor(hydro_fields f, hydro_params p, int step,
		double energydensity) {

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

  fftw_complex **outcpts = (fftw_complex **)malloc(6*sizeof(fftw_complex *));
  
  for(i=0;i<TENSOR_CPTS;i++) {
    outcpts[i] = fftw_alloc_complex(alloc_local);
  }

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
	  check += f.udotij[i][x][y][z];
	  trim[(x-1)*p.slicey*p.Lz + (y-1)*p.Lz + z] = f.udotij[i][x][y][z];
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


    /*
    if(!i) {
      fprintf(stderr,"total outcpts %d is %6.10lf\n", p.rank,
	      reduce_sum(total,p));
    }
    */

    // now to deal with out
    //    fprintf(stderr,"Done Tensor FFT %d/%d\n", i, TENSOR_CPTS);

  }

  /*
  if(!p.rank)
    fprintf(stderr,"slice 111: %.12lf\n", outcpts[0][p.Ly*p.Lz + p.Lz + 1][0]);
  */


  // At this point, we should have a normed FFT

  // The brains of the operation:
  // Turn the FFT'd tensor into projected power spectrum
  gwproject(p, x_start, slab, slice, outcpts);



  double rhogw = 0.0;

  for(x=0; x<slab; x++) {
    for(y=0; y<p.Ly; y++) {
      for(z=0; z<p.Lz; z++) {
	rhogw += slice[x*p.Ly*p.Lz + y*p.Lz + z];
	
      }
    }
  }


  // we take G=1.0...

  double gwen = reduce_sum(rhogw, p)
    /(1.0*((double)(p.a*p.a*p.a*p.Lx*p.Ly*p.Lz*p.dx*p.dx*p.dx)));

  // To get energy density expect to have to divide by V again
  // as it seems that keeping p.Lx*p.Ly*p.Lz*p.dx*p.dx*p.dx constant
  // then gwen is constant.
  printf0(p,
	  "Unnormalised GW energy [density] claimed %6.10lf\n", gwen);

  
  // **** now calculate and output the power spectrum ****



#ifdef DUMPFFT
  // Temporary until we get the binning sorted out
  char fftdest[200];

  sprintf(fftdest,"fft-%d.txt",p.rank);

  FILE *fp = fopen(fftdest,"w");
 

  for(x=0;x<slab;x++) {
    for(y=0;y<p.Ly;y++) {
      for(z=0;z<p.Lz;z++) {

	//	if((x+x_start)==0 && y==0 && z==0)
	//	  continue;

        if(((x+((int)x_start))>p.Lx/2) || (y> p.Ly/2) || (z>p.Lz/2))
          continue;

	kmode = sqrt(
		     ((double)((x+((int)x_start))
			       *(x+((int)x_start))))/((double)(p.Lx*p.Lx))
		     + ((double)(y*y))/((double)(p.Ly*p.Ly))
		     + ((double)(z*z))/((double)(p.Lz*p.Lz))
		     )*2.0*M_PI;

	/*
	modsq = out[x*p.Ly*p.Lz + y*p.Lz + z][0]
	  *out[x*p.Ly*p.Lz + y*p.Lz + z][0]
	  + out[x*p.Ly*p.Lz + y*p.Lz + z][1]
	  *out[x*p.Ly*p.Lz + y*p.Lz + z][1];
	  */

	fprintf(fp, "%lf %d %d %d %lf\n", kmode,
		x + ((int)x_start), y, z, slice[x*p.Ly*p.Lz + y*p.Lz + z]);


      }

      fflush(fp);

    }
  }

  close(fp);

#else // not DUMPFFT
  // calculate the power spectrum on the fly

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
  double comovingk, thisf, thisdiff, thisomega;

  // not sure (includes h^2)
  double omegarad = 3.5e-5;

  double doffrac = pow(1.0/100.0,1.0/3.0);

  // spokesman does the final analysis
  if(!p.rank) {

  char fftdest[200];

  sprintf(fftdest,"power-step%d.txt", step);

  FILE *fp = fopen(fftdest,"w");

    for(i=0;i<nbins;i++) {

      // thisk = sqrt((nx/Lx)^2 + (ny/Ly)^2 + (nz/Lz)^2) 2*pi
      // so divide by dx
      comovingk = thisk/(p.a*p.dx);

      // formula taken as a given, looks correct
      thisf = comovingk*6.0e10; // /(p.H)

      // d(rhogw)/d(logk)
      // (p.H^(-3) is comoving volume; 8*4*pi is the solid angle normalisation)
      thisdiff = (comovingk*comovingk*comovingk/(32.0*M_PI))
	*bins[i]/(p.Lx*p.Ly*p.Lz*p.dx*p.dx*p.dx);

      // omegaGW(k), as above but corrections for degrees of freedom,
      // energy density and radiation density
      thisomega = omegarad*doffrac*(1.0/energydensity)*thisdiff;

      fprintf(fp, "%lf %lf %d %g %g\n",
	      thisk, bins[i], counts[i], thisf, thisomega);

      thisk = thisk + dk;
    }

  fclose(fp);

  }



#endif


  // Tidy up

  free(slice);
  free(trim);

  fftw_destroy_plan(plan);
  
  fftw_free(in);
  fftw_free(out);

  for(i=0;i<TENSOR_CPTS;i++)
    fftw_free(outcpts[i]);


  fftw_mpi_cleanup();

  double end = clock();

  printf0(p, "FFT stuff took %lf\n",
	  ((double) (end - start)) / CLOCKS_PER_SEC);

}




#endif // FFT
