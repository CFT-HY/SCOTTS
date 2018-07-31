/** @file gw.c
 *
 * Fourier transform and GW power spectrum calculation.
 *
 * For the power specrum of a single field see fft.c
 * For the power spectrum of the fluid, see velps.c
 *
 * Projectors and expressions are principally taken from
 * "Gravitational wave background from reheating after hybrid inflation"
 * Garcia-Bellido, Figueroa and Sastre
 * PRD 77, 043517 (2008)
 */
#include "hydro.h"


#ifdef FFT


/** Projector \f$ P_{ij} \f$, where `T` corresponds to a combination
 * of \f$ i \f$ and \f$ j \f$.
 *
 * \f[
 * P_{ij}=\delta_{ij} -\hat{k}_{i}\hat{k}_{j}
 * \f]
 * 
 * 
 */
float proj(int T, float kx, float ky, float kz) {

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


/** The lambda projector object.
 *
 * `i`, `j`, `l`, `m` are the indices (NB: one-based)
 * and `kx`, `ky`, `kz` are momentum components.
 * 
 * The lambda projector object is given by
 * \f[ 
 * \Lambda_{ij,lm}(\mathbf{k}) = P_{im}(\mathbf{k}) P_{jl}(\mathbf{k}) - 
 *                             \dfrac{1}{2} P_{im}(\mathbf{k})P_{jl}(\mathbf{k})
 * \f]
 */
float lambda(int i, int j, int l, int m, float kx, float ky, float kz) {
  
  float total = 0.0;

  int cpt1, cpt2;



  if(i > j){
    cpt1 = i*10 + j;
  }
  else{
    cpt1 = j*10 + i;
  }
  if(l > m){
    cpt2 = l*10 + m;
  }
  else{
    cpt2 = m*10 + l;
  }

  total -= 0.5*proj(cpt1, kx, ky, kz)*proj(cpt2, kx, ky, kz);


  if(i > l){
    cpt1 = i*10 + l;
  }
  else{
    cpt1 = l*10 + i;
  }

  if(j > m){
    cpt2 = j*10 + m;
  }
  else{
    cpt2 = m*10 + j;
  }

  total += proj(cpt1, kx, ky, kz)*proj(cpt2, kx, ky, kz);

  return total;

}


/** Helper function to map indices `i` and `j` to the unique degrees of
 * encoded in our tensor object.
 *
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



/** Used to obtain \f$ |\dot{h}_{ij}|^2 \f$ in momentum space from `udot`.
 *
 *
 * Project from the six unique degrees of freedom in `udot`, (obtained
 * here in momentum space from a distributed shared memory FFT) to
 * obtain \f$ |\dot{h}_{ij}|^2 \f$ in momentum space, cf eq. 30 in the
 * paper cited at the start.
 *
 * `slab` is the thickness in the x-direction of the local domain
 *
 * `x_start` is the physical start of this region
 *
 * `y` and `z` have the full spatial extent `p.Ly` and `p.Lz`
 */ 
void gwproject(hydro_params p, int x_start, int slab,
	       float *product, fftwf_complex **udot) {

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





/** Calculates the gravitational wave power spectrum from `udotij`.
 *
 * Take the `udotij`'s, which are the un-projected
 * linearised metric perturbations, FFT them, project them in momentum
 * space and then calculate a power spectrum for them.
 *
 * The projector calculations (above) are slightly fiddly, so if you
 * decide to calculate the power spectrum several times during the
 * simulation be prepared to pre-calculate the above.
 *
 * Question: how is the output normalised, what exactly is this?
 */
float fft_tensor(hydro_fields f, hydro_params p, int step,
		float energydensity) {

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

  float *trim = (float *)malloc(p.slicex*p.slicey*p.Lz*sizeof(float));

  printf0(p, "Starting GW FFT calculation.\n");

  float start = clock();

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
	  //	  fprintf(stderr, "therefore node %d is responsible for x=%d\n", i, x);
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
  // slab is the thickness the 'used' nodes wil FFT
  // x_thickness is the thickness the current node will FFT
  slab = (int) x_thickness;

  if(((int)x_thickness) == 0) {
    fprintf(stderr, "Rank %d nothing to do for FFT: dx=0, x_start=%d!\n",
	    p.rank, (int)x_start);
    slab = 1;
  }
  /* else if(((int)x_thickness) < p.Lx/p.size) {
    fprintf(stderr,
	    "Rank %d giving up in FFT: FFTW told me to use a silly layout!\n",
	    p.rank);
    die(-42);
  }
  */


  fftwf_complex *in = fftwf_alloc_complex(alloc_local);
  fftwf_complex *out = fftwf_alloc_complex(alloc_local);

  fftwf_complex **outcpts = (fftwf_complex **)malloc(6*sizeof(fftwf_complex *));
  
  for(i=0;i<TENSOR_CPTS;i++) {
    outcpts[i] = fftwf_alloc_complex(alloc_local);
  }

  float *slice = (float *)malloc(x_thickness*p.Ly*p.Lz*sizeof(float));


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
	  check += f.udotij[i][x][y][z];
	  trim[(x-1)*p.slicey*p.Lz + (y-1)*p.Lz + z] = f.udotij[i][x][y][z];
	}
      }
    }

    //    fprintf(stderr,"Check: uij %lf\n", reduce_sum(check, p));
    
    for(x=0; x<p.Lx; x++) {
      
      // Check whether we are the target node for this slab
      if(map[x] == p.rank) {
	
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

	  //	  fprintf(stderr, "%d: got pencil with x=%d\n", p.rank, x);

	}
	
	
      } else if((x >= p.shiftx) && (x < (p.shiftx + p.slicex))) {

	/*       
	fprintf(stderr, "%d: my pencil, x=%d to dest %d/%d\n",
		p.rank, x, x, slab);
	*/	

	MPI_Send(&trim[(x-p.shiftx)*p.slicey*p.Lz],
		 p.slicey*p.Lz,
		 MPI_FLOAT,
		 map[x],
		 x*ny + p.myposy,
		 MPI_COMM_WORLD);
      }
      
    }

    printf0(p, "GW: Done slicing for cpt %d\n", i);

    float total = 0.0;

    for(x=0;x<x_thickness;x++) {
      for(y=0;y<p.Ly;y++) {
	for(z=0;z<p.Lz;z++) {
	  in[x*p.Ly*p.Lz + y*p.Lz + z][0] = slice[x*p.Ly*p.Lz + y*p.Lz + z];
	  in[x*p.Ly*p.Lz + y*p.Lz + z][1] = 0.0;
	  total += in[x*p.Ly*p.Lz + y*p.Lz + z][0];
	}
      }
    }


    printf0(p, "Executing GW FFT, cpt %d\n", i);
    /*
     * NB: fftw calculates the unnormalised FFT,
     * namely:
     * out[k] = \sum[j] in[j] exp(-2*pi*(j.k)*i/L)
     *
     * our fft_norm makes this sum turn into the integral
     * when the continuum limit is taken.
     */
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


    /*
    if(!i) {
      fprintf(stderr,"total outcpts %d is %6.10lf\n", p.rank,
	      reduce_sum(total,p));
    }
    */

    // now to deal with out
    //    fprintf(stderr,"Done Tensor FFT %d/%d\n", i, TENSOR_CPTS);

  }


  printf0(p, "GW projection...\n");
  // At this point, we should have a normed FFT

  // The brains of the operation:
  // Turn the FFT'd tensor into projected power spectrum
  gwproject(p, x_start, x_thickness, slice, outcpts);


  printf0(p, "Calculating total energy\n");

  float rhogw = 0.0;

  for(x=0; x<x_thickness; x++) {
    for(y=0; y<p.Ly; y++) {
      for(z=0; z<p.Lz; z++) {
	rhogw += slice[x*p.Ly*p.Lz + y*p.Lz + z];
	
      }
    }
  }


  // we take G=1.0...

  float gwen = reduce_sum(rhogw, p);
  //    /(1.0*((float)(p.Lx*p.Ly*p.Lz*p.dx*p.dx*p.dx)));

  // To get energy density expect to have to divide by V again
  // as it seems that keeping p.Lx*p.Ly*p.Lz*p.dx*p.dx*p.dx constant
  // then gwen is constant.
  printf0(p,
	  "Unnormalised GW energy [density] claimed %6.10lf\n", gwen);

  

  // Bin the power spectrum on the fly

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
  //  float comovingk, thisf, thisdiff, thisomega;

  // not sure (includes h^2)
  //  float omegarad = 3.5e-5;

  //  float doffrac = pow(1.0/100.0,1.0/3.0);

  // spokesman does the final analysis
  if(!p.rank) {

  char fftdest[200];

  sprintf(fftdest,"power-step%d.txt", step);

  FILE *fp = fopen(fftdest,"w");

    for(i=0;i<nbins;i++) {

      // so divide by dx
      //      comovingk = thisk/(p.a*p.dx);
      //      comovingk = thisk/p.dx;

      // formula taken as a given, looks correct
      //      thisf = comovingk*6.0e10; // /(p.H)
      //      thisf = comovingk;

      // d(rhogw)/d(logk)
      // (p.H^(-3) is comoving volume
      // and 8*4*pi is the solid angle normalisation)
      //      thisdiff = (comovingk*comovingk*comovingk/(32.0*M_PI))
      //	*bins[i]/(p.Lx*p.Ly*p.Lz*p.dx*p.dx*p.dx);
      //      thisomega = comovingk*comovingk*comovingk*bins[i]
      //	/(p.Lx*p.Ly*p.Lz*p.dx*p.dx*p.dx);

      // omegaGW(k), as above but corrections for degrees of freedom,
      // energy density and radiation density
      //      thisomega = omegarad*doffrac*(1.0/energydensity)*thisdiff;

      fprintf(fp, "%lf %g %d\n",
	      thisk/(p.dx*p.a), ((float)(i+1))*bins[i], counts[i]);

      thisk = thisk + dk;
    }

  fclose(fp);

  }

  free(bins);
  free(counts);


  // Tidy up

  free(slice);
  free(trim);

  fftwf_destroy_plan(plan);
  
  fftwf_free(in);
  fftwf_free(out);

  for(i=0;i<TENSOR_CPTS;i++)
    fftwf_free(outcpts[i]);

  free(outcpts);

  fftwf_mpi_cleanup();

  float end = clock();

  printf0(p, "GW FFT calculation took %lf\n",
	  ((float) (end - start)) / CLOCKS_PER_SEC);

  return gwen;
}




#endif // FFT
