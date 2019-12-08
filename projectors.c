/** @file projectors.c
 *
 * File containing projectors for projecting vector and tensor fields
 * in fourier space.
 *
 * - Fourier transforms themselves are contained in `fft_func.c`
 * - The tensor power spectra binning is in tensorps.c,
 * - Vector field power spectra binning is in vectorps.c,
 * - Scalar field power spectra binning is in scalarps.c.
 * - Unequal time correlators are in uetc.c.
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
  case 12:
    return -1.0*kx*ky;
    break;
  case 13:
    return -1.0*kx*kz;
    break;
  case 23:
    return -1.0*ky*kz;
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

/** Used to project tensors into the transverse traceless gauge for
 * the power spectrum.  e.g populates `product` with 
 * \f$ p(\mathbf{k}) = \Lambda(\mathbf{k})_{ij,lm} T_{ij}(\mathbf{k})
 * T^*_{lm}(\mathbf{k}) \f$ where \f$ p \f$ is `product` and `T` is
 * `tensor`.
 *
 *
 * Project from the six unique degrees of freedom in `tensor`,
 * (obtained here in momentum space from a distributed shared memory
 * FFT) to obtain the modulus squared of the transverse traceless
 * projection of `tensor` in momentum space. Uses the projector
 * function `lambda` described above.
 *
 * `slab` is the thickness in the x-direction of the local domain
 *
 * `x_start` is the physical start of this region
 *
 * `y` and `z` have the full spatial extent `p.Ly` and `p.Lz`
 */ 
void tens_proj(hydro_params p, int x_start, int slab,
	       float *product, fftwf_complex **tensor) {

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
		  *(tensor[indexof(i,j)][x*p.Ly*p.Lz + y*p.Lz + z][0]
		   *tensor[indexof(l,m)][x*p.Ly*p.Lz + y*p.Lz + z][0]
		    + tensor[indexof(i,j)][x*p.Ly*p.Lz + y*p.Lz + z][1]
		    *tensor[indexof(l,m)][x*p.Ly*p.Lz + y*p.Lz + z][1]);

	      }
	    }
	  }
	}
      }
    }
  }
}

/** Like `tens_proj`, but for two separate tensors, and populates
 * instead two arrays with the real and complex part of the projected
 * tensor combination.
 *
 *
 * In more detail, this populates two float arrays: `product_re` \f$
 * p_{r} \f$ and `product_im` \f$ p_{im} \f$ with:
 *
 * \f[ p_{r}(\mathbf{k}) + i p_{im}(\mathbf{k}) = \Lambda_{ij,lm}
 * X_{ij}(\mathbf{k}) * Y^*_{lm}(\mathbf{k}) \f]
 *
 * and \f$X\f$ and \f$Y\f$ are `tensora` and `tensorb` respectively.
 */
void uetc_tens_proj(hydro_params p, int x_start, int slab,
		  float *product_re, float *product_im,
		  fftwf_complex **tensora, fftwf_complex **tensorb) {

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

	// Sum over indices on lambda_{ij,lm} and tensora_{ij}(k) tensorb_{lm}(k)
	for(i=1; i<=3; i++) {
	  for(j=1; j<=3; j++) {
	    for(l=1; l<=3; l++) {
	      for(m=1; m<=3; m++) {

		product_re[x*p.Ly*p.Lz + y*p.Lz + z] +=
		  lambda(i, j, l, m, kx, ky, kz)
		  *(tensora[indexof(i,j)][x*p.Ly*p.Lz + y*p.Lz + z][0]
		   *tensorb[indexof(l,m)][x*p.Ly*p.Lz + y*p.Lz + z][0]
		    + tensora[indexof(i,j)][x*p.Ly*p.Lz + y*p.Lz + z][1]
		    *tensorb[indexof(l,m)][x*p.Ly*p.Lz + y*p.Lz + z][1]);

		product_im[x*p.Ly*p.Lz + y*p.Lz + z] +=
		  lambda(i, j, l, m, kx, ky, kz)
		  *(-tensora[indexof(i,j)][x*p.Ly*p.Lz + y*p.Lz + z][0]
		   *tensorb[indexof(l,m)][x*p.Ly*p.Lz + y*p.Lz + z][1]
		    + tensora[indexof(i,j)][x*p.Ly*p.Lz + y*p.Lz + z][1]
		    *tensorb[indexof(l,m)][x*p.Ly*p.Lz + y*p.Lz + z][0]);
	      }
	    }
	  }
	}

      }
    }
  }
}

/** Split the momentum space vector `vk` into the transverse and
 * longitudinal ('div') components.
 *
 * Populates the arrays:
 *
 * `product_rot` 
 * \f[ p_{rot}(\mathbf{k}) 
 * = P_{ij}(\mathbf{k}P_{il}(\mathbf{k}) v_j(\mathbf{k}) v^*_l(\mathbf{k}) \f]
 *
 * `product_div` 
 * \f[ p_{div}(\mathbf{k}) 
 * = (\delta_{ij} - P_{ij}(\mathbf{k}))(\delta_{im} - P_{im}(\mathbf{k}))
 * v_i(\mathbf{k}) v^*_m(\mathbf{k}) \f]
 *
 * `product_tot` \f[ p_{tot}(\mathbf{k}) = v_i(\mathbf{k}) v^*_i(\mathbf{k}) \f]
 *
 *  
 */
void split_vector(hydro_params p, int x_start, int slab,
		     float *product, float *product_div, float *product_tot,
		     fftwf_complex **vk) {

  int i, j;
  int x, y, z;
  float kx, ky, kz;

  float res_r, res_i;
  float resid_r, resid_i;
  float tot_r, tot_i;

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


        kx = s_x*sqrt((2.0 - 2.0*cos(((float)(true_x))*2.0*M_PI/(((float)p.Lx))))/(p.dx*p.dx));
	ky = s_y*sqrt((2.0 - 2.0*cos(((float)(true_y))*2.0*M_PI/(((float)p.Ly))))/(p.dx*p.dx));
        kz = s_z*sqrt((2.0 - 2.0*cos(((float)(true_z))*2.0*M_PI/(((float)p.Lz))))/(p.dx*p.dx));


	/*
	kx = ((float)true_x)*2.0*M_PI/((float)p.Lx);
	ky = ((float)true_y)*2.0*M_PI/((float)p.Ly);
	kz = ((float)true_z)*2.0*M_PI/((float)p.Lz);
	*/

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
	    res_r += proj(i*10 + j, kx, ky, kz)
	      *vk[j-1][x*p.Ly*p.Lz + y*p.Lz + z][0];
	    res_i += proj(i*10 + j, kx, ky, kz)
	      *vk[j-1][x*p.Ly*p.Lz + y*p.Lz + z][1];

	    // And longitudinal...
	    if(i == j) {
	    resid_r += (1.0 - proj(i*10 + j, kx, ky, kz))
	      *vk[j-1][x*p.Ly*p.Lz + y*p.Lz + z][0];
	    resid_i += (1.0 - proj(i*10 + j, kx, ky, kz))
	      *vk[j-1][x*p.Ly*p.Lz + y*p.Lz + z][1];
	    } else {
	    resid_r += (-1.0*proj(i*10 + j, kx, ky, kz))
	      *vk[j-1][x*p.Ly*p.Lz + y*p.Lz + z][0];
	    resid_i += (-1.0*proj(i*10 + j, kx, ky, kz))
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

#endif //FFT
