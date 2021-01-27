/** @file initps.c
*
* Random Gaussian initial conditions with specified power spectrum.
* Allows to specifiy either longitudinal waves (div), vortical (rot) or
*
* The main entry point is the function init_ps()
*/
#include "hydro.h"

#if defined(FFT) && ! defined(SCALAR)

/** Lagrange Interpolating formula, for three fixed points
 * Finds the polynomial with the lowest order that interpolates between a set of points
 * \f[ \text{quadratic}(x) = \sum_j y_j \left[\prod_{i \neq j}
 * \frac{(x - x_i)}{x_j - x_i} \right] \f]
 * Previously used in spectrum_interp()
 */
float quadratic(float x1, float x2, float x3,float y1, float y2, float y3, float x) {

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


/** DEBUG : function that writes the power spectra right after initialization
 *
 */
void debug_write_power(hydro_params p, fftwf_complex **in, ptrdiff_t x_start, ptrdiff_t x_thickness) {

	int x, y, z;
	int i, j;
	int true_x, true_y, true_z;
	float kx, ky, kz;

	char fftdest[200];

	// float *product = (float *)malloc(x_thickness*p.Ly*p.Lz*sizeof(float));
	// float *product_div = (float *)malloc(x_thickness*p.Ly*p.Lz*sizeof(float));
	float *product_tot = (float *)malloc(x_thickness*p.Ly*p.Lz*sizeof(float));

	float res_r, res_i, resid_r, resid_i, tot_r, tot_i;

	MPI_Barrier(MPI_COMM_WORLD);
	fprintf(stderr,"rank: %d starting debug PS... (thickness %d, start %d)\n",
	p.rank, x_thickness, x_start);

	MPI_Barrier(MPI_COMM_WORLD);

	for(x=0; x<x_thickness; x++) {
		for(y=0; y<p.Ly; y++) {
			for(z=0; z<p.Lz; z++) {

				if(x+x_start > p.Lx/2) {
					true_x = -(p.Lx - (x+x_start));
				} else {
					true_x = x+x_start;
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

				kx = 2.0*sin(((float)(true_x))*M_PI/(((float)p.Lx)))/p.dx;
				ky = 2.0*sin(((float)(true_y))*M_PI/(((float)p.Ly)))/p.dx;
				kz = 2.0*sin(((float)(true_z))*M_PI/(((float)p.Lz)))/p.dx;

				// product[x*p.Ly*p.Lz + y*p.Lz + z] = 0.0;
				// product_div[x*p.Ly*p.Lz + y*p.Lz + z] = 0.0;
				product_tot[x*p.Ly*p.Lz + y*p.Lz + z] = 0.0;


				for(i=1; i<=3; i++) {

					res_r = 0.0;
					res_i = 0.0;

					resid_r = 0.0;
					resid_i = 0.0;

					tot_r = in[i-1][x*p.Ly*p.Lz + y*p.Lz + z][0];
					tot_i = in[i-1][x*p.Ly*p.Lz + y*p.Lz + z][1];

					/*
					for(j=1; j<=3; j++) {

					// Transverse components
					res_r += proj(i*10 + j, kx, ky, kz)
					*in[j-1][x*p.Ly*p.Lz + y*p.Lz + z][0];
					res_i += proj(i*10 + j, kx, ky, kz)
					*in[j-1][x*p.Ly*p.Lz + y*p.Lz + z][1];

					// And longitudinal...
					if(i == j) {
					resid_r += (1.0 - proj(i*10 + j, kx, ky, kz))
					*in[j-1][x*p.Ly*p.Lz + y*p.Lz + z][0];
					resid_i += (1.0 - proj(i*10 + j, kx, ky, kz))
					*in[j-1][x*p.Ly*p.Lz + y*p.Lz + z][1];
					} else {
					resid_r += (-1.0*proj(i*10 + j, kx, ky, kz))
					*in[j-1][x*p.Ly*p.Lz + y*p.Lz + z][0];
					resid_i += (-1.0*proj(i*10 + j, kx, ky, kz))
					*in[j-1][x*p.Ly*p.Lz + y*p.Lz + z][1];
					}
					}
					*/

					/*
					product[x*p.Ly*p.Lz + y*p.Lz + z] +=
					res_r*res_r + res_i*res_i;

					product_div[x*p.Ly*p.Lz + y*p.Lz + z] +=
					resid_r*resid_r + resid_i*resid_i;
					*/

					product_tot[x*p.Ly*p.Lz + y*p.Lz + z] +=
					tot_r*tot_r + tot_i*tot_i;

				}

			}
		}
	}

	sprintf(fftdest, "test-ps-init.txt");
	MPI_Barrier(MPI_COMM_WORLD);
	if(!p.rank)
	fprintf(stderr,"writing debug PS...\n");
	histogram(p, product_tot, fftdest, x_thickness, x_start);
}


void UtoZ(hydro_fields f, hydro_params p) {

	int i, x, y, z;
	float ubarx, ubary, ubarz;
	float utildex, utildey, utildez;
	float Wfacex, Wfacey, Wfacez;


	// This assumes U has been initialised (but not necessarily haloed)

	// First, then, halo U
	halo_field(f.U[0], p);
	halo_field(f.U[1], p);
	halo_field(f.U[2], p);

	// Calculate Lorentz factor W.
	for(x = 1; x <= p.slicex; x++) {
	  for(y = 1; y <= p.slicey; y++) {
	    for(z = 0; z < p.Lz; z++) {


	      utildex = (f.U[0][x][y][z]
			 + f.U[0][x+1][y][z]
			 + f.U[0][x][y+1][z]
			 + f.U[0][x][y][((z+1)%p.Lz)]
			 + f.U[0][x+1][y][((z+1)%p.Lz)]
			 + f.U[0][x][y+1][((z+1)%p.Lz)]
			 + f.U[0][x+1][y+1][z]
			 + f.U[0][x+1][y+1][((z+1)%p.Lz)]
			 )/8.0;

	      utildey = (f.U[1][x][y][z]
			 + f.U[1][x+1][y][z]
			 + f.U[1][x][y+1][z]
			 + f.U[1][x][y][((z+1)%p.Lz)]
			 + f.U[1][x+1][y][((z+1)%p.Lz)]
			 + f.U[1][x][y+1][((z+1)%p.Lz)]
			 + f.U[1][x+1][y+1][z]
			 + f.U[1][x+1][y+1][((z+1)%p.Lz)]
			 )/8.0;

	      utildez = (f.U[2][x][y][z]
			 + f.U[2][x+1][y][z]
			 + f.U[2][x][y+1][z]
			 + f.U[2][x][y][((z+1)%p.Lz)]
			 + f.U[2][x+1][y][((z+1)%p.Lz)]
			 + f.U[2][x][y+1][((z+1)%p.Lz)]
			 + f.U[2][x+1][y+1][z]
			 + f.U[2][x+1][y+1][((z+1)%p.Lz)]
			 )/8.0;


	      //    if(x < 3 || x > p.N-3)
	      //      fprintf(stderr,"utildex %d = %lf\n", x, utildex);

	      f.W[x][y][z] = sqrt(1.0
				  + utildex*utildex
				  + utildey*utildey
				  + utildez*utildez);
	    }
	  }
	}


	// Work out V from U ensuring accurate centring.

	for(x = 1; x <= p.slicex; x++) {
	  for(y = 1; y <= p.slicey; y++) {
	    for(z = 0; z < p.Lz; z++) {

	      ubarx = (f.U[0][x][y][z]
		       + f.U[0][x][y+1][z]
		       + f.U[0][x][y][((z+1)%p.Lz)]
		       + f.U[0][x][y+1][((z+1)%p.Lz)]
		       )/4.0;

	      ubary = (f.U[1][x][y][z]
		       + f.U[1][x][y+1][z]
		       + f.U[1][x][y][((z+1)%p.Lz)]
		       + f.U[1][x][y+1][((z+1)%p.Lz)]
		       )/4.0;

	      ubarz = (f.U[2][x][y][z]
		       + f.U[2][x][y+1][z]
		       + f.U[2][x][y][((z+1)%p.Lz)]
		       + f.U[2][x][y+1][((z+1)%p.Lz)]
		       )/4.0;

	      Wfacex = sqrt(1.0 + ubarx*ubarx + ubary*ubary + ubarz*ubarz);

	      f.V[0][x][y][z] = ubarx/Wfacex;
	    }
	  }
	}




	// y-cpt
	for(x = 1; x <= p.slicex; x++) {
	  for(y = 1; y <= p.slicey; y++) {
	    for(z = 0; z < p.Lz; z++) {

	      ubarx = (f.U[0][x][y][z]
		       + f.U[0][x+1][y][z]
		       + f.U[0][x][y][((z+1)%p.Lz)]
		       + f.U[0][x+1][y][((z+1)%p.Lz)]
		       )/4.0;

	      ubary = (f.U[1][x][y][z]
		       + f.U[1][x+1][y][z]
		       + f.U[1][x][y][((z+1)%p.Lz)]
		       + f.U[1][x+1][y][((z+1)%p.Lz)]
		       )/4.0;

	      ubarz = (f.U[2][x][y][z]
		       + f.U[2][x+1][y][z]
		       + f.U[2][x][y][((z+1)%p.Lz)]
		       + f.U[2][x+1][y][((z+1)%p.Lz)]
		       )/4.0;

	      Wfacey = sqrt(1.0 + ubarx*ubarx + ubary*ubary + ubarz*ubarz);

	      f.V[1][x][y][z] = ubary/Wfacey;

	    }
	  }
	}

	// z-cpt
	for(x = 1; x <= p.slicex; x++) {
	  for(y = 1; y <= p.slicey; y++) {
	    for(z = 0; z < p.Lz; z++) {

	      ubarx = (f.U[0][x][y][z]
		       + f.U[0][x][y+1][z]
		       + f.U[0][x+1][y][z]
		       + f.U[0][x+1][y+1][z]
		       )/4.0;

	      ubary = (f.U[1][x][y][z]
		       + f.U[1][x][y+1][z]
		       + f.U[1][x+1][y][z]
		       + f.U[1][x+1][y+1][z]
		       )/4.0;

	      ubarz = (f.U[2][x][y][z]
		       + f.U[2][x][y+1][z]
		       + f.U[2][x+1][y][z]
		       + f.U[2][x+1][y+1][z]
		       )/4.0;

	      Wfacez = sqrt(1.0 + ubarx*ubarx + ubary*ubary + ubarz*ubarz);

	      f.V[2][x][y][z] = ubarz/Wfacez;
	    }
	  }
	}






	// Halo W and V
	halo_field(f.W, p);
	halo_field(f.V[0], p);
	halo_field(f.V[1], p);
	halo_field(f.V[2], p);

	float sigmabar;

	// Following assumes that f.kappa has been set
	// (i.e that eq_of_state has been called.)
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


/** Initialize the energy density
 *
 * For a rotational-free flow, the energy density is linked to the velocity field
 * For a divergence-free flow, the energy density is set to be constant
 */
void init_energy(hydro_params p, hydro_fields f, ptrdiff_t x_start, ptrdiff_t x_thickness, int* map, ptrdiff_t alloc_local,float *k_bins, float *pow_bins){
	int x, y, z;
	// TODO : INITPSFILE_ALL

	if(p.initpsfile_type == INITPSFILE_ROT) {
		for(x=1; x<=p.slicex; x++) {
			for(y=1; y<=p.slicey; y++) {
				for(z=0; z<p.Lz; z++) {
					f.E[x][y][z] = 1;
				}
			}
		}

	} else if(p.initpsfile_type == INITPSFILE_DIV) {

		ptrdiff_t n0 = p.Lx;
		ptrdiff_t n1 = p.Ly;
		ptrdiff_t n2 = p.Lz;
		MPI_Status status;

		fftwf_complex *in = fftwf_alloc_complex(alloc_local);
		fftwf_complex *out = fftwf_alloc_complex(alloc_local);
		fftwf_complex *swap_in = fftwf_alloc_complex(alloc_local);
		float *slice = (float *)malloc(x_thickness*p.Ly*p.Lz*sizeof(float));
		float *trim = (float *)malloc(p.slicex*p.slicey*p.Lz*sizeof(float));

		fftwf_plan plan = fftwf_mpi_plan_dft_3d(p.Lx, p.Ly, p.Lz,
			in, out, MPI_COMM_WORLD,
			FFTW_FORWARD, FFTW_ESTIMATE);

		float kx,ky,kz,kmode,ksq;
		int true_x, true_y, true_z;
		fftwf_complex kiVki, vec_i;
		int nx = p.Lx/p.slicex;
		int ny = p.Ly/p.slicey;
		int ry;

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

					kx = 2.0*sin(((float)(true_x))*M_PI/(((float)p.Lx)))/p.dx;
					ky = 2.0*sin(((float)(true_y))*M_PI/(((float)p.Ly)))/p.dx;
					kz = 2.0*sin(((float)(true_z))*M_PI/(((float)p.Lz)))/p.dx;

					kmode = sqrt(kx*kx+ky*ky+kz*kz);
					ksq = kx*kx+ky*ky+kz*kz;

					spectrum_interp(ksq, p,&(vec_i),k_bins, pow_bins, p.initpsbins);
					kiVki[0] = kx * vec_i[0];
					kiVki[1] = kx * vec_i[1];
					spectrum_interp(ksq, p,&(vec_i),k_bins, pow_bins, p.initpsbins);
					kiVki[0] += ky * vec_i[0];
					kiVki[1] += ky * vec_i[1];
					spectrum_interp(ksq, p,&(vec_i),k_bins, pow_bins, p.initpsbins);
					kiVki[0] += kz * vec_i[0];
					kiVki[1] += kz * vec_i[1];

					if (kmode != 0){
						in[(x-x_start)*p.Ly*p.Lz + y*p.Lz + z][0] = 4 / sqrt(3) * kiVki[1] / kmode;
						in[(x-x_start)*p.Ly*p.Lz + y*p.Lz + z][1] = - 4 / sqrt(3) * kiVki[0] / kmode;

					}

				}
			}
		}
		MPI_Barrier(MPI_COMM_WORLD);

		// These are the minimum and maximum (INCLUSIVE) x-slice coordinates needed by this node
		int min_needed = p.Lx-(x_start+x_thickness)+1;
		int max_needed = p.Lx-x_start;

		// At this point each node knows what it ought to have, but we need to 'hermitianise' it
		// Need to do swapping
		int j;

		for(j=0; j<=p.Lx; j++) {
			if(j>= x_start && j < x_start+x_thickness) {
				if(map[((p.Lx-j) % p.Lx)] == p.rank) {
				} else {
					MPI_Send(&(((float *)in)[((j-x_start)%p.Lx)*p.Ly*p.Lz*2]),
					2*p.Ly*p.Lz, MPI_FLOAT,
					map[((p.Lx-j) % p.Lx)], j, MPI_COMM_WORLD);
				}
			}

			if(j>=min_needed && j<=max_needed) {
				if(map[(j % p.Lx)] == p.rank) {
					memcpy(&(((fftwf_complex *)swap_in)[(j-min_needed)*p.Ly*p.Lz]),
					&(((fftwf_complex *)in)[((j-x_start)%p.Lx)*p.Ly*p.Lz]),
					p.Ly*p.Lz*sizeof(fftwf_complex));


				} else {
					MPI_Recv(&(((float *)swap_in)[(j-min_needed)*p.Ly*p.Lz*2]),
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
						in[(x-x_start)*p.Ly*p.Lz + y*p.Lz + z][0] \
							= sqrt(swap_in[swap_row*p.Ly*p.Lz
							+ ((p.Ly-y)%p.Ly)*p.Lz + ((p.Lz-z)%p.Lz)][0]
							*swap_in[swap_row*p.Ly*p.Lz
							+ ((p.Ly-y)%p.Ly)*p.Lz + ((p.Lz-z)%p.Lz)][0]
							+ swap_in[swap_row*p.Ly*p.Lz
							+ ((p.Ly-y)%p.Ly)*p.Lz + ((p.Lz-z)%p.Lz)][1]
							*swap_in[swap_row*p.Ly*p.Lz
							+ ((p.Ly-y)%p.Ly)*p.Lz + ((p.Lz-z)%p.Lz)][1]);

						in[(x-x_start)*p.Ly*p.Lz + y*p.Lz + z][1] = 0.0;


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
						in[(x-x_start)*p.Ly*p.Lz + y*p.Lz + z][0] = \
							in[(x-x_start)*p.Ly*p.Lz
							+ ((p.Ly-y)%p.Ly)*p.Lz + ((p.Lz-z)%p.Lz)][0];
						in[(x-x_start)*p.Ly*p.Lz + y*p.Lz + z][1] = \
							-1.0*in[(x-x_start)*p.Ly*p.Lz
							+ ((p.Ly-y)%p.Ly)*p.Lz + ((p.Lz-z)%p.Lz)][1];

					} else if(x >= p.Lx/2) {
						/* in this case we conjugate the stuff in the lower
						slices, using swap_row to work out where it was put.
						*/

						in[(x-x_start)*p.Ly*p.Lz + y*p.Lz + z][0] = \
						1.0*swap_in[swap_row*p.Ly*p.Lz
						+ ((p.Ly-y)%p.Ly)*p.Lz + ((p.Lz-z)%p.Lz)][0];
						in[(x-x_start)*p.Ly*p.Lz + y*p.Lz + z][1] = \
						-1.0*swap_in[swap_row*p.Ly*p.Lz
						+ ((p.Ly-y)%p.Ly)*p.Lz + ((p.Lz-z)%p.Lz)][1];

					}
					// otherwise do nothing


				}
			}
		}

		MPI_Barrier(MPI_COMM_WORLD);

		// Calculate the inverse fourier transform for the energy
		fftwf_execute(plan);
		MPI_Barrier(MPI_COMM_WORLD);

		// prepare slice
		for(x=0; x<x_thickness; x++) {
			for(y=0; y<p.Ly; y++) {
				for(z=0; z<p.Lz; z++) {
					slice[x*p.Ly*p.Lz + y*p.Lz + z] = out[x*p.Ly*p.Lz + y*p.Lz + z][0];
				}
			}
		}

		MPI_Barrier(MPI_COMM_WORLD);

		for(x=0; x<p.Lx; x++) {
			// Check whether we are the source node for this slab
			if(map[x] == p.rank) {
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
					f.E[x][y][z] = exp(trim[(x-1)*p.slicey*p.Lz + (y-1)*p.Lz + z]);
				}
			}
		}

		free(slice);
		free(trim);
		fftwf_destroy_plan(plan);
		fftwf_free(in);
		fftwf_free(out);
	}
	halo_field(f.E, p);

}

/** Projects a 3-dimensional complex vector field in Fourier space onto either
 * longitudinal or rotational components (or not at all)
 *
 * - The rotational projection is given by
 * \f$ v_i^\perp = P_{ij} v_j \f$
 * where \f$ P_{ij} = \delta_{ij} - \hat{k}_i \hat{k}_j \f$
 * and \f$ \hat{k} \f$ is a unit vector in the \f$ k \f$ direction.
 * - The longitudinal projection is given by \f$ v_i^\parallel = (\delta_{ij}
 * - P_{ij}) v_j = \hat{k}_i \hat{k}_j v_j \f$
 *
 * **Note that** each degree of freedom of the Fourier velocity field follows
 * the target power spectrum. Hence for the total field, one has to renormalize
 * by dividing the amplitude by :
 * - \f$ \sqrt{3} \f$ in case of no projection
 * - \f$ \sqrt{2} \f$ for the rotational projection
 * - \f$ 1 \f$ for the divergent projection
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
	float dofs = sqrt(3.0);

	if(p.initpsfile_type == INITPSFILE_ROT) {
		dofs = sqrt(2.0);
	} else if(p.initpsfile_type == INITPSFILE_DIV) {
		dofs = sqrt(1.0);
	}

	// Renormalize the power by taking into account the various dofs
	for(x=0;x<x_thickness;x++) {
		for(y=0;y<p.Ly;y++) {
			for(z=0;z<p.Lz;z++) {
				in[0][x*p.Ly*p.Lz + y*p.Lz + z][0] /= dofs;
				in[1][x*p.Ly*p.Lz + y*p.Lz + z][0] /= dofs;
				in[2][x*p.Ly*p.Lz + y*p.Lz + z][0] /= dofs;

				in[0][x*p.Ly*p.Lz + y*p.Lz + z][1] /= dofs;
				in[1][x*p.Ly*p.Lz + y*p.Lz + z][1] /= dofs;
				in[2][x*p.Ly*p.Lz + y*p.Lz + z][1] /= dofs;
			}
		}
	}

	// TODO : this for-loop was present when I took other the code
	// Is it still necessary ? One run seems to be enough
	for(reruns = 0; reruns < times; reruns++) {
		resid = 0.0;
		stuff = 0.0;

		// Iterate over the slice
		for(x=0;x<x_thickness;x++) {
			for(y=0;y<p.Ly;y++) {
				for(z=0;z<p.Lz;z++) {

					// Computes the value of the wave-vector
					if(x+shift_x > p.Lx/2) {
						true_x = -(p.Lx - (x+shift_x));
					} else {
						true_x = x+shift_x;
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

					kx = 2.0*sin(((float)(true_x))*M_PI/(((float)p.Lx)))/p.dx;
					ky = 2.0*sin(((float)(true_y))*M_PI/(((float)p.Ly)))/p.dx;
					kz = 2.0*sin(((float)(true_z))*M_PI/(((float)p.Lz)))/p.dx;

					in_proj_re[0] = 0.0;
					in_proj_re[1] = 0.0;
					in_proj_re[2] = 0.0;

					in_proj_im[0] = 0.0;
					in_proj_im[1] = 0.0;
					in_proj_im[2] = 0.0;

					if(p.initpsfile_type == INITPSFILE_ALL) {
					// do nothing
					} else {

						for(i=1; i<=3; i++) {

							res_re = 0.0;
							res_im = 0.0;

							for(j=1; j<=3; j++) {

								if(p.initpsfile_type == INITPSFILE_ROT) {
									//  Rot?
									// v_i^\perp = P_{ij} v_j
									// where P_{ij} = \delta_{ij} - \hat{k}_i \hat{k}_j
									// and $\hat{k}$ is a unit vector in the $k$ direction.
									in_proj_re[i-1] += proj(i*10 + j, kx, ky, kz)*in[j-1][x*p.Ly*p.Lz + y*p.Lz + z][0];
									in_proj_im[i-1] += proj(i*10 + j, kx, ky, kz)*in[j-1][x*p.Ly*p.Lz + y*p.Lz + z][1];


									// this is delta_{ij} - P_{ij}V_j
									// TODO : This if-statement could be replaced with a simpler
									// res_re += ((i==j)-proj(i*10 + j, kx, ky, kz))*in[j-1][x*p.Ly*p.Lz + y*p.Lz + z][0];
									// res_rm += ((i==j)-proj(i*10 + j, kx, ky, kz))*in[j-1][x*p.Ly*p.Lz + y*p.Lz + z][1];
									if(i == j) {
										res_re += (1.0-proj(i*10 + j, kx, ky, kz))*in[j-1][x*p.Ly*p.Lz + y*p.Lz + z][0];
										res_im += (1.0-proj(i*10 + j, kx, ky, kz))*in[j-1][x*p.Ly*p.Lz + y*p.Lz + z][1];
									} else {
										res_re += (-1.0*proj(i*10 + j, kx, ky, kz))*in[j-1][x*p.Ly*p.Lz + y*p.Lz + z][0];
										res_im += (-1.0*proj(i*10 + j, kx, ky, kz))*in[j-1][x*p.Ly*p.Lz + y*p.Lz + z][1];
									}


								} else if(p.initpsfile_type == INITPSFILE_DIV) {

									// Div?
									// v_i^\parallel = (\delta_{ij} - P_{ij}) v_j
									if(i == j) {
										in_proj_re[i-1] += (1.0-proj(i*10 + j, kx, ky, kz))*in[j-1][x*p.Ly*p.Lz + y*p.Lz + z][0];
										in_proj_im[i-1] += (1.0-proj(i*10 + j, kx, ky, kz))*in[j-1][x*p.Ly*p.Lz + y*p.Lz + z][1];
									} else {
										in_proj_re[i-1] += (-1.0*proj(i*10 + j, kx, ky, kz))*in[j-1][x*p.Ly*p.Lz + y*p.Lz + z][0];
										in_proj_im[i-1] += (-1.0*proj(i*10 + j, kx, ky, kz))*in[j-1][x*p.Ly*p.Lz + y*p.Lz + z][1];
									}
									// #endif

									res_re += proj(i*10 + j, kx, ky, kz)*in[j-1][x*p.Ly*p.Lz + y*p.Lz + z][0];
									res_im += proj(i*10 + j, kx, ky, kz)*in[j-1][x*p.Ly*p.Lz + y*p.Lz + z][1];

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
		}

		if(!p.rank)
		fprintf(stderr,"run %d, stuff is %g, resid is now %g\n", reruns, stuff, resid);
	}
}


/** Box-Muller Gaussian random number generator
 *
 * Produces random numbers following Gaussian statistics starting with
 * \f$ u_1, u_2 \f$ uniformly distributed in \f$[0, 1]\f$ :
 *
 * \f[ Z_0 = \sqrt{-2 \ln u_0} \cos(2\pi u_1)\f]
 * \f[ Z_1 = \sqrt{-2 \ln u_0} \sin(2\pi u_1)\f]
 *
 * \f$Z_0, Z_1\f$ are independent random variables with a standard normal distribution
 * We take \f$Z_1\f$ and rescale using the mean and standard deviation
 *
 * We throw away \f$ Z_0 \f$
 */
float get_normal(float mean, float dev) {
	float u1, u2;

	u1 = drand48();
	u2 = drand48();

	return mean + dev*(sqrt(-2.0*log(u1))*sin(2.0*M_PI*u2));
}


/** Do a linear interpolation to compute the initial power spectrum.
*
* A key assumption is that the k-modes are _evenly distributed_ throughout
* each bin, which isn't true! And increasing the order of the interpolation
* will not really deal with this systematic error (I have tried quadratic).
*
* If this causes undue concern, there are two options as to how to fix it:
* 1. When saving power spectrums, store the mean $k$ in each bin
*    as well as the central $k$ for each bin, and use the mean $k$ (at least)
*    when using power spectra as inputs.
* 2. Recompute the mean $k$ here when building out the power spectrum again.
*    That would make it harder to change lattice volumes between runs, as
*    the relevant mode counts per bin for making the interpolation
*    of the spectral density are for the original simulation
*    that generated the input PS.
*/
void spectrum_interp(float ksq, hydro_params p, fftwf_complex *res, float *k_bins, float *pow_bins, int n_bins) {

	float phase;

	float L = p.Lx;


	float amp, amp_last, amp_this, amp_next, grad;
	int i;


	amp_last = 0.0;
	amp_this = 0.0;
	amp_next = 0.0;
	amp = 0.0;

	// bin coordinates are midpoints
	// TODO : Could be replaced by :
	// k_bins[1] - k_bins[0]
	float dk = 2.0*k_bins[0];

	float kmode = sqrt(fabs(ksq));
	int whichbin = (int)floor(kmode/dk);

	// Interpolates the value of pow_bins at kmode at first order
	if(whichbin > 0 && whichbin < (n_bins-1)) {
		amp_this = sqrt(pow_bins[whichbin]);
		amp_next = sqrt(pow_bins[whichbin+1]);
		grad = (amp_next - amp_this)/dk;
		amp = amp_this + (kmode - ((float)(whichbin))*dk)*grad;
	}

	// Draws a random gaussian number with variance amp
	(*res)[0] = get_normal(0.0, amp);
	(*res)[1] = (*res)[0];

	phase = drand48();

	// Random phase
	(*res)[0] = (*res)[0]*cos(2.0*M_PI*phase);
	(*res)[1] = (*res)[1]*sin(2.0*M_PI*phase);
}


/** OBSOLETE : Hard-coded analytical formula to initialize the power_spectrum
 *
 * Among other things, does not accept external parameters
 */
void spectrum(float ksq, hydro_params p, fftwf_complex *res){

	float phase;
	float L = p.Lx;


	if(fabs(ksq) < 0.000001 || fabs(ksq) > p.initcutoff) {
		(*res)[0] = 0.0;
		(*res)[1] = 0.0;
	}
	else {


		//  (*res)[0] = p.initnorm*exp(-0.25*sqrt(ksq)*p.dx*((float)L));
		(*res)[0] = get_normal(0.0, p.initnorm*exp(-0.25*sqrt(ksq)*p.dx*((float)L)));
		// (*res)[0] = get_normal(0.0, exp(-1.0*sqrt(ksq)*p.dx*p.initlength));

		// (*res)[0] = get_normal(0.0, cherian_spectrum(ksq,p));
		(*res)[1] = (*res)[0];

		phase = drand48();

		// Random phase
		(*res)[0] = (*res)[0]*cos(2.0*M_PI*phase);
		(*res)[1] = (*res)[1]*sin(2.0*M_PI*phase);

	}
}


/** Initialize the simulation with a gaussian velocity field following
 * a specified power spectrum
 *
 * The code is organized as follows :
 * 1. Set the scalar field in the broken phase / true vacuum
 * 2. Allocate a 3-dimensional complex grid for FFT spread accross the multiple cores
 * 3. Builds a map storing which core contains which slice
 * 4. Reads the input file
 * 5. Fills in random gaussian velocity fields in Fourier space
 *
 * **Note that** the FFTs are performed using a library called FFTW which distributes
 * 3-dimensional grids on multiple cores very differently. Hence the very tedious
 * process of mapping/swapping between the FFTW distribution and the distribution
 * used in the rest of the code
 */
void init_ps(hydro_fields f, hydro_params p, float ****field) {

	printf0(p, "Loading initial power spectrum from %s\n", p.initpsfile);

	ptrdiff_t x_thickness, x_start, alloc_local;
	ptrdiff_t n0 = p.Lx;
	ptrdiff_t n1 = p.Ly;
	ptrdiff_t n2 = p.Lz;


	MPI_Status status;

	float kmode, modsq;
	int x, y, z;
	int i, j;
	float *trim = (float *)malloc(p.slicex*p.slicey*p.Lz*sizeof(float));

	/*
	 * 1. Initialize the field in the broken phase
	 */
	for(x = 1; x <= p.slicex; x++) {
		for(y = 1; y <= p.slicey; y++) {
			for(z = 0; z < p.Lz; z++) {
				f.phi[x][y][z] = p.phi_0;
			}
		}
	}


	/*
	 * 2. Allocates a 3-dimensional complex grid for FFT
	 */
	fftwf_mpi_init();
	alloc_local = fftwf_mpi_local_size_3d(n0, n1, n2,
	   MPI_COMM_WORLD,
	   &x_thickness, &x_start);


	/*
	 *3. Builds a map storing which core contains which slice
	 */
	int *thicknesses = malloc(p.size*sizeof(int));
	int *starts = malloc(p.size*sizeof(int));
	int *map = malloc(p.Lx*sizeof(int));

	// If not the main core, then send your caracteristics to the main core
	if(p.rank) {
		MPI_Send(&x_thickness, 1, MPI_INTEGER, 0, p.rank, MPI_COMM_WORLD);
		MPI_Send(&x_start, 1, MPI_INTEGER, 0, p.rank, MPI_COMM_WORLD);

	}
	// If the main core, receive all the caracteristics
	else {
		thicknesses[0] = x_thickness;
		starts[0] = x_start;

		// Receive MPI messages from other cores
		for(i=1; i<p.size; i++) {
			MPI_Recv(&thicknesses[i], 1, MPI_INTEGER, i, i, MPI_COMM_WORLD, &status);
			MPI_Recv(&starts[i], 1, MPI_INTEGER, i, i, MPI_COMM_WORLD, &status);
		}

		// Computes a map accross the x-axis for the location of each slice
		for(x=0; x<p.Lx; x++) {
			map[x] = -1;
			for(i=0; i<p.size; i++) {
				if(starts[i] <= x && (starts[i] + thicknesses[i]) > x) {
					map[x] = i;
					break;
				}
			}
			// Error message if no one takes responsibility for this slice
			if(map[x] == -1) {
				fprintf(stderr,"cannot find a node responsible for x=%d!\n", x);
				die(-99);
			}

		}
	}

	// Broadcasts the map to the other cores
	printf0(p, "broadcasting our results\n");
	MPI_Bcast(map, p.Lx, MPI_INTEGER, 0, MPI_COMM_WORLD);

	/*
	 * 4. Read the input file and fills in the 3d grids accordingly
	 */
	printf0(p, "now to layout\n");

	// 3-dimensional vector field
	fftwf_complex **in = (fftwf_complex **)malloc(3*sizeof(fftwf_complex *));
		in[0] = fftwf_alloc_complex(alloc_local);
		in[1] = fftwf_alloc_complex(alloc_local);
		in[2] = fftwf_alloc_complex(alloc_local);

	// 3-dimensional vector field
	fftwf_complex **swap_in = (fftwf_complex **)malloc(3*sizeof(fftwf_complex *));
	swap_in[0] = fftwf_alloc_complex(x_thickness*p.Ly*p.Lz);
	swap_in[1] = fftwf_alloc_complex(x_thickness*p.Ly*p.Lz);
	swap_in[2] = fftwf_alloc_complex(x_thickness*p.Ly*p.Lz);

	// 3-dimensional vector field
	fftwf_complex **out = (fftwf_complex **)malloc(3*sizeof(fftwf_complex *));
	out[0] = fftwf_alloc_complex(alloc_local);
	out[1] = fftwf_alloc_complex(alloc_local);
	out[2] = fftwf_alloc_complex(alloc_local);

	int nx = p.Lx/p.slicex;
	int ny = p.Ly/p.slicey;

	int ry;


	float ksq;
	int true_x, true_y, true_z;
	// Different values for the wave-vector
	float *k_bins = (float *)malloc(p.initpsbins*sizeof(float));
	// PMean power per cell
	float *pow_bins = (float *)malloc(p.initpsbins*sizeof(float));
	// Volume of the shell [k, k+dk]
	int in_bin;
	int items_read;
	float fudge;

	// Opens the input file
	FILE *fp = fopen(p.initpsfile, "r");

	// Reads the file line by line
	for(i = 0; i < p.initpsbins; i++) {

		// If the file ends too soon, raise an error
		if(feof(fp)) {
			printf0(p, "Fewer bins in file (%d) than initpsbins suggested (%d)\n",
				i, p.initpsbins);
			die(100);
		}

		// Reads the line, and stores the three columns
		items_read = fscanf(fp, "%f%g%d", &k_bins[i], &pow_bins[i], &in_bin);

		// If the number of columns is wrong, raise an error
		if(items_read != 3) {
			printf0(p,
				"Initpsfile %s: row %d: did not read a full line (3 items), "
				"just %d\n",
				p.initpsfile, i, items_read);
			die(100);
		}


		// If the volume if 0, then the power is 0
		if(in_bin == 0) {
			pow_bins[i] = 0.0;
		}
		// Else, divide by the volume to obtain the mean power per cell
		else {
			pow_bins[i] /= ((float)(((long int)in_bin)*((long int)(i+1))));
		}

		if(isnan(pow_bins[i])) {
			printf0(p, "Bin %d is a NaN\n", i);
		}
	}
	fclose(fp);

	// Corrects for the volume of the shell if dk is different between
	// the input file and the grid
	fudge = (2.0*M_PI/((float)p.Lx))/((k_bins[1]-k_bins[0])*p.dx);
	if(fabs(fudge - 1.0) > 1e-6) {
		printf0(p, "Applying fudge factor %g^3 to power spectrum: "
			"seems like dk or volume or dx is different\n",
			fudge);

		for(i = 0; i < p.initpsbins; i++) {
			pow_bins[i] *= fudge*fudge*fudge;
		}
	}

	/*
	 * 5. Fills in random gaussian velocity fields in Fourier space
	 *
	 * Spans the slice contained in this core
	 */
	for(x=x_start;x<x_start+x_thickness;x++) {
		for(y=0;y<p.Ly;y++) {
			for(z=0;z<p.Lz;z++) {

				/** TODO : Since we only need the absolute value of true_x
				 * the if-statement could be replaced by :
				 * true_x = min(x, p.Lx - x)
				 * though a min function is not defined in vanilla C
				 */
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

				// Computes the magnitude of the wavevector on site (x, y, z)
				ksq  = (2.0 - 2.0*cos(((float)true_x)*2.0*M_PI/(((float)p.Lx))))/(p.dx*p.dx);
				ksq += (2.0 - 2.0*cos(((float)true_y)*2.0*M_PI/(((float)p.Ly))))/(p.dx*p.dx);
				ksq += (2.0 - 2.0*cos(((float)true_z)*2.0*M_PI/(((float)p.Lz))))/(p.dx*p.dx);


				for(i=0; i<3; i++) {
					// Draws a gaussian random number with variance given by the
					// interpolated spectrum
					spectrum_interp(ksq, p,
						&(in[i][(x-x_start)*p.Ly*p.Lz + y*p.Lz + z]),
						k_bins, pow_bins, p.initpsbins);
				}
			}
		}
	}

	// Projects the velocity field on longitudinal or vortical comp.
	// Important even for INITPSFILE_ALL
	project_down(p, in, x_start, x_thickness, 3);

	// Debug function that writes the power spectrum
	debug_write_power(p, in, x_start, x_thickness);

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


	//  printf0(p, "Done FFTing, example value %g %g\n", out[0][0][0], out[0][0][1]);

	/*
	for(x=x_start; x<(x_start+x_thickness); x++) {
	fprintf(stderr, "node %d, x %d:\t->\t%10.3g %10.3g\t\n", p.rank, x, out[0][(x-x_start)*p.Ly*p.Lz][0], out[0][(x-x_start)*p.Ly*p.Lz][1]);
	}

	*/

	float *slice = (float *)malloc(x_thickness*p.Ly*p.Lz*sizeof(float));


	//  float maximag = 0.0;
	for(i=0; i<3; i++) {

		// prepare slice

		for(x=0; x<x_thickness; x++) {
			for(y=0; y<p.Ly; y++) {
				for(z=0; z<p.Lz; z++) {
					/*	  if(fabs(out[i][x*p.Ly*p.Lz + y*p.Lz + z][1]) > fabs(maximag)) {
					maximag = out[i][x*p.Ly*p.Lz + y*p.Lz + z][1];
					} */
					slice[x*p.Ly*p.Lz + y*p.Lz + z] = out[i][x*p.Ly*p.Lz + y*p.Lz + z][0];
				}
			}
		}

		//    fprintf(stderr,"maximag is %g\n", maximag);
		MPI_Barrier(MPI_COMM_WORLD);


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

	printf0(p,"Start the energy density initialization\n");
	init_energy(p, f, x_start, x_thickness, map, alloc_local,k_bins,pow_bins);
	printf0(p,"Energy density initialized\n");

	free(map);
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
