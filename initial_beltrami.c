#include "hydro.h"


/*
 *  Beltrami velocity field for an incompressible fluid
 *  Depends on a velocity V0 and a number of 'turns' n
 *  Field is initilized along the z axis
 */
void initial_beltrami(hydro_fields f, hydro_params p,double U0, int n){

    double k0 = 2.0*M_PI*n/p.Lz;
    double eps = 1.0;
    int x,y,z;

    for(x=1; x<=p.slicex; x++) {
		for(y=1; y<=p.slicey; y++) {
			for(z=0; z<p.Lz; z++) {
                f.U[0][x][y][z] = U0 * cos(k0*z);
                f.U[1][x][y][z] = U0 * sin(k0*z);
                f.U[2][x][y][z] = 0;

                f.W[x][y][z] = 1.0 + U0*U0;

                f.V[0][x][y][z] = f.U[0][x][y][z]/(1.0+U0*U0);
                f.V[1][x][y][z] = f.U[1][x][y][z]/(1.0+U0*U0);
                f.V[2][x][y][z] = f.U[2][x][y][z]/(1.0+U0*U0);

                f.E[x][y][z] = eps*f.W[x][y][z];
                f.p[x][y][z] = (1.0/3.0)*eps;
                f.kappa[x][y][z] = 1.0 + (1.0/3.0)/f.W[x][y][z];

			}
		}
	}

	halo_field(f.W, p);
	halo_field(f.E, p);
	halo_field(f.kappa, p);
    halo_field(f.V[0], p);
	halo_field(f.V[1], p);
	halo_field(f.V[2], p);
    halo_field(f.U[0], p);
	halo_field(f.U[1], p);
	halo_field(f.U[2], p);

	double sigmabar;

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
