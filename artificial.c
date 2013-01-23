/* artificial.c
 *
 * Artificial viscosity.
 */
#include "hydro.h"


void artificial_viscosity(hydro_fields f, hydro_params p) {

  int x, y, z;

  /*  
  double ***gradu = make_field(p);

  double Ux, Uxp, Uy, Uyp, Uz, Uzp;


  for(x = 1; x <= p.slicex; x++) {
    for(y = 1; y <= p.slicey; y++) {
      for(z = 0; z < p.Lz; z++) {


	Ux = 0.25*(f.U[0][x][y][z] + f.U[0][x][y+1][z] + f.U[0][x][y][(z+1)%p.Lz]
		   + f.U[0][x][y+1][(z+1)%p.Lz]);

	Uxp = 0.25*(f.U[0][x+1][y][z] + f.U[0][x+1][y+1][z] + f.U[0][x+1][y][(z+1)%p.Lz]
		    + f.U[0][x+1][y+1][(z+1)%p.Lz]);


	Uy = 0.25*(f.U[1][x][y][z] + f.U[1][x+1][y][z] + f.U[1][x][y][(z+1)%p.Lz]
		   + f.U[1][x+1][y][(z+1)%p.Lz]);

	Uyp = 0.25*(f.U[1][x][y+1][z] + f.U[1][x+1][y+1][z] + f.U[1][x][y+1][(z+1)%p.Lz]
		    + f.U[1][x+1][y+1][(z+1)%p.Lz]);


	Uz = 0.25*(f.U[2][x][y][z] + f.U[2][x+1][y][z] + f.U[2][x][y+1][z] + f.U[2][x+1][y+1][z]);

	Uzp = 0.25*(f.U[2][x][y][(z+1)%p.Lz] f.U[2][x+1][y][(z+1)%p.Lz]
		    + f.U[2][x][y+1][(z+1)%p.Lz] + f.U[2][x+1][y+1][(z+1)%p.Lz]);

	gradu[x][y][z] = (Uxp - Ux + Uyp - Uy + Uzp - Uz)/p.dx;
	
      }
    }
  }


  double Umin, Umax, Udix, Udiv;


  for(x = 1; x <= p.slicex; x++) {
    for(y = 1; y <= p.slicey; y++) {
      for(z = 0; z < p.Lz; z++) {
	Umin = minof3(f.U[0][x-1][y][z], f.U[0][x][y][z], f.U[0][x+1][y][z]);
	Umax = maxof3(f.U[0][x-1][y][z], f.U[0][x][y][z], f.U[0][x+1][y][z]);
	Udix = minof2(Umax-U[0][x][y][z], U[0][x][y][z]-Umin);

	Udiv = 0.0;

	if((U[0][x-1][y][z] - U[0][x][y][z])*(U[0][x][y][z]-U[0][x-1][y][z]) > 0.0) {
	  Udiv = Udix;
	} else {
	  Udiv = 0.0;
	}
  */

  double dvx, dvy, dvz;

  double ****Q = make_vector(p);

  for(x = 1; x <= p.slicex; x++) {
    for(y = 1; y <= p.slicey; y++) {
      for(z = 0; z < p.Lz; z++) {

	dvx = f.V[0][x][y][z] - f.V[0][x+1][y][z];
	dvy = f.V[1][x][y][z] - f.V[1][x][y+1][z];
	dvz = f.V[2][x][y][z] - f.V[2][x][y][(z+1)%p.Lz];

	if(dvx < 0)
	  Q[0][x][y][z] = p.Cav*f.E[x][y][z]*dvx*dvx;
	else
	  Q[0][x][y][z] = 0.0;

	if(dvy < 0)
	  Q[1][x][y][z] = p.Cav*f.E[x][y][z]*dvy*dvy;
	else
	  Q[1][x][y][z] = 0.0;

	if(dvz < 0)
	  Q[2][x][y][z] = p.Cav*f.E[x][y][z]*dvz*dvz;
	else
	  Q[2][x][y][z] = 0.0;

	f.E[x][y][z] = f.E[x][y][z] - p.dt*Q[0][x][y][z]*dvx/p.dx;
	f.E[x][y][z] = f.E[x][y][z] - p.dt*Q[1][x][y][z]*dvy/p.dx;
	f.E[x][y][z] = f.E[x][y][z] - p.dt*Q[2][x][y][z]*dvz/p.dx;
      }
    }
  }

  halo_field(Q[0], p);
  halo_field(Q[1], p);
  halo_field(Q[2], p);


  for(x = 1; x <= p.slicex; x++) {
    for(y = 1; y <= p.slicey; y++) {
      for(z = 0; z < p.Lz; z++) {
	f.Z[0][x][y][z] = f.Z[0][x][y][z] - p.dt*(Q[0][x+1][y][z] - Q[0][x][y][z])/p.dx;
	f.Z[1][x][y][z] = f.Z[1][x][y][z] - p.dt*(Q[1][x][y+1][z] - Q[1][x][y][z])/p.dx;
	f.Z[2][x][y][z] = f.Z[2][x][y][z] - p.dt*(Q[2][x][y][(z+1)%p.Lz] - Q[2][x][y][z])/p.dx;
      }
    }
  }

  free_vector(p,Q);


}

