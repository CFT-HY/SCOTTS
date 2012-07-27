
/* transport.c
 *
 * Do advection. Two versions, depending on whether the
 * field being advected lives on lattice sites or between them.
 */
#include "hydro.h"



/* Donor cell advection for E in direction dir */
void donor_E_dir(hydro_fields f, int **nb, hydro_params p, int dir) {

  int x, y, z;

  //  halo_field(f.V[dir],p);
  //  halo_field(f.E,p);

  int dx = 0;
  int dy = 0;
  int dz = 0;

  if(dir == 0)
    dx = 1;
  if(dir == 1)
    dy = 1;
  if(dir == 2)
    dz = 1;


  for(x = 0; x < p.Lx; x++) {
    for(y = 0; y < p.Ly; y++) {
      for(z = 0; z < p.Lz; z++) {
    if(f.V[dir][x][y][z] >= 0.0)
      f.F[dir][x][y][z] = f.V[dir][x][y][z]
	*f.E[(x-dx+p.Lx)%p.Lz][(y-dy+p.Ly)%p.Ly][(z-dz+p.Lz)%p.Lz];
    else
      f.F[dir][x][y][z] = f.V[dir][x][y][z]*f.E[x][y][z];
      }
    }
  }

  halo_field(f.F[dir],p);


  for(x = 0; x < p.Lx; x++) {
    for(y = 0; y < p.Ly; y++) {
      for(z = 0; z < p.Lz; z++) {
	f.E[x][y][z] = f.E[x][y][z] 
	  + p.dt*(f.F[dir][x][y][z]
		  - f.F[dir][(x+dx+p.Lx)%p.Lx][(y+dy+p.Ly)%p.Ly][(z+dz+p.Lz)%p.Lz])/p.dx;
      }
    }
  }  

  halo_field(f.E, p);


}


// straight from fortran
void donor_Z_dir(hydro_fields f, int **nb, hydro_params p, int dir) {

  double vc;
  int x, y, z;


  double *Vbody_root = (double *)malloc(p.fieldN*sizeof(double));
  double ***Vbody = make_field(p, Vbody_root);



  // Regenerate fluxes of inertial density
  
  for(x = 0; x < p.Lx; x++) {
    for(y = 0; y < p.Ly; y++) {
      for(z = 0; z < p.Lz; z++) {


	Vbody[x][y][z] =  0.5*(f.V[0][x][y][z] + f.V[0][(x+1)%p.Lx][y][z]);


	  if(Vbody[x][y][z] >= 0.0)
	    f.F[0][x][y][z] = Vbody[x][y][z]*f.Z[dir][x][y][z];
	  else
	    f.F[0][x][y][z] = Vbody[x][y][z]*f.Z[dir][(x+1)%p.Lx][y][z];



	  Vbody[x][y][z] =  0.5*(f.V[1][x][y][z] + f.V[1][x][(y+1)%p.Ly][z]);


	  if(Vbody[x][y][z] >= 0.0)
	    f.F[1][x][y][z] = Vbody[x][y][z]*f.Z[dir][x][y][z];
	  else
	    f.F[1][x][y][z] = Vbody[x][y][z]*f.Z[dir][x][(y+1)%p.Ly][z];


	  Vbody[x][y][z] =  0.5*(f.V[2][x][y][z] + f.V[2][x][y][(z+1)%p.Lz]);


	  if(Vbody[x][y][z] >= 0.0)
	    f.F[2][x][y][z] = Vbody[x][y][z]*f.Z[dir][x][y][z];
	  else
	    f.F[2][x][y][z] = Vbody[x][y][z]*f.Z[dir][x][y][(z+1)%p.Lz];



      }
    }
  }

  

  halo_field(f.F[0],p);
  halo_field(f.F[1],p);
  halo_field(f.F[2],p);

  // Eq 2.11
  for(x = 0; x < p.Lx; x++) {
    for(y = 0; y < p.Ly; y++) {
      for(z = 0; z < p.Lz; z++) {


    f.Z[dir][x][y][z] = f.Z[dir][x][y][z] - p.dt*(f.F[0][x][y][z] - f.F[0][(x-1+p.Lx)%p.Lx][y][z])/p.dx;
    f.Z[dir][x][y][z] = f.Z[dir][x][y][z] - p.dt*(f.F[1][x][y][z] - f.F[1][x][(y-1+p.Ly)%p.Ly][z])/p.dx;
    f.Z[dir][x][y][z] = f.Z[dir][x][y][z] - p.dt*(f.F[2][x][y][z] - f.F[2][x][y][(z-1+p.Lz)%p.Lz])/p.dx;
      }
    }
  }
  halo_field(f.Z[dir], p);

  free_field(p, Vbody, Vbody_root);

}


void advect_E(hydro_fields f, int **nb, hydro_params p) {

  int order; 
  order = 0; // for deterministic results across sites: lrand48() % 3;

  donor_E_dir(f, nb, p, order);
  donor_E_dir(f, nb, p, (order + 1) % 3);
  donor_E_dir(f, nb, p, (order + 2) % 3);

}




void advect_Z(hydro_fields f, int **nb, hydro_params p) {

  int order; 
  order = 0; // lrand48() % 3;

  donor_Z_dir(f, nb, p, order);
  donor_Z_dir(f, nb, p, (order + 1) % 3);
  donor_Z_dir(f, nb, p, (order + 2) % 3);

}
















