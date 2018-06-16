
/** @file transport.c
 *
 * This file contains the functions that perform advection. 
 *
 * We have two versions, depending on whether the field being advected
 * lives on lattice sites or between them.
 *
 * Currently we use simple donor cell advection. We may need to look
 * into fancier advection methods to avoid getting strong shocks.
 */
#include "hydro.h"



/** Donor cell advection for `f.E` in direction `dir`.
 *
 * Does nothing ifdef `SCALAR`. 
 *
 * Simple donor cell i.e the value of `f.E` at the boundary is taken
 * to be the same as the value in the body of the cell. It may be more
 * consistent to swap to piecewise-linear advection.
 */
void donor_E_dir(hydro_fields f, hydro_params p, int dir) {

#ifndef SCALAR
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


  for(x = 1; x <= p.slicex; x++) {
    for(y = 1; y <= p.slicey; y++) {
      for(z = 0; z < p.Lz; z++) {
    if(f.V[dir][x][y][z] >= 0.0)
      f.F[dir][x][y][z] = f.V[dir][x][y][z]
	*f.E[x-dx][y-dy][(z-dz+p.Lz)%p.Lz];
    else
      f.F[dir][x][y][z] = f.V[dir][x][y][z]*f.E[x][y][z];
      }
    }
  }

  halo_field(f.F[dir],p);


  for(x = 1; x <= p.slicex; x++) {
    for(y = 1; y <= p.slicey; y++) {
      for(z = 0; z < p.Lz; z++) {
	f.E[x][y][z] = f.E[x][y][z] 
	  + p.dt*(f.F[dir][x][y][z]
		  - f.F[dir][x+dx][y+dy][(z+dz+p.Lz)%p.Lz])/p.dx;
      }
    }
  }  

  halo_field(f.E, p);
#endif

}


/** Donor cell advection for `f.Z` in direction `dir`.
 *
 * Does nothing ifdef `SCALAR`. 
 *
 * 
 * Simple donor cell i.e the value of `f.Z` at the boundary is taken
 * to be the same as the value in the body of the cell. It may be more
 * consistent to swap to piecewise-linear advection.
 *
 * Needs to calculate velocity in the body of the cell as that is
 * where the interfaces for the momentum flux live.
 *
 * Note here that `dir` corresponds to the directional component \f$ i
 * \f$ of the momentum density \f$ Z_i \f$, not the direction that we
 * perform advection. This is perhaps confusing and could be changed.
 */
void donor_Z_dir(hydro_fields f, hydro_params p, int dir) {

#ifndef SCALAR
  Real vc;
  int x, y, z;


  Real ***Vbody = make_field(p);



  // Regenerate fluxes of inertial density
  
  for(x = 1; x <= p.slicex; x++) {
    for(y = 1; y <= p.slicey; y++) {
      for(z = 0; z < p.Lz; z++) {


	Vbody[x][y][z] =  0.5*(f.V[0][x][y][z] + f.V[0][x+1][y][z]);


	  if(Vbody[x][y][z] >= 0.0)
	    f.F[0][x][y][z] = Vbody[x][y][z]*f.Z[dir][x][y][z];
	  else
	    f.F[0][x][y][z] = Vbody[x][y][z]*f.Z[dir][x+1][y][z];



	  Vbody[x][y][z] =  0.5*(f.V[1][x][y][z] + f.V[1][x][y+1][z]);


	  if(Vbody[x][y][z] >= 0.0)
	    f.F[1][x][y][z] = Vbody[x][y][z]*f.Z[dir][x][y][z];
	  else
	    f.F[1][x][y][z] = Vbody[x][y][z]*f.Z[dir][x][y+1][z];


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
  for(x = 1; x <= p.slicex; x++) {
    for(y = 1; y <= p.slicey; y++) {
      for(z = 0; z < p.Lz; z++) {


	f.Z[dir][x][y][z] = f.Z[dir][x][y][z]
	  - p.dt*(f.F[0][x][y][z] - f.F[0][x-1][y][z])/p.dx;
	f.Z[dir][x][y][z] = f.Z[dir][x][y][z]
	  - p.dt*(f.F[1][x][y][z] - f.F[1][x][y-1][z])/p.dx;
	f.Z[dir][x][y][z] = f.Z[dir][x][y][z]
	  - p.dt*(f.F[2][x][y][z] - f.F[2][x][y][(z-1+p.Lz)%p.Lz])/p.dx;
      }
    }
  }
  halo_field(f.Z[dir], p);

  free_field(p, Vbody);

#endif
}



/** Donor cell advection for `E` in all directions.
 *
 * Currently the order in which each direction is advected is
 * fixed. This can be changed to be random by switching to the code in
 * the comment. Note that this is not what the Grant and Wilson
 * suggest, rather they state the update should be split into two half
 * steps, and the order should be reversed for each half.
 */
void advect_E(hydro_fields f, hydro_params p) {

  int order; 
  order = 0; // for deterministic results across sites: lrand48() % 3;

  donor_E_dir(f, p, order);
  donor_E_dir(f, p, (order + 1) % 3);
  donor_E_dir(f, p, (order + 2) % 3);

}



/** Donor cell advection for `f.Z` in all directions.
 *
 * Currently the order in which component of `f.Z` is advected is
 * fixed. This can be changed to be random by switching to the code in
 * the comment. Note that this is not what the Grant and Wilson
 * suggest, rather they state the update should be split into two half
 * steps, and the order of advection should be reversed for each half.
 */
void advect_Z(hydro_fields f, hydro_params p) {

  int order; 
  order = 0; // lrand48() % 3;

  donor_Z_dir(f, p, order);
  donor_Z_dir(f, p, (order + 1) % 3);
  donor_Z_dir(f, p, (order + 2) % 3);

}
















