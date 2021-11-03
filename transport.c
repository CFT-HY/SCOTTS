
/** @file transport.c
 *
 * This file contains the functions that perform advection. 
 *
 * We have two versions, depending on whether the field being advected
 * lives on lattice sites or between them.
 *
 * Currently we use simple donor cell advection, with a compiler flag
 * enabling Van Leer advection. We may need to look into fancier
 * advection methods to avoid getting strong shocks.
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
 * Advects each component of Z in the direction `dir`.
 */
void donor_Z_dir(hydro_fields f, hydro_params p, int dir) {

#ifndef SCALAR

  int x, y, z, i;

  int dx = 0;
  int dy = 0;
  int dz = 0;

  if(dir == 0){
    dx = 1;
  }
  if(dir == 1){
    dy = 1;
  }
  if(dir == 2){
    dz = 1;
  }

  float ***Vface = make_field(p);

  
  for(x = 1; x <= p.slicex; x++) {
    for(y = 1; y <= p.slicey; y++) {
      for(z = 0; z < p.Lz; z++) {

	// We need to find the flux between nodes which requires us to find the
	// velocity on the appropriate 'face' between nodes.

	if(dir == 0){
	  Vface[x][y][z] =  0.125*((f.V[0][x][y][z] + f.V[0][x-1][y][z])
				   + (f.V[0][x][y-1][z] + f.V[0][x-1][y-1][z])
				   + (f.V[0][x][y][(z-1+p.Lz)%p.Lz]
				      + f.V[0][x-1][y][(z-1+p.Lz)%p.Lz])
				   + (f.V[0][x][y-1][(z-1+p.Lz)%p.Lz]
				      + f.V[0][x-1][y-1][(z-1+p.Lz)%p.Lz]));
	}
	else if(dir == 1){
	  Vface[x][y][z] = 0.125*((f.V[1][x][y][z] + f.V[1][x][y-1][z])
				  + (f.V[1][x-1][y][z] + f.V[1][x-1][y-1][z])
				  + (f.V[1][x][y][(z-1+p.Lz)%p.Lz]
				     + f.V[1][x][y-1][(z-1+p.Lz)%p.Lz])
				  + (f.V[1][x-1][y][(z-1+p.Lz)%p.Lz]
				     + f.V[1][x-1][y-1][(z-1+p.Lz)%p.Lz]));
	}
	else if(dir == 2){
	  Vface[x][y][z] = 0.125*((  f.V[2][x][y][z]
				     + f.V[2][x][y][(z-1+p.Lz)%p.Lz])
				  + (f.V[2][x-1][y][z]
				     + f.V[2][x-1][y][(z-1+p.Lz)%p.Lz])
				  + (f.V[2][x][y-1][z]
				     + f.V[2][x][y-1][(z-1+p.Lz)%p.Lz])
				  + (f.V[2][x-1][y-1][z]
				     + f.V[2][x-1][y-1][(z-1+p.Lz)%p.Lz]));
	}


	for(i = 0; i < 3; i++){
	  if(Vface[x][y][z] >= 0.0){
	    f.F[i][x][y][z] = Vface[x][y][z]*f.Z[i][x-dx][y-dy][(z-dz+p.Lz)%p.Lz];
	  }
	  else{
	    f.F[i][x][y][z] = Vface[x][y][z]*f.Z[i][x][y][z];
	  }
	}
      }
    }
  }

  halo_field(f.F[0],p);
  halo_field(f.F[1],p);
  halo_field(f.F[2],p);

  for(x = 1; x <= p.slicex; x++) {
    for(y = 1; y <= p.slicey; y++) {
      for(z = 0; z < p.Lz; z++) {
	for(i = 0; i < 3; i++){
	  f.Z[i][x][y][z] = (f.Z[i][x][y][z]
			     - p.dt*(f.F[i][x+dx][y+dy][(z+dz)%p.Lz]
				     - f.F[i][x][y][z])/p.dx);
	}
      }
    }
  }
  halo_field(f.Z[0], p);
  halo_field(f.Z[1], p);
  halo_field(f.Z[2], p);

  
  free_field(p, Vface);


#endif
}


/** Advection for `f.E` in direction `dir` using Van Leer.
 *
 * Does nothing ifdef `SCALAR`. 
 *
 * Advects using the Van Leer flux limiter described in 
 * Anninos and Fragile (2002).
 */
void van_leer_E_dir(hydro_fields f, hydro_params p, int dir) {
#ifndef SCALAR
  int x, y, z;
  float r;
  float ***delta = make_field(p);
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

  // This loop constructs equivalent of (3-14) from Anninos and Fragile (2002).
  for(x = 1; x <= p.slicex; x++) {
    for(y = 1; y <= p.slicey; y++) {
      for(z = 0; z < p.Lz; z++) {
	r = ((f.E[x+dx][y+dy][(z+dz)%p.Lz] - f.E[x][y][z])
	     *(f.E[x][y][z]-f.E[x-dx][y-dy][(z-dz+p.Lz)%p.Lz]));
	if(r>0) {
	  delta[x][y][z] =  2*r/(p.dx*(f.E[x+dx][y+dy][(z+dz)%p.Lz]
				       - f.E[x-dx][y-dy][(z-dz+p.Lz)%p.Lz]));
	}
	else {
	  delta[x][y][z] = 0;
	}
      }
    }
  }

  halo_field(delta, p);

  // This loop constructs the terms in the brackets of (3-12) using (3-13)
  // from Anninos and Fragile (2002).
  for(x = 1; x <= p.slicex; x++) {
    for(y = 1; y <= p.slicey; y++) {
      for(z = 0; z < p.Lz; z++) {
	if(f.V[dir][x][y][z]>0){
	  f.F[dir][x][y][z] = (f.V[dir][x][y][z]
			       *(f.E[x-dx][y-dy][(z-dz+p.Lz)%p.Lz]
				 + 0.5*(p.dx - f.V[dir][x][y][z]*p.dt)
				 *delta[x-dx][y-dy][(z-dz+p.Lz)%p.Lz]));
	}
	else{
	  f.F[dir][x][y][z] = (f.V[dir][x][y][z]
			       *(f.E[x][y][z]
				 - 0.5*(p.dx + f.V[dir][x][y][z]*p.dt)
				 *delta[x][y][z]));
	}
      }
    }
  }

  halo_field(f.F[dir], p);

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

  free_field(p, delta);
#endif


  
}


/** Advection for `f.Z` in direction `dir` using Van Leer.
 *
 * Does nothing ifdef `SCALAR`. 
 *
 * Advects using the Van Leer flux limiter described in 
 * Anninos and Fragile (2002).
 * This advects each component of Z in the direction `dir`.
 */
void van_leer_Z_dir(hydro_fields f, hydro_params p, int dir) {

#ifndef SCALAR
  int x, y, z, i;

  float r;
  int dx = 0;
  int dy = 0;
  int dz = 0;

  if(dir == 0){
    dx = 1;
  }
  if(dir == 1){
    dy = 1;
  }
  if(dir == 2){
    dz = 1;
  }

  float ***Vface = make_field(p);
  float ****delta = make_vector(p);
  for(x = 1; x <= p.slicex; x++) {
    for(y = 1; y <= p.slicey; y++) {
      for(z = 0; z < p.Lz; z++) {
	for(i = 0; i < 3; i++){
	  r = ((f.Z[i][x+dx][y+dy][(z+dz)%p.Lz] - f.Z[i][x][y][z])
	       *(f.Z[i][x][y][z] - f.Z[i][x-dx][y-dy][(z-dz+p.Lz)%p.Lz]));

	  if(r>0){
	    delta[i][x][y][z] = 2*r/(p.dx*(f.Z[i][x+dx][y+dy][(z+dz)%p.Lz]
					   - f.Z[i][x-dx][y-dy][(z-dz+p.Lz)%p.Lz]));
	  }
	  else{
	    delta[i][x][y][z] = 0;
	  }

	}
      }
    }
  }
  
  halo_field(delta[0], p);
  halo_field(delta[1], p);
  halo_field(delta[2], p);

  
  for(x = 1; x <= p.slicex; x++) {
    for(y = 1; y <= p.slicey; y++) {
      for(z = 0; z < p.Lz; z++) {
	// Get corrected gradients (3-14) in Anninos and Fragile (2002)

	// We need to find the flux between nodes which requires us to find the
	// velocity on the appropriate 'face' between nodes.

	if(dir == 0){
	  Vface[x][y][z] =  0.125*((f.V[0][x][y][z] + f.V[0][x-1][y][z])
				   + (f.V[0][x][y-1][z] + f.V[0][x-1][y-1][z])
				   + (f.V[0][x][y][(z-1+p.Lz)%p.Lz]
				      + f.V[0][x-1][y][(z-1+p.Lz)%p.Lz])
				   + (f.V[0][x][y-1][(z-1+p.Lz)%p.Lz]
				      + f.V[0][x-1][y-1][(z-1+p.Lz)%p.Lz]));
	}
	else if(dir == 1){
	  Vface[x][y][z] = 0.125*((f.V[1][x][y][z] + f.V[1][x][y-1][z])
				  + (f.V[1][x-1][y][z] + f.V[1][x-1][y-1][z])
				  + (f.V[1][x][y][(z-1+p.Lz)%p.Lz]
				     + f.V[1][x][y-1][(z-1+p.Lz)%p.Lz])
				  + (f.V[1][x-1][y][(z-1+p.Lz)%p.Lz]
				     + f.V[1][x-1][y-1][(z-1+p.Lz)%p.Lz]));
	}
	else if(dir == 2){
	  Vface[x][y][z] = 0.125*((  f.V[2][x][y][z]
				     + f.V[2][x][y][(z-1+p.Lz)%p.Lz])
				  + (f.V[2][x-1][y][z]
				     + f.V[2][x-1][y][(z-1+p.Lz)%p.Lz])
				  + (f.V[2][x][y-1][z]
				     + f.V[2][x][y-1][(z-1+p.Lz)%p.Lz])
				  + (f.V[2][x-1][y-1][z]
				     + f.V[2][x-1][y-1][(z-1+p.Lz)%p.Lz]));
	}


	for(i = 0; i < 3; i++){
	  if(Vface[x][y][z] >= 0.0){
	    f.F[i][x][y][z] = (Vface[x][y][z]*
			       (f.Z[i][x-dx][y-dy][(z-dz+p.Lz)%p.Lz]
				+ 0.5*(p.dx - Vface[x][y][z]*p.dt)
				*delta[i][x-dx][y-dy][(z-dz+p.Lz)%p.Lz]));
	  }
	  else{
	    f.F[i][x][y][z] = (Vface[x][y][z]*
			       (f.Z[i][x][y][z]
				- 0.5*(p.dx + Vface[x][y][z]*p.dt)
				*delta[i][x][y][z]));
	  }
	}
      }
    }
  }

  halo_field(f.F[0],p);
  halo_field(f.F[1],p);
  halo_field(f.F[2],p);

  // Eq (3-12)
  for(x = 1; x <= p.slicex; x++) {
    for(y = 1; y <= p.slicey; y++) {
      for(z = 0; z < p.Lz; z++) {
	for(i = 0; i < 3; i++){
	  f.Z[i][x][y][z] = (f.Z[i][x][y][z]
			     - p.dt*(f.F[i][x+dx][y+dy][(z+dz)%p.Lz]
				     - f.F[i][x][y][z])/p.dx);
	}
      }
    }
  }
  halo_field(f.Z[0], p);
  halo_field(f.Z[1], p);
  halo_field(f.Z[2], p);

  
  free_field(p, Vface);

  free_vector(p, delta);

#endif
}


/** Advection for `E` in all directions.
 *
 * Currently the order in which each direction is advected is
 * fixed. This can be changed to be random by switching to the code in
 * the comment. Note that this is not what the Grant and Wilson
 * suggest, rather they state the update should be split into two half
 * steps, and the order should be reversed for each half.
 */
void advect_E(hydro_fields f, hydro_params p, int adv_order) {

  int directions[18] = {0,1,2,
			2,1,0,
			2,0,1,
			1,0,2,
			1,2,0,
			0,2,1};
  
#ifdef VANLEER
  van_leer_E_dir(f, p, directions[(adv_order % 6)*3]);    
  van_leer_E_dir(f, p, directions[(adv_order % 6)*3 + 1]);
  van_leer_E_dir(f, p, directions[(adv_order % 6)*3 + 2]);
#else
  donor_E_dir(f, p, directions[(adv_order % 6)*3]);    
  donor_E_dir(f, p, directions[(adv_order % 6)*3 + 1]);
  donor_E_dir(f, p, directions[(adv_order % 6)*3 + 2]);
#endif //VANLEER
}



/** Advection for `f.Z` in all directions.
 *
 * Currently the order in which component of `f.Z` is advected is
 * fixed. This can be changed to be random by switching to the code in
 * the comment. Note that this is not what the Grant and Wilson
 * suggest, rather they state the update should be split into two half
 * steps, and the order of advection should be reversed for each half.
 */
void advect_Z(hydro_fields f, hydro_params p, int adv_order) {

  int directions[18] = {0,1,2,
			2,1,0,
			2,0,1,
			1,0,2,
			1,2,0,
			0,2,1};
  
#ifdef VANLEER
  van_leer_Z_dir(f, p, directions[(adv_order % 6)*3]);
  van_leer_Z_dir(f, p, directions[(adv_order % 6)*3 + 1]);
  van_leer_Z_dir(f, p, directions[(adv_order % 6)*3 + 2]);
#else
  donor_Z_dir(f, p, directions[(adv_order % 6)*3]);    
  donor_Z_dir(f, p, directions[(adv_order % 6)*3 + 1]);
  donor_Z_dir(f, p, directions[(adv_order % 6)*3 + 2]);
#endif //VANLEER
}


void advect(hydro_fields f, hydro_params p) {
  p.dt = p.dt/2.0;
  advect_E(f, p, 0);
  advect_Z(f, p, 0);
  advect_E(f, p, 1);
  advect_Z(f, p, 1);
  p.dt = p.dt*2.0;
}















