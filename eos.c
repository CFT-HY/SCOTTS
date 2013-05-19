/* eos.c
 * 
 * Calculations of the equation of state and T
 */
#include "hydro.h"


/* find_Ta(hydro_fields f, hydro_params p)
 *
 * Find temperature by 'solving' equation of state.
 * Directly copied from the 1+1 fortran code.
 */
void find_Ta(hydro_fields f, hydro_params p) {

#ifndef SCALAR

  int x, y, z;

  double Tfix;

  for(x=1; x<=p.slicex; x++) {
    for(y=1; y<=p.slicey; y++) {
      for(z=0; z<p.Lz; z++) {
	/*
	 * NaN's happen when Tfix goes negative,
	 * they then spread out from here.
	 */

	Tfix = 0.25*p.gamma*p.gamma*f.phi[x][y][z]*f.phi[x][y][z]
	  *f.phi[x][y][z]*f.phi[x][y][z]
	  - 12.0*p.gdeg*(0.25*p.lambda*f.phi[x][y][z]*f.phi[x][y][z]
		      *f.phi[x][y][z]*f.phi[x][y][z]
		      - 0.5*p.gamma*p.T0*p.T0*f.phi[x][y][z]*f.phi[x][y][z]
		      - f.E[x][y][z]/f.W[x][y][z]);
	
	// Tfix not used (and ruins energy conservation anyway)
	//    if(Tfix < 0)
	//      Tfix = 0.0;
	
	f.T[x][y][z] = sqrt((1.0/(6.0*p.gdeg))
			    * (0.5*p.gamma*f.phi[x][y][z]*f.phi[x][y][z]
			       + sqrt(Tfix)));

      }
    }
  }

  // Not needed...
  //  halo_field(f.T, p);

#endif // SCALAR
}



/* void eq_of_state(hydro_fields f, hydro_params p)
 *
 * By 'solving' the equation of state, determine first T and then
 * pressure p and the adiabatic parameter kappa.
 */
void eq_of_state(hydro_fields f, hydro_params p) {

#ifndef SCALAR

  int x, y, z;

  /*
   * Getting some space on the heap takes time -- but then, not sure 
   * Vpot call is most efficient either (cannot fuse loops between find_Ta,
   * Vpot and this function. Will get rid of both simultaneously.
   */
  double ***Vnew = make_field(p);

  double tolE = 1e-3;

  find_Ta(f, p);



  Vpot(p, f.T[0][0], f.phi[0][0], Vnew[0][0]);


  for(x=1; x<=p.slicex; x++) {
    for(y=1; y<=p.slicey; y++) {
      for(z=0; z<p.Lz; z++) {

	if(f.E[x][y][z] 
	   < tolE*f.W[x][y][z]*3.0*p.gdeg*f.T[x][y][z]
	   *f.T[x][y][z]*f.T[x][y][z]*f.T[x][y][z]) {
	  fprintf(stderr,"(at %d,%d,%d) E getting dangerously small "
		  "due to -ve V cont.\n",
		  p.shiftx+x-1,p.shifty+y-1,z);
	  die(100);
	}

	// pressure P is radiative pressure less the potential
	f.p[x][y][z] = p.gdeg*f.T[x][y][z]*f.T[x][y][z]
	  *f.T[x][y][z]*f.T[x][y][z] - Vnew[x][y][z];
	f.kappa[x][y][z] = 1.0 + f.W[x][y][z]*f.p[x][y][z]/f.E[x][y][z];
	
      }
    }
  }
  
  //  free(Vnew);
  free_field(p, Vnew);

  halo_field(f.p, p);
  halo_field(f.kappa, p);
  // halo_field(f.T, p);

#endif // SCALAR
}
