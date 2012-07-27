/* eos.c
 * 
 * Calculations of the equation of state and T
 */
#include "hydro.h"


// Straight fortran port
void find_Ta(hydro_fields f, hydro_params p) {

  int x, y, z;

  double Tfix;

  for(x=0; x<p.Lx; x++) {
    for(y=0; y<p.Ly; y++) {
      for(z=0; z<p.Lz; z++) {
	/*
	 * NaN's happen when Tfix goes negative,
	 * they then spread out from here.
	 */
	Tfix = 0.25*p.gamma*p.gamma*f.phi[x][y][z]*f.phi[x][y][z]*f.phi[x][y][z]*f.phi[x][y][z]
	  - 12.0*p.a*(0.25*p.lambda*f.phi[x][y][z]*f.phi[x][y][z]*f.phi[x][y][z]*f.phi[x][y][z]
		      - 0.5*p.gamma*p.T0*p.T0*f.phi[x][y][z]*f.phi[x][y][z]
		      - f.E[x][y][z]/f.W[x][y][z]);
	
	//    if(Tfix < 0)
	//      Tfix = 0.0;
	
	f.T[x][y][z] = sqrt((1.0/(6.0*p.a))
			    * (0.5*p.gamma*f.phi[x][y][z]*f.phi[x][y][z] + sqrt(Tfix)));
      }
    }
  }

  //  halo_field(f.T, p);
}


// Straight fortran port
void eq_of_state(hydro_fields f, hydro_params p) {

  int x, y, z;

  /*
   * Getting some space on the heap takes time -- but then, not sure 
   * Vpot call is most efficient either (cannot fuse loops between find_Ta,
   * Vpot and this function. Will get rid of both simultaneously.
   */
  double *Vnew = (double *) malloc(p.N*sizeof(double));

  double tolE = 1e-3;

  find_Ta(f, p);



  Vpot(p, f.T[0][0], f.phi[0][0], Vnew);


  for(x=0; x<p.Lx; x++) {
    for(y=0; y<p.Ly; y++) {
      for(z=0; z<p.Lz; z++) {

	if(f.E[x][y][z] 
	   < tolE*f.W[x][y][z]*3.0*p.a*f.T[x][y][z]*f.T[x][y][z]*f.T[x][y][z]*f.T[x][y][z]) {
	  fprintf(stderr,"E getting dangerously small due to -ve V cont.\n");
	  exit(100);
	}
	
	f.p[x][y][z] = p.a*f.T[x][y][z]*f.T[x][y][z]*f.T[x][y][z]*f.T[x][y][z] - Vnew[x*p.Ly*p.Lz + y*p.Lz + z];
	f.kappa[x][y][z] = 1.0 + f.W[x][y][z]*f.p[x][y][z]/f.E[x][y][z];
	
      }
    }
  }
  
  free(Vnew);

  halo_field(f.p, p);
  halo_field(f.kappa, p);
  // halo_field(f.T, p);

}
