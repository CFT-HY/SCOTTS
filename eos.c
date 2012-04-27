/* eos.c
 * 
 * Calculations of the equation of state and T
 */
#include "hydro.h"


// Straight fortran port
void find_Ta(hydro_fields f, hydro_params p) {

  int x;

  double Tfix;

  for(x=0; x<p.N; x++) {
    /*
     * NaN's happen when Tfix goes negative,
     * they then spread out from here.
     */
    Tfix = 0.25*p.gamma*p.gamma*f.phi[x]*f.phi[x]*f.phi[x]*f.phi[x]
      - 12.0*p.a*(0.25*p.lambda*f.phi[x]*f.phi[x]*f.phi[x]*f.phi[x]
		- 0.5*p.gamma*p.T0*p.T0*f.phi[x]*f.phi[x]
		- f.E[x]/f.W[x]);

    //    if(Tfix < 0)
    //      Tfix = 0.0;

    f.T[x] = sqrt((1.0/(6.0*p.a)) * (0.5*p.gamma*f.phi[x]*f.phi[x] + sqrt(Tfix)));
  }
}


// Straight fortran port
void eq_of_state(hydro_fields f, hydro_params p) {

  int x;

  /*
   * Getting some space on the heap takes time -- but then, not sure 
   * Vpot call is most efficient either (cannot fuse loops between find_Ta,
   * Vpot and this function. Will get rid of both simultaneously.
   */
  double *Vnew = (double *) malloc(p.N*sizeof(double));

  double tolE = 1e-3;

  find_Ta(f, p);

  Vpot(p, f.T, f.phi, Vnew);

  for(x = 0; x < p.N; x++) {
    if(f.E[x] < tolE*f.W[x]*3.0*p.a*f.T[x]*f.T[x]*f.T[x]*f.T[x]) {
      fprintf(stderr,"E getting dangerously small due to -ve V cont.\n");
      exit(100);
    }

    f.p[x] = p.a*f.T[x]*f.T[x]*f.T[x]*f.T[x] - Vnew[x];
    f.kappa[x] = 1.0 + f.W[x]*f.p[x]/f.E[x];
  }
  
  free(Vnew);

}
