/* eos.c
 * 
 * Calculations of the equation of state and T
 */
#include "hydro.h"


// Straight fortran port
void find_Ta(double *E, double *gb, double *phi, double *T,
	     hydro_params params) {

  int x;

  double Tfix;

  for(x=0; x<params.N; x++) {
    /*
     * NaN's happen when Tfix goes negative,
     * they then spread out from here.
     */
    Tfix = 0.25*params.gamma*params.gamma*phi[x]*phi[x]*phi[x]*phi[x]
      - 12.0*params.a*(0.25*params.lambda*phi[x]*phi[x]*phi[x]*phi[x]
		- 0.5*params.gamma*params.T0*params.T0*phi[x]*phi[x]
		- E[x]/gb[x]);

    //    if(Tfix < 0)
    //      Tfix = 0.0;

    T[x] = sqrt(1.0/6.0/params.a * (0.5*params.gamma 
				    *phi[x]*phi[x] + sqrt(Tfix)));
  }
}


// Straight fortran port
void eq_of_state(double *E, double *gb, double *phi,
		 double *T, double *p, double *kappa, hydro_params params) {

  int x;

  /*
   * Getting some space on the heap takes time -- but then, not sure 
   * Vpot call is most efficient either (cannot fuse loops between find_Ta,
   * Vpot and this function. Will get rid of both simultaneously.
   */
  double *Vnew = (double *) malloc(params.N*sizeof(double));

  double tolE = 1e-3;

  find_Ta(E, gb, phi, T, params);

  Vpot(params, T, phi, Vnew);

  for(x=0; x<params.N; x++) {
    if(E[x] < tolE*gb[x]*3.0*params.a*T[x]*T[x]*T[x]*T[x]) {
      fprintf(stderr,"E getting dangerously small due to -ve V cont.\n");
      exit(100);
    }

    p[x] = params.a*T[x]*T[x]*T[x]*T[x] - Vnew[x];
    kappa[x] = 1.0 + gb[x]*p[x]/E[x];
  }
  
  free(Vnew);

}
