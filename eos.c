/* eos.c
 * 
 * Calculations of the equation of state and T
 */
#include "hydro.h"


// Straight fortran port
void find_Ta(double *E, double *gb, double *phi, double a,
	     double alpha, double gamma, double lambda, double T0,
	     int N, double *T) {

  int x;

  double Tfix;

  for(x=0; x<N; x++) {
    Tfix = 0.25*gamma*gamma*phi[x]*phi[x]*phi[x]*phi[x]
      - 12.0*a*(0.25*lambda*phi[x]*phi[x]*phi[x]*phi[x]
		- 0.5*gamma*T0*T0*phi[x]*phi[x]
		- E[x]/gb[x]);

    T[x] = sqrt(1.0/6.0/a * (0.5*gamma *phi[x]*phi[x] + sqrt(Tfix)));
  }
}


// Straight fortran port
void eq_of_state(double a, int N, double *E, double *gb, double *phi,
		 double alpha, double gamma, double lambda,
		 double T0, double *T, double *p, double *kappa) {

  int x;

  /*
   * Getting some space on the heap takes time -- but then, not sure 
   * Vpot call is most efficient either (cannot fuse loops between find_Ta,
   * Vpot and this function. Will get rid of both simultaneously.
   */
  double *Vnew = (double *) malloc(N*sizeof(double));

  double tolE = 1e-3;

  find_Ta(E, gb, phi, a, alpha, gamma, lambda, T0, N, T);

  Vpot(N, alpha, gamma, lambda, T, T0, phi, Vnew);

  for(x=0; x<N; x++) {
    if(E[x] < tolE*gb[x]*3.0*a*T[x]*T[x]*T[x]*T[x]) {
      fprintf(stderr,"E getting dangerously small due to -ve V cont.\n");
      exit(100);
    }

    p[x] = a*T[x]*T[x]*T[x]*T[x] - Vnew[x];
    kappa[x] = 1.0 + gb[x]*p[x]/E[x];
  }
  
  free(Vnew);

}
