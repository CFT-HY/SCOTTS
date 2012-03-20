/* energy.c
 *
 * Calculation of total energy and related quantities.
 */
#include "hydro.h"


// straight from fortran, slightly verbose way of doing it!
double total_energy(double *kappa, double *E, double *gb,
		    double *xe, double *xc, double *phi, double *pifull,
		    double *T, int **nb, hydro_params params) {

  int x;

  double vol;

  double Etot = 0.0;
  double restE = 0.0;
  double kinE = 0.0;
  double kinphi = 0.0;
  double grdphi = 0.0;

  for(x=0;x<params.N;x++) {
    vol = xe[x]*xe[x]*xe[x] - xe[nb[x][1]]*xe[nb[x][1]]*xe[nb[x][1]];

    // rest energy
    restE += E[x]/gb[x]*vol;

    // kinetic energy
    kinE += kappa[x]*E[x]/gb[x]*(gb[x]*gb[x]-1.0)*vol;

    // momentum squared (scalar field kinetic energy)
    kinphi += 0.5*pifull[x]*pifull[x]*vol;

    // gradient term
    grdphi += 0.125*((phi[nb[x][0]] - phi[nb[x][1]])/params.dx)
      *((phi[nb[x][0]] - phi[nb[x][1]])/params.dx)*vol;

  }

  // Cast here
  vol = (double)(xe[params.N-1]*xe[params.N-1]*xe[params.N-1]);

  Etot = (restE+kinE+kinphi+grdphi);

  return Etot/vol;



}
    
			    
