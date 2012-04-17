/* energy.c
 *
 * Calculation of total energy and related quantities.
 */
#include "hydro.h"


// straight from fortran, slightly verbose way of doing it!
double total_energy(hydro_fields f, int **nb, hydro_params p) {

  int x;

  double vol;

  double Etot = 0.0;
  double restE = 0.0;
  double kinE = 0.0;
  double kinphi = 0.0;
  double grdphi = 0.0;

  for(x = 0; x < p.N; x++) {
    vol = f.xe[x]*f.xe[x]*f.xe[x] 
      - f.xe[nb[x][1]]*f.xe[nb[x][1]]*f.xe[nb[x][1]];

    // rest energy
    restE += f.E[x]/f.gb[x]*vol;

    // kinetic energy
    kinE += f.kappa[x]*f.E[x]/f.gb[x]*(f.gb[x]*f.gb[x]-1.0)*vol;

    // momentum squared (scalar field kinetic energy)
    kinphi += 0.5*f.pifull[x]*f.pifull[x]*vol;

    // gradient term
    grdphi += 0.125*((f.phi[nb[x][0]] - f.phi[nb[x][1]])/p.dx)
      *((f.phi[nb[x][0]] - f.phi[nb[x][1]])/p.dx)*vol;

  }

  // Cast here
  vol = (double)(f.xe[p.N-1]*f.xe[p.N-1]*f.xe[p.N-1]);

  Etot = (restE+kinE+kinphi+grdphi);

  return Etot/vol;



}
    
			    
