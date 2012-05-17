/* energy.c
 *
 * Calculation of total energy and related quantities.
 */
#include "hydro.h"


double field_energy(hydro_fields f, int **nb, hydro_params p) {

  int x;

  double vol;

  double Etot = 0.0;

  for(x = 0; x < p.N; x++) {
    vol = p.dx*p.dx*p.dx;

    Etot += 0.5*f.pifull[x]*f.pifull[x]*vol;

    Etot += 0.5*((f.phi[x] - f.phi[nb[x][1]])/p.dx)
      *((f.phi[x] - f.phi[nb[x][1]])/p.dx)*vol;

    Etot += 0.5*((f.phi[x] - f.phi[nb[x][3]])/p.dx)
      *((f.phi[x] - f.phi[nb[x][3]])/p.dx)*vol;

    Etot += 0.5*((f.phi[x] - f.phi[nb[x][5]])/p.dx)
      *((f.phi[x] - f.phi[nb[x][5]])/p.dx)*vol;

    Etot += Vf(p, f.T[x], f.phi[x])*vol;

  }

  // Cast here
  //  vol = (double)(f.xe[p.N-1]);
  vol = ((double)(p.Lx*p.Ly*p.Lz))*p.dx*p.dx*p.dx;

  return Etot/vol;


}

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
    vol = p.dx*p.dx*p.dx;

    // rest energy
    restE += (f.E[x]/f.W[x])*vol;

    // kinetic energy
    kinE += f.kappa[x]*(f.E[x]/f.W[x])*(f.W[x]*f.W[x]-1.0)*vol;

    // momentum squared (scalar field kinetic energy)
    kinphi += 0.5*f.pifull[x]*f.pifull[x]*vol;

    // gradient term
    grdphi += 0.125*((f.phi[nb[x][0]] - f.phi[nb[x][1]])/p.dx)
      *((f.phi[nb[x][0]] - f.phi[nb[x][1]])/p.dx)*vol;

    grdphi += 0.125*((f.phi[nb[x][2]] - f.phi[nb[x][3]])/p.dx)
      *((f.phi[nb[x][2]] - f.phi[nb[x][3]])/p.dx)*vol;

    grdphi += 0.125*((f.phi[nb[x][4]] - f.phi[nb[x][5]])/p.dx)
      *((f.phi[nb[x][4]] - f.phi[nb[x][5]])/p.dx)*vol;

  }

  // Cast here
  //  vol = (double)(f.xe[p.N-1]);
  vol = ((double)(p.Lx*p.Ly*p.Lz))*p.dx*p.dx*p.dx;

  Etot = (restE+kinE+kinphi+grdphi);

  return Etot/vol;



}
    
			    
