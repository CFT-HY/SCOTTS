/* energy.c
 *
 * Calculation of total energy and related quantities.
 */
#include "hydro.h"


double field_energy(hydro_fields f, hydro_params p) {

  int x, y, z;

  double vol;

  double Etot = 0.0;


  vol = p.dx*p.dx*p.dx;

  for(x = 1; x <= p.Lx; x++) {
    for(y = 1; y <= p.Ly; y++) {
      for(z = 0; z < p.Lz; z++) {

	
	Etot += 0.5*f.pifull[x][y][z]*f.pifull[x][y][z]*vol;
	
	Etot += 0.5*((f.phi[x][y][z] - f.phi[x-1][y][z])/p.dx)
	  *((f.phi[x][y][z] - f.phi[x-1][y][z])/p.dx)*vol;
	
	Etot += 0.5*((f.phi[x][y][z] - f.phi[x][y-1][z])/p.dx)
	  *((f.phi[x][y][z] - f.phi[x][y-1][z])/p.dx)*vol;
	
	Etot += 0.5*((f.phi[x][y][z] - f.phi[x][y][(z-1+p.Lz)%p.Lz])/p.dx)
	  *((f.phi[x][y][z] - f.phi[x][y][(z-1+p.Lz)%p.Lz])/p.dx)*vol;
	
	Etot += Vf(p, f.T[x][y][z], f.phi[x][y][z])*vol;
      }
    }
  }
  

  // Cast here
  //  vol = (double)(f.xe[p.N-1]);
  //  vol = ((double)(p.N))*p.dx*p.dx*p.dx;

  return Etot;


}

// straight from fortran, slightly verbose way of doing it!
double total_energy(hydro_fields f, hydro_params p) {

  int x, y, z;

  double vol;

  double Etot = 0.0;
  double restE = 0.0;
  double kinE = 0.0;
  double kinphi = 0.0;
  double grdphi = 0.0;

  vol = p.dx*p.dx*p.dx;

  for(x = 1; x <= p.Lx; x++) {
    for(y = 1; y <= p.Ly; y++) {
      for(z = 0; z < p.Lz; z++) {


	
	// rest energy
	restE += (f.E[x][y][z]/f.W[x][y][z])*vol;
	
	// kinetic energy
	kinE += f.kappa[x][y][z]*(f.E[x][y][z]/f.W[x][y][z])*(f.W[x][y][z]*f.W[x][y][z]-1.0)*vol;
	
	// momentum squared (scalar field kinetic energy)
	kinphi += 0.5*f.pifull[x][y][z]*f.pifull[x][y][z]*vol;
	
	// gradient term
	
	grdphi += 0.125*((f.phi[x+1][y][z] - f.phi[x-1][y][z])/p.dx)
	  *((f.phi[x+1][y][z] - f.phi[x-1][y][z])/p.dx)*vol;
	
	grdphi += 0.125*((f.phi[x][y+1][z] - f.phi[x][y-1][z])/p.dx)
	  *((f.phi[x][y+1][z] - f.phi[x][y-1][z])/p.dx)*vol;
	
	grdphi += 0.125*((f.phi[x][y][(z+1)%p.Lz] - f.phi[x][y][(z-1+p.Lz)%p.Lz])/p.dx)
	  *((f.phi[x][y][(z+1)%p.Lz] - f.phi[x][y][(z-1+p.Lz)%p.Lz])/p.dx)*vol;


      }
      
    }
  }

  // Cast here
  //  vol = (double)(f.xe[p.N-1]);
  //  vol = ((double)(p.Lx*p.Ly*p.Lz))*p.dx*p.dx*p.dx;

  //  fprintf(stderr,"%lf/%lf/%lf/%lf\n",restE,kinE,kinphi,grdphi);
  
  Etot = (restE+kinE+kinphi+grdphi);
  
  return Etot; // /vol;



}
    
			    
