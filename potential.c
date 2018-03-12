/* potential.c
 *
 * The potential, and various derivatives thereof.
 */
#include "hydro.h"


/* float Vf(hydro_params p, float T, float this_phi)
 *
 * Potential.
 */
float Vf(hydro_params p, float T, float this_phi) {
  return  0.5*p.gamma*(T*T - p.T0*p.T0)*this_phi*this_phi
    - p.alpha*T*this_phi*this_phi*this_phi/3.0
    + 0.25*p.lambda*this_phi*this_phi*this_phi*this_phi;
}


/* float Vdf(hydro_params p, float T, float this_phi)
 *
 * First derivative wrt phi
 */
float Vdf(hydro_params p, float T, float this_phi) {
  return p.gamma*(T*T - p.T0*p.T0)*this_phi
    - p.alpha*T*this_phi*this_phi
    + p.lambda*this_phi*this_phi*this_phi;
}


/* float VTf(hydro_params p, float T, float this_phi)
 *
 * First derivative wrt T
 */
float VTf(hydro_params p, float T, float this_phi) {
  return p.gamma*T*this_phi*this_phi 
    - p.alpha*this_phi*this_phi*this_phi/3.0;
}


/* float VTTf(hydro_params p, float T, float this_phi)
 *
 * Second derivative wrt T
 */
float VTTf(hydro_params p, float T, float this_phi) {
  return p.gamma*this_phi*this_phi;
}


/* void Vpot(hydro_params p, float *T, float *phi, float *Vprecalc)
 *
 * Calculate potential for all sites. Except that doing it this way
 * does not encourage the compiler to fuse the loops.
 */
void Vpot(hydro_params p,
	  float ***T,
	  float ***phi,
	  float ***Vprecalc) {

  int x, y, z;

  for(x = 1; x <= p.slicex; x++) {
    for(y = 1; y <= p.slicey; y++) {
      for(z = 0; z < p.Lz; z++) {
	Vprecalc[x][y][z] = Vf(p, T[x][y][z], phi[x][y][z]);
      }
    }
  }
}


/* void Vdpot(hydro_params p, float *T, float *phi, float *Vprecalc)
 *
 * Calculate first deriviative of potential wrt phi for all sites.
 */
void Vdpot(hydro_params p,
	   float ***T,
	   float ***phi,
	   float ***Vprecalc) {

  int x, y, z;

  for(x = 1; x <= p.slicex; x++) {
    for(y = 1; y <= p.slicey; y++) {
      for(z = 0; z < p.Lz; z++) {
	Vprecalc[x][y][z] = Vdf(p, T[x][y][z], phi[x][y][z]);
      }
    }
  }

}
