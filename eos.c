/** @file eos.c
 * 
 * Calculations of the equation of state and T.
 * 
 * Note that find_Ta() can be called seperately from eq_of_state() but
 * eq_of_state() at the start. 
 *
 * We find `f.T` by solving
 * \f[
 * \epsilon = 3 a T^4 + V(\phi, T) - T \frac{\partial V}{\partial T}
 * \f]
 * for \f$ T \f$. 
 *
 * This leads to
 * \f[
 * T^2= \dfrac{ \frac{1}{2} \gamma \phi^2 
 *              \pm \sqrt{\frac{1}{4} \gamma^2\phi^4 +12a(\epsilon + 
 *		    \frac{1}{2}\gamma T_0^2 \phi^2 -\frac{1}{4}\lambda\phi^4)}}
 *	{6a}\text.
 * \f]
 *
 * Then `f.p` is found from
 * \f[
 * P = a T^4 + V(\phi, T)\text.
 * \f]
 *
 * The adiabatic constant `f.kappa` is given by
 * \f[
 * \kappa= 1 + \frac{PW}{E}\text.
 * \f]
 */
#include "hydro.h"


/** Find temperature by 'solving' equation of state.
 *
 * Root finding algorithm from internal energy equation. 
 * Note that if `Tfix` goes negative then at that site `f.T` will become
 * NaN's which then spread.
 *
 * Directly copied from the 1+1 fortran code.
 */
void find_Ta(hydro_fields f, hydro_params p) {

#ifndef SCALAR

  int x, y, z;

  float Tfix;

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

	if(Tfix < 0){
	  fprintf(stderr,"Tfix -ve at (%d,%d,%d) \n"
		  "Tfix = %f \n W = %f \n phi = %f \n E=%f \n",
		  p.shiftx+x-1,p.shifty+y-1,z,Tfix,f.W[x][y][z],
		  f.phi[x][y][z],f.E[x][y][z]);
	}
	
	f.T[x][y][z] = sqrt((1.0/(6.0*p.gdeg))
			    * (0.5*p.gamma*f.phi[x][y][z]*f.phi[x][y][z]
			       + sqrt(Tfix)));

      }
    }
  }

  // Not needed... except for phi^2/T damping
  halo_field(f.T, p);

#endif // SCALAR
}




/** By 'solving' the equation of state, determine first T and then
 * pressure p and the adiabatic parameter kappa. 
 * 
 * Calls find_Ta() to obtain the temperature.
 */
void eq_of_state(hydro_fields f, hydro_params p) {

#ifndef SCALAR

  int x, y, z;

  /*
   * Getting some space on the heap takes time -- but then, not sure 
   * Vpot call is most efficient either (cannot fuse loops between find_Ta,
   * Vpot and this function. Will get rid of both simultaneously.
   */
  float ***Vnew = make_field(p);

  float tolE = 1e-3;

  find_Ta(f, p);



  Vpot(p, f.T, f.phi, Vnew);


  for(x=1; x<=p.slicex; x++) {
    for(y=1; y<=p.slicey; y++) {
      for(z=0; z<p.Lz; z++) {
	/*
	if(f.E[x][y][z] 
	   < tolE*f.W[x][y][z]*3.0*p.gdeg*f.T[x][y][z]
	   *f.T[x][y][z]*f.T[x][y][z]*f.T[x][y][z]) {
	  fprintf(stderr,"(at %d,%d,%d) E getting dangerously small "
		  "due to -ve V cont.\n",
		  p.shiftx+x-1,p.shifty+y-1,z);
	  die(100);
	}
	*/
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
