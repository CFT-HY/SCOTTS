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
 *		    \frac{1}{2}\gamma T_0^2 \phi^2 
 *                   -\frac{1}{4}\lambda\phi^4 + V_0)}}
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
 *
 * Alternatively in the BAG model we find `f.T` by solving 
 * \f[
 * \epsilon = 3 a(\phi) T^4 + V_B(\phi) 
 * \f]
 * and so
 * \f[
 *    T=\left( \dfrac{ \epsilon - V_B(\phi)}{3 a(\phi)}\right)^{1/4}
 * \f]
 * Then the pressure `f.p` is found from
 * \f[
 * P = a(\phi) T^4 - V_B(\phi)
 * \f]
 */
#include "hydro.h"


/** Find temperature by 'solving' equation of state.
 *
 * Root finding algorithm from internal energy equation.  Note that if
 * `Tfix` goes negative then at that site `f.T` will become NaN's. If
 * this happens we write out a slice through one of the locations
 * where `Tfix` became negative and kill the run.
 *
 */
void find_Ta(hydro_fields f, hydro_params p, int step) {

#ifndef SCALAR

  int x, y, z;

  float Tfix;
  int TfixError = 0;
  int sliceError = -1;

  for(x=1; x<=p.slicex; x++) {
    for(y=1; y<=p.slicey; y++) {
      for(z=0; z<p.Lz; z++) {

#ifdef BAG
    // Second term in power is just the potential,
    // should we just create an array here?
	Tfix=((f.E[x][y][z]/f.W[x][y][z]
	       - ( 0.5*p.gamma*f.phi[x][y][z]*f.phi[x][y][z]
		   - (p.alpha*f.phi[x][y][z]*f.phi[x][y][z]
		      *f.phi[x][y][z]/3.0)
		   + (0.25*p.lambda*f.phi[x][y][z]*f.phi[x][y][z]
		      *f.phi[x][y][z]*f.phi[x][y][z])
		   - p.V0))
	      /(3.*(p.gdeg - p.V0*f.phi[x][y][z]*f.phi[x][y][z]
		    *(2*f.phi[x][y][z]/p.phi_0 - 3.)
		    /(p.phi_0*p.phi_0))));
	if(Tfix < 0){
	  fprintf(stderr,"Tfix -ve at (%d,%d,%d) \n"
		  "Tfix = %f \n W = %f \n phi = %f \n E=%f \n",
		  p.shiftx+x-1,p.shifty+y-1,z,Tfix,f.W[x][y][z],
		  f.phi[x][y][z],f.E[x][y][z]);
	  TfixError = 1;
	  sliceError = p.shiftx + x - 1;
	}
	
	f.T[x][y][z]=pow(Tfix , 0.25);
#else

	/*
	 * NaN's happen when Tfix goes negative,
	 * they then spread out from here.
	 */ 
	Tfix = 0.25*p.gamma*p.gamma*f.phi[x][y][z]*f.phi[x][y][z]
	  *f.phi[x][y][z]*f.phi[x][y][z]
	  - 12.0*p.gdeg*(0.25*p.lambda*f.phi[x][y][z]*f.phi[x][y][z]
		      *f.phi[x][y][z]*f.phi[x][y][z]
			 - p.V0
			 - 0.5*p.gamma*p.T0*p.T0*f.phi[x][y][z]*f.phi[x][y][z]
			 - f.E[x][y][z]/f.W[x][y][z]);

	if(Tfix < 0){
	  fprintf(stderr,"Tfix -ve at (%d,%d,%d) \n"
		  "Tfix = %f \n W = %f \n phi = %f \n E=%f \n",
		  p.shiftx+x-1,p.shifty+y-1,z,Tfix,f.W[x][y][z],
		  f.phi[x][y][z],f.E[x][y][z]);
	  TfixError = 1;
	  sliceError = p.shiftx + x - 1;
	}
	
	f.T[x][y][z] = sqrt((1.0/(6.0*p.gdeg))
			    * (0.5*p.gamma*f.phi[x][y][z]*f.phi[x][y][z]
			       + sqrt(Tfix)));
#endif // BAG
      }
    }
  }

  // Not needed... except for phi^2/T damping
  halo_field(f.T, p);

  TfixError = reduce_max_int(TfixError, p);
  sliceError = reduce_max_int(sliceError, p);
  
  if(TfixError == 1){

#ifdef SILO
    printf0(p,"Writing out a slice where Tfix -ve, then aborting run\n");
    write_silo_slice_step(f, p, step, sliceError);
    MPI_Barrier(MPI_COMM_WORLD);
#else 
	printf0(p,"Aborting run as Tfix -ve\n");
#endif //SILO
    die(10);
  }


#endif // SCALAR
}




/** By 'solving' the equation of state, determine first T and then
 * pressure p and the adiabatic parameter kappa. 
 * 
 * Calls find_Ta() to obtain the temperature.
 */
void eq_of_state(hydro_fields f, hydro_params p, int step) {

#ifndef SCALAR

  int x, y, z;

  /*
   * Getting some space on the heap takes time -- but then, not sure 
   * Vpot call is most efficient either (cannot fuse loops between find_Ta,
   * Vpot and this function. Will get rid of both simultaneously.
   */
  float ***Vnew = make_field(p);

  float tolE = 1e-3;

  find_Ta(f, p, step);



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
