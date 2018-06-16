/** @file potential.c
 *
 * The potential, and various derivatives thereof.
 *
 * This file contains methods for calculating the potential \f$ V(\phi,T) \f$
 * and the derivatives \f$ \dfrac{\partial V}{\partial T} \f$,
 * \f$ \dfrac{\partial^2 V}{\partial T^2} \f$
 * and \f$ \dfrac{\partial V}{\partial \phi} \f$.
 *
 * The functions Vf(), VTf(), VTTf(), Vdf() correspondingly compute
 * the potential, it's first derivative wrt \f$ T \f$, it's second
 * derivative wrt \f$ T \f$, and it's derivative wrt \f$ \phi \f$ for
 * single values of \f$ \phi \f$ and \f$ T \f$.
 *
 * Vpot() calculates and stores in a supplied array the potential at all sites.
 *
 * Likewise Vdpot() computes and stores in a supplied array the derivative of
 * the potential wrt \f$ \phi \f$ at all sites.
 */
#include "hydro.h"


/** Potential given \f$ \phi \f$ and \f$ T \f$.
 *
 * Returns the potential \f$ V(\phi,T) \f$ for a given \f$ T \f$ and
 * \f$ \phi \f$.
 *
 * \f[ V(\phi,T) = \frac{1}{2} \gamma(T^2 - T_0^2)\phi^2 -
 * \frac{1}{3}\alpha T \phi^3 + \frac{1}{4}\lambda\phi^4 + 
 * \frac{1}{4}\frac{\gamma^2 T_0^4}{\lambda} \text.  \f]
 */
Real Vf(hydro_params p, Real T, Real this_phi) {
#ifdef TINDEP
  return 0.5*p.gamma*(p.Tconst*p.Tconst - p.T0*p.T0)*this_phi*this_phi
    - p.alpha*p.Tconst*this_phi*this_phi*this_phi/3.0
    + 0.25*p.lambda*this_phi*this_phi*this_phi*this_phi
    + p.V0;
#else
  return  0.5*p.gamma*(T*T - p.T0*p.T0)*this_phi*this_phi
    - p.alpha*T*this_phi*this_phi*this_phi/3.0
    + 0.25*p.lambda*this_phi*this_phi*this_phi*this_phi
    + p.V0;
#endif // TINDEP
}


/** First derivative of potential wrt \f$ \phi \f$ for given \f$ T \f$
 * and \f$ \phi \f$.
 *
 * Returns the derivative of the potential \f$ \dfrac{\partial
 * V(\phi,T)}{\partial \phi} \f$ for a given \f$ T \f$ and \f$ \phi
 * \f$.
 *
 * \f[ \dfrac{\partial V(\phi,T)}{\partial \phi} = \gamma(T^2 -
 * T_0^2)\phi - \alpha T \phi^2 + \lambda\phi^3\text.  \f]
 */
Real Vdf(hydro_params p, Real T, Real this_phi) {
#ifdef TINDEP
  return p.gamma*(p.Tconst*p.Tconst - p.T0*p.T0)*this_phi
    - p.alpha*p.Tconst*this_phi*this_phi
    + p.lambda*this_phi*this_phi*this_phi;
#else
  return p.gamma*(T*T - p.T0*p.T0)*this_phi
    - p.alpha*T*this_phi*this_phi
    + p.lambda*this_phi*this_phi*this_phi;
#endif
}


/** First derivative of potential wrt \f$ T \f$ for given \f$ T \f$ and
 * \f$ \phi \f$.
 *
 * Returns the derivative of the potential \f$ \dfrac{\partial
 * V(\phi,T)}{\partial T} \f$ for a given \f$ T \f$ and \f$ \phi \f$.
 *
 * \f[ V(\phi,T) = \gamma T \phi^2 - \frac{1}{3}\alpha \phi^3 \text.
 * \f]
 */
Real VTf(hydro_params p, Real T, Real this_phi) {
#ifdef TINDEP
  return 0.;
#else
  return p.gamma*T*this_phi*this_phi 
    - p.alpha*this_phi*this_phi*this_phi/3.0;
#endif // TINDEP
}


/** Second derivative of potential wrt \f$ T \f$ for given \f$ T \f$
 * and \f$ \phi \f$.
 *
 * Returns the second derivative of the potential \f$
 * \dfrac{\partial^2 V(\phi,T)}{\partial T^2} \f$ for a given \f$ T
 * \f$ and \f$ \phi \f$.
 *
 * \f[ V(\phi,T) = \gamma \phi^2 \text.  \f]
 */
Real VTTf(hydro_params p, Real T, Real this_phi) {
#ifdef TINDEP
  return 0.;
#else
  return p.gamma*this_phi*this_phi;
#endif // TINDEP
}


/** Calculate potential for all sites and populate `Vprecalc`.
 * 
 * Except that doing it this way
 * does not encourage the compiler to fuse the loops.
 * 
 * Supply \f$ T \f$ and \f$ \phi \f$ at all sites 
 * using `T` and `phi` arrays.
 */
void Vpot(hydro_params p,
	  Real ***T,
	  Real ***phi,
	  Real ***Vprecalc) {

  int x, y, z;

  for(x = 1; x <= p.slicex; x++) {
    for(y = 1; y <= p.slicey; y++) {
      for(z = 0; z < p.Lz; z++) {
	Vprecalc[x][y][z] = Vf(p, T[x][y][z], phi[x][y][z]);
      }
    }
  }
}

/** Calculate first derivative of potential wrt \f$ \phi \f$ for all
 * sites and populate `Vprecalc`.
 * 
 * Except that doing it this way
 * does not encourage the compiler to fuse the loops.
 * 
 * Supply \f$ T \f$ and \f$ \phi \f$ at all sites 
 * using `T` and `phi` arrays.
 */
void Vdpot(hydro_params p,
	   Real ***T,
	   Real ***phi,
	   Real ***Vprecalc) {

  int x, y, z;

  for(x = 1; x <= p.slicex; x++) {
    for(y = 1; y <= p.slicey; y++) {
      for(z = 0; z < p.Lz; z++) {
	Vprecalc[x][y][z] = Vdf(p, T[x][y][z], phi[x][y][z]);
      }
    }
  }

}
