/* initial.c
 *
 * Initial conditions for hydro evolution.
 */
#include "hydro.h"

/* void create_1D_bubble(...)
 *
 * Straight from fortran. Not very useful since we're spherical.
 */
void create_1D_bubble(int N, double xxwall, double *xe, double *xc,
		      double *phi, double *pifull,
		      double *T, double *E, double *Z,
		      double *v, double *gb, double dx, double a,
		      double alpha, double gamma, double lambda,
		      double Tconst, double T0) {
  int x;

  const double Tc = 1.0;



  fprintf(stderr, "T0 = %lf\n", T0);
  double cstre = 1.0;
  double Rtensp = 1.0;
  double rMboso = sqrt(gamma*(Tconst*Tconst - T0*T0));
  double deofV = alpha*Tconst;

  double lbar = 9.0/2.0 *lambda*rMboso*rMboso/(deofV*deofV);


  for(x=0; x<N; x++) {
    xe[x] = -1.0*xxwall + ((double)x)*dx;
    xc[x] = -1.0*xxwall + ((double)x-0.5)*dx;

    phi[x] = cstre*3.0*(rMboso*rMboso/(2.0*deofV))
      *psibar(rMboso*(xc[x]+xxwall+0.5*dx)/Rtensp,lbar);


    pifull[x] = 0.0;
    T[x] = Tconst;

    E[x] = 3*a*T[x]*T[x]*T[x]*T[x] + Vf(alpha,gamma,lambda,T[x],T0,phi[x])
      - T[x]*VTf(alpha,gamma,lambda,T[x],T0,phi[x]);

    Z[x] = 0.0;
    v[x] = 0.0;
    gb[x] = 1.0;

  }

}


/* double psibar(double x, double lbar)
 *
 * Scale-invariant 1D bubble.
 * Straight from fortran.
 */
double psibar(double x, double lbar) {

  double thx2 = tanh(x)*tanh(x);
  double reps = 1e-8;

  if(fabs(thx2) < reps)
    return (2.0/lbar)*(1.0-sqrt(1.0-lbar));
  else if(fabs(thx2-1) < reps)
    return 0.0;

  double fx = (1.0-lbar/thx2)/(1.0 - 1.0/thx2);


  return (2.0/fx)*(1-sqrt(1.0-fx));
}


/* void create_gaussian_bubble(...)
 *
 * Create critical spherical gaussian bubble.
 * Straight from fortran.
 */
void create_gaussian_bubble(int N, double xxwall, double *xe, double *xc,
			    double *phi, double *pifull,
			    double *T, double *E, double *Z,
			    double *v, double *gb, double dx, double a,
			    double alpha, double gamma, double lambda,
			    double Tconst, double T0) {
  

  double sigmlo = 2.0*sqrt(2.0)/81.0 
    * alpha*alpha*alpha /(lambda*lambda*sqrt(lambda));

  fprintf(stderr, \
	  "** Initial conditions magic:\n" \
	  "** sigmlo %g\n", sigmlo);
  
  double phimin =  ( alpha*Tconst 
		    + sqrt((alpha*Tconst)*(alpha*Tconst) - 
			   4*lambda*gamma*(Tconst*Tconst-T0*T0)) )/ (2*lambda);

  double cstrab = 1.0*phimin;

  double Rlapab = 2.0*sigmlo/(-1.0*Vf(alpha,gamma,lambda,Tconst,T0,phimin));

  double Rtenab = Rlapab;

  fprintf(stderr, \
	  "** phimin %g, cstrab %g, Rtenab %g\n", \
	  phimin, cstrab, Rtenab);

  int x;

  for(x=0; x<N; x++) {
    xe[x] = -1.0*xxwall + (x-0.0)*dx;
    xc[x] = -1.0*xxwall + (x-0.5)*dx;

    phi[x] = cstrab*exp(-1.0*(xc[x] + xxwall + 0.5*dx)
			*(xc[x] + xxwall + 0.5*dx)
			/2.0/(Rtenab*Rtenab) );

    pifull[x] = 0.0;

    T[x] = Tconst;

    E[x] = 3.0*a*T[x]*T[x]*T[x]*T[x] + Vf(alpha,gamma,lambda,T[x],T0,phi[x])
      - T[x]*VTf(alpha,gamma,lambda,T[x],T0,phi[x]);


    // For debugging purposes
    // fprintf(stderr,"xc = %lf fi = %lf E = %lf\n", xc[x], phi[x],E[x]);

    Z[x] = 0.0;
    v[x] = 0.0;

    gb[x] = 1.0;

  }
}
