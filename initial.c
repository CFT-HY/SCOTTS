/* initial.c
 *
 * Initial conditions for hydro evolution.
 */
#include "hydro.h"

/* void create_1D_bubble(...)
 *
 * Straight from fortran. Not very useful since we're spherical.
 */
void create_1D_bubble(double *xe, double *xc,
		      double *phi, double *pifull,
		      double *T, double *E, double *Z,
		      double *v, double *gb, hydro_params params) {

  int x;

  const double Tc = 1.0;



  fprintf(stderr, "T0 = %lf\n", params.T0);
  double cstre = 1.0;
  double Rtensp = 1.0;
  double rMboso = sqrt(params.gamma*(params.Tconst*params.Tconst
			      - params.T0*params.T0));
  double deofV = params.alpha*params.Tconst;

  double lbar = 9.0/2.0 *params.lambda*rMboso*rMboso/(deofV*deofV);


  for(x=0; x<params.N; x++) {
    xe[x] = -1.0*params.xxwall + ((double)x)*params.dx;
    xc[x] = -1.0*params.xxwall + ((double)x-0.5)*params.dx;

    phi[x] = cstre*3.0*(rMboso*rMboso/(2.0*deofV))
      *psibar(rMboso*(xc[x]+params.xxwall+0.5*params.dx)/Rtensp,lbar);


    pifull[x] = 0.0;
    T[x] = params.Tconst;

    E[x] = 3.0*params.a*T[x]*T[x]*T[x]*T[x]
      + Vf(params,T[x],phi[x])
      - T[x]*VTf(params,T[x],phi[x]);

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
void create_gaussian_bubble(double *xe, double *xc,
			    double *phi, double *pifull,
			    double *T, double *E, double *Z,
			    double *v, double *gb, 
			    hydro_params params) {
  

  double sigmlo = 2.0*sqrt(2.0)/81.0 
    * params.alpha*params.alpha*params.alpha
    /(params.lambda*params.lambda*sqrt(params.lambda));

  fprintf(stderr, \
	  "** Initial conditions magic:\n" \
	  "** sigmlo %g\n", sigmlo);
  
  double phimin =  ( params.alpha*params.Tconst 
		    + sqrt((params.alpha*params.Tconst)
			   *(params.alpha*params.Tconst) - 
			   4*params.lambda*params.gamma*(params.Tconst
					   *params.Tconst -
					   params.T0*params.T0)) )
    / (2*params.lambda);

  double cstrab = 1.0*phimin;

  double Rlapab = 2.0*sigmlo/(-1.0*Vf(params,params.Tconst,phimin));

  double Rtenab = Rlapab;

  fprintf(stderr, \
	  "** phimin %g, cstrab %g, Rtenab %g\n", \
	  phimin, cstrab, Rtenab);

  int x;

  for(x=0; x<params.N; x++) {
    xe[x] = -1.0*params.xxwall + (x-0.0)*params.dx;
    xc[x] = -1.0*params.xxwall + (x-0.5)*params.dx;

    phi[x] = cstrab*exp(-1.0*(xc[x] + params.xxwall + 0.5*params.dx)
			*(xc[x] + params.xxwall + 0.5*params.dx)
			/2.0/(Rtenab*Rtenab) );

    pifull[x] = 0.0;

    T[x] = params.Tconst;

    E[x] = 3.0*params.a*T[x]*T[x]*T[x]*T[x]
      + Vf(params,T[x],phi[x])
      - T[x]*VTf(params,T[x],phi[x]);


    // For debugging purposes
    // fprintf(stderr,"xc = %lf fi = %lf E = %lf\n", xc[x], phi[x],E[x]);

    Z[x] = 0.0;
    v[x] = 0.0;

    gb[x] = 1.0;

  }
}
