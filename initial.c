/* initial.c
 *
 * Initial conditions for hydro evolution.
 */
#include "hydro.h"

/* void create_1D_bubble(...)
 *
 * Straight from fortran. Not very useful since we're spherical.
 */
/*
void create_1D_bubble(hydro_fields f, hydro_params p) {

  int x;

  const double Tc = 1.0;



  fprintf(stderr, "T0 = %lf\n", p.T0);
  double cstre = 1.0;
  double Rtensp = 1.0;
  double rMboso = sqrt(p.gamma*(p.Tconst*p.Tconst
			      - p.T0*p.T0));
  double deofV = p.alpha*p.Tconst;

  double lbar = 9.0/2.0 *p.lambda*rMboso*rMboso/(deofV*deofV);


  for(x=0; x<p.N; x++) {
    f.xe[x] = -1.0*p.xxwall + ((double)x)*p.dx;
    f.xc[x] = -1.0*p.xxwall + ((double)x-0.5)*p.dx;

    f.phi[x] = cstre*3.0*(rMboso*rMboso/(2.0*deofV))
      *psibar(rMboso*(f.xc[x]+p.xxwall+0.5*p.dx)/Rtensp,lbar);


    f.pifull[x] = 0.0;
    f.T[x] = p.Tconst;

    f.E[x] = 3.0*p.a*f.T[x]*f.T[x]*f.T[x]*f.T[x]
      + Vf(p, f.T[x], f.phi[x])
      - f.T[x]*VTf(p, f.T[x], f.phi[x]);

    f.Z[x] = 0.0;
    f.v[x] = 0.0;
    f.W[x] = 1.0;

  }

}
*/

/* double psibar(double x, double lbar)
 *
 * Scale-invariant 1D bubble.
 * Straight from fortran.
 */
/*
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
*/


/* void create_gaussian_bubble(...)
 *
 * Create critical spherical gaussian bubble.
 * Straight from fortran.
 */
/*
void create_gaussian_bubble(hydro_fields f, hydro_params p) {
  

  double sigmlo = 2.0*sqrt(2.0)/81.0*p.alpha*p.alpha*p.alpha
    /(p.lambda*p.lambda*sqrt(p.lambda));

  fprintf(stderr, \
	  "** Initial conditions magic:\n" \
	  "** sigmlo %g\n", sigmlo);
  
  double phimin =  ( p.alpha*p.Tconst 
		    + sqrt((p.alpha*p.Tconst)
			   *(p.alpha*p.Tconst) - 
			   4*p.lambda*p.gamma*(p.Tconst
					   *p.Tconst -
					   p.T0*p.T0)) )
    / (2*p.lambda);

  double cstrab = 1.0*phimin;

  double Rlapab = 2.0*sigmlo/(-1.0*Vf(p, p.Tconst, phimin));

  double Rtenab = Rlapab;

  fprintf(stderr, \
	  "** phimin %g, cstrab %g, Rtenab %g\n", \
	  phimin, cstrab, Rtenab);

  int x;

  for(x = 0; x < p.N; x++) {
    f.xe[x] = -1.0*p.xxwall + (x-0.0)*p.dx;
    f.xc[x] = -1.0*p.xxwall + (x-0.5)*p.dx;

    f.phi[x] = cstrab*exp(-1.0*(f.xc[x] + p.xxwall + 0.5*p.dx)
			*(f.xc[x] + p.xxwall + 0.5*p.dx)
			/2.0/(Rtenab*Rtenab) );

    f.pifull[x] = 0.0;

    f.T[x] = p.Tconst;

    f.E[x] = 3.0*p.a*f.T[x]*f.T[x]*f.T[x]*f.T[x]
      + Vf(p, f.T[x], f.phi[x])
      - f.T[x]*VTf(p, f.T[x], f.phi[x]);

    // For debugging purposes
    // fprintf(stderr,"xc = %lf fi = %lf E = %lf\n", xc[x], phi[x],E[x]);

    f.Z[x] = 0.0;
    f.v[x] = 0.0;

    f.W[x] = 1.0;

  }
}
*/





void initial_3D(hydro_fields f, hydro_params p) {
  

  double sigmlo = 2.0*sqrt(2.0)/81.0*p.alpha*p.alpha*p.alpha
    /(p.lambda*p.lambda*sqrt(p.lambda));

  //  srand48(time());

  fprintf(stderr, \
	  "** Initial conditions magic:\n" \
	  "** sigmlo %g\n", sigmlo);
  
  double phimin =  ( p.alpha*p.Tconst 
		    + sqrt((p.alpha*p.Tconst)
			   *(p.alpha*p.Tconst) - 
			   4*p.lambda*p.gamma*(p.Tconst
					   *p.Tconst -
					   p.T0*p.T0)) )
    / (2*p.lambda);

  double cstrab = 1.0*phimin;

  double Rlapab = 2.0*sigmlo/(-1.0*Vf(p, p.Tconst, phimin));

  double Rtenab = Rlapab;

  fprintf(stderr, \
	  "** phimin %g, cstrab %g, Rtenab %g\n", \
	  phimin, cstrab, Rtenab);

  double Er, El, rE;

  rE = 50.0;

  Er = 1.0/(1.0+rE);
  El = rE*Er;



  int x, y, z;

  for(x = 0; x < p.Lx; x++) {
    for(y = 0; y < p.Ly; y++) {
      for(z = 0; z < p.Lz; z++) {

	f.phi[iix(x,y,z,p)] = 0.0; // cstrab*(1.0 + 0.1*drand48());

	f.pifull[iix(x,y,z,p)] = 0.0;
	
	f.T[iix(x,y,z,p)] = 0.0; // p.Tconst;

	
	/*			
	if( ((x+y+z) < (p.Lx+p.Ly+p.Lz)/2))
	  f.E[iix(x,y,z,p)] = El;
	else
	  f.E[iix(x,y,z,p)] = Er;
	*/

	if( sqrt((x-p.Lx/2)*(x-p.Lx/2)+(y-p.Ly/2)*(y-p.Ly/2)) < 40)
	  f.E[iix(x,y,z,p)] = El;
	else
	  f.E[iix(x,y,z,p)] = Er;

	
	//	f.E[iix(x,y,z,p)] = 1.0 + drand48()*1.5;

		  // f.E[iix(x,y,z,p)] = drand48(); // 1.0;
	  //	else
	  //	  f.E[iix(x,y,z,p)] = 0.1;
	  /* 3.0*p.a*f.T[x]*f.T[x]*f.T[x]*f.T[x]
	  + Vf(p, f.T[x], f.phi[x])
	  - f.T[x]*VTf(p, f.T[x], f.phi[x]); */
	
	
    // For debugging purposes
    // fprintf(stderr,"xc = %lf fi = %lf E = %lf\n", xc[x], phi[x],E[x]);

	f.Zx[iix(x,y,z,p)] = 0.0;	
	f.Zy[iix(x,y,z,p)] = 0.0;
	f.Zz[iix(x,y,z,p)] = 0.0;
	
	
	f.W[iix(x,y,z,p)] = 1.0;
	
      }
    }
  }
}

    





/*
void create_shock_tube(hydro_fields f, hydro_params p) {
  
  double Er, El;

  fprintf(stderr, "** shock tube\n");

  int x;

  double rE = 50.0;

  Er = 1.0/(1.0+rE);
  El = rE*Er;

  for(x = 0; x < p.N; x++) {
    f.xe[x] = (x-0.0)*p.dx;
    f.xc[x] = (x-0.5)*p.dx;

    if(x<p.N/2)
      f.E[x] = El;
    else
      f.E[x] = Er;

    f.phi[x] = 0.0;
    f.pi[x] = 0.0;
    f.pifull[x] = 0.0;
    f.Z[x] = 0.0;
    f.v[x] = 0.0;
    f.W[x] = 1.0;

  }
}
*/

