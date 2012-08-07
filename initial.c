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




void initial_scalar_bubble(hydro_fields f, hydro_params p) {
  double sigmlo = 2.0*sqrt(2.0)/81.0*p.alpha*p.alpha*p.alpha
    /(p.lambda*p.lambda*sqrt(p.lambda));

  if(!p.rank)
    fprintf(stderr,			   \
	    "** Initial conditions magic:\n"	\
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

  if(!p.rank)
    fprintf(stderr,					\
	    "** phimin %g, cstrab %g, Rtenab %g\n",	\
	    phimin, cstrab, Rtenab);

  int x, y, z;
  int i;

  //  fprintf(stderr,"my shiftx is %d and shifty is %d\n", p.shiftx, p.shifty);

  for(x=1;x<=p.slicex;x++) {
    for(y=1;y<=p.slicey;y++) {
      for(z=0;z<p.Lz;z++) {


    
	f.phi[x][y][z] = cstrab*exp(-1.0*
			  p.dx*p.dx*( ((double)(p.shiftx+x-1-p.Lx/2))*((double)(p.shiftx+x-1-p.Lx/2))
				      + ((double)(p.shifty+y-1-p.Ly/2))*((double)(p.shifty+y-1-p.Ly/2))
				      + ((double)(z-p.Lz/2))*((double)(z-p.Lz/2)))
			  /2.0/(Rtenab*Rtenab) );
	
    f.pifull[x][y][z] = 0.0;
	
    f.T[x][y][z] = p.Tconst;
	
    f.E[x][y][z] = 3.0*p.a*f.T[x][y][z]*f.T[x][y][z]
      *f.T[x][y][z]*f.T[x][y][z]
      + Vf(p, f.T[x][y][z], f.phi[x][y][z])
      - f.T[x][y][z]*VTf(p, f.T[x][y][z], f.phi[x][y][z]);
 
    f.Z[0][x][y][z] = 0.0;
    f.V[0][x][y][z] = 0.0;
    f.W[x][y][z] = 1.0;
      }
    }
  }
}





void initial_blank(hydro_fields f, hydro_params p) {
  double sigmlo = 2.0*sqrt(2.0)/81.0*p.alpha*p.alpha*p.alpha
    /(p.lambda*p.lambda*sqrt(p.lambda));

  if(!p.rank)
    fprintf(stderr,			   \
	    "** Initial conditions magic:\n"	\
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

  if(!p.rank)
    fprintf(stderr,					\
	    "** phimin %g, cstrab %g, Rtenab %g\n",	\
	    phimin, cstrab, Rtenab);

  int x, y, z;
  int i;

  //  fprintf(stderr,"my shiftx is %d and shifty is %d\n", p.shiftx, p.shifty);

  for(x=1;x<=p.slicex;x++) {
    for(y=1;y<=p.slicey;y++) {
      for(z=0;z<p.Lz;z++) {


    
	f.phi[x][y][z] = 0.0; 

	/*
cstrab*exp(-1.0*
			  p.dx*p.dx*( ((double)(p.shiftx+x-1-p.Lx/2))*((double)(p.shiftx+x-1-p.Lx/2))
				      + ((double)(p.shifty+y-1-p.Ly/2))*((double)(p.shifty+y-1-p.Ly/2))
				      + ((double)(z-p.Lz/2))*((double)(z-p.Lz/2)))
			  /2.0/(Rtenab*Rtenab) );
	*/
	
    f.pifull[x][y][z] = 0.0;
	
    f.T[x][y][z] = p.Tconst;
	
    f.E[x][y][z] = 3.0*p.a*f.T[x][y][z]*f.T[x][y][z]
      *f.T[x][y][z]*f.T[x][y][z]
      + Vf(p, f.T[x][y][z], f.phi[x][y][z])
      - f.T[x][y][z]*VTf(p, f.T[x][y][z], f.phi[x][y][z]);
 
    f.Z[0][x][y][z] = 0.0;
    f.V[0][x][y][z] = 0.0;
    f.W[x][y][z] = 1.0;
      }
    }
  }
}


int safe_distance(hydro_fields f, hydro_params p) {


  double sigmlo = 2.0*sqrt(2.0)/81.0*p.alpha*p.alpha*p.alpha
    /(p.lambda*p.lambda*sqrt(p.lambda));



  double phimin =  ( p.alpha*p.Tconst
		     + sqrt((p.alpha*p.Tconst)
			    *(p.alpha*p.Tconst) -
			    4*p.lambda*p.gamma*(p.Tconst
                                           *p.Tconst -
						p.T0*p.T0)) )    / (2*p.lambda);

  double cstrab = 1.0*phimin;

  double Rlapab = 2.0*sigmlo/(-1.0*Vf(p, p.Tconst, phimin));

  double Rtenab = Rlapab;


  /*
  if(!p.rank)
    fprintf(stderr,					\
	    "** phimin %g, cstrab %g, Rtenab %g\n",	\
	    phimin, cstrab, Rtenab);
  */

  double sigma = 1.0/(sqrt(p.dx*p.dx/2.0/(Rtenab*Rtenab) ));


  return (int)(round(sigma));
 	
}




int can_nucleate(hydro_fields f, hydro_params p, int x0, int y0, int z0) {

  int safedist = safe_distance(f,p);
  double threshphi = 0.0001;

  int shortx, shorty, shortz;
  int longx, longy, longz;
  int dx, dy, dz;

  int is_good = 1;
  int x, y, z;

  for(x=1;x<=p.slicex;x++) {
    for(y=1;y<=p.slicey;y++) {
      for(z=0;z<p.Lz;z++) {

	shortx = abs(p.shiftx+x-1-x0);
	shorty = abs(p.shifty+y-1-y0);
	shortz = abs(z - z0);

	longx = (p.Lx - shortx);
	longy = (p.Ly - shorty);
	longz = (p.Lz - shortz);

	if(shortx < longx)
	  dx = shortx;
	else
	  dx = longx;

	if(shorty < longy)
	  dy = shorty;
	else
	  dy = longy;

	if(shortz < longz)
	  dz = shortz;
	else
	  dz = longz;

	if(((dx*dx + dy*dy + dz*dz) < safedist*safedist) && (fabs(f.phi[x][y][z]) > threshphi)) {
	  fprintf(stderr,"Dead at (%d,%d,%d) because dx=%d, dy=%d, dz=%d, safedist=%d and phi %lf\n",
		  x0,y0,z0,dx,dy,dz,safedist,f.phi[x][y][z]);
	  is_good = 0;
	  break;
	}
      }

      if(!is_good)
	break;
    }

      if(!is_good)
	break;
  }

  return reduce_and(is_good, p);
}


void nucleate_at(hydro_fields f, hydro_params p, int x0, int y0, int z0) {


  double sigmlo = 2.0*sqrt(2.0)/81.0*p.alpha*p.alpha*p.alpha
    /(p.lambda*p.lambda*sqrt(p.lambda));



  double phimin =  ( p.alpha*p.Tconst
		     + sqrt((p.alpha*p.Tconst)
			    *(p.alpha*p.Tconst) -
			    4*p.lambda*p.gamma*(p.Tconst
                                           *p.Tconst -
						p.T0*p.T0)) )    / (2*p.lambda);

  double cstrab = 1.0*phimin;

  double Rlapab = 2.0*sigmlo/(-1.0*Vf(p, p.Tconst, phimin));

  double Rtenab = Rlapab;


  /*
  if(!p.rank)
    fprintf(stderr,					\
	    "** phimin %g, cstrab %g, Rtenab %g\n",	\
	    phimin, cstrab, Rtenab);
  */

  int x, y, z;
  int i;

  if(!p.rank)
    fprintf(stderr,"Nucleating at (%d,%d,%d)\n",x0,y0,z0);

  //  fprintf(stderr,"my shiftx is %d and shifty is %d\n", p.shiftx, p.shifty);

  int shortx, shorty, shortz;
  int longx, longy, longz;
  int dx, dy, dz;

  for(x=1;x<=p.slicex;x++) {
    for(y=1;y<=p.slicey;y++) {
      for(z=0;z<p.Lz;z++) {

	shortx = abs(p.shiftx+x-1-x0);
	shorty = abs(p.shifty+y-1-y0);
	shortz = abs(z - z0);

	longx = (p.Lx - shortx);
	longy = (p.Ly - shorty);
	longz = (p.Lz - shortz);

	if(shortx < longx)
	  dx = shortx;
	else
	  dx = longx;

	if(shorty < longy)
	  dy = shorty;
	else
	  dy = longy;

	if(shortz < longz)
	  dz = shortz;
	else
	  dz = longz;
    
	f.phi[x][y][z] += cstrab*exp(-1.0*
				     p.dx*p.dx*( (double) ((dx*dx) + (double)(dy*dy) + (double)(dz*dz)))
			  /2.0/(Rtenab*Rtenab) );
 	
      }
    }
  }

  //  fprintf(stderr,"Done\n");
}





void initial_3D(hydro_fields f, hydro_params p) {
  

  double sigmlo = 2.0*sqrt(2.0)/81.0*p.alpha*p.alpha*p.alpha
    /(p.lambda*p.lambda*sqrt(p.lambda));

  //  srand48(time());

  if(!p.rank)
    fprintf(stderr,			   \
	    "** Initial conditions magic:\n"	\
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

  if(!p.rank)
    fprintf(stderr,				  \
	    "** phimin %g, cstrab %g, Rtenab %g\n",	\
	    phimin, cstrab, Rtenab);

  double Er, El, rE;

  rE = 50.0;

  Er = 1.0/(1.0+rE);
  El = rE*Er;



  int x, y, z;
  int i;


  for(x=1;x<=p.slicex;x++) {
    for(y=1;y<=p.slicey;y++) {
      for(z=0;z<p.Lz;z++) {

    f.phi[x][y][z] = 0.0; // cstrab*(1.0 + 0.1*drand48());

    f.pifull[x][y][z] = 0.0;
	
    f.T[x][y][z] = 0.0; // p.Tconst;

    if( (x+p.shiftx-1 + y+p.shifty-1) < p.Lx/2  || (x+p.shiftx-1 + y+p.shifty-1) > 3*p.Lx/2) //  sqrt((x-p.Lx/2)*(x-p.Lx/2)+(y-p.Ly/2)*(y-p.Ly/2)) < 40)
      f.E[x][y][z] = El;
    else
      f.E[x][y][z] = Er;
	
	
    // For debugging purposes
    // fprintf(stderr,"xc = %lf fi = %lf E = %lf\n", xc[x], phi[x],E[x]);

    f.Z[0][x][y][z] = 0.0;	
    f.Z[1][x][y][z] = 0.0;
    f.Z[2][x][y][z] = 0.0;
    
    
    f.W[x][y][z] = 1.0;
      }
    }
  }
}

/*



void initial_step(hydro_fields f, hydro_params p, int **inverse) {
  

  double sigmlo = 2.0*sqrt(2.0)/81.0*p.alpha*p.alpha*p.alpha
    /(p.lambda*p.lambda*sqrt(p.lambda));

  //  srand48(time());

  if(!p.rank)
    fprintf(stderr,			   \
	    "** Initial conditions magic:\n"	\
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

  if(!p.rank)
    fprintf(stderr,				  \
	    "** phimin %g, cstrab %g, Rtenab %g\n",	\
	    phimin, cstrab, Rtenab);
  
  double Er, El, rE;

  rE = 50.0;

  Er = 1.0/(1.0+rE);
  El = rE*Er;



  int x, y, z;
  int i;

  for(i=0; i < p.N; i++) {

    x = inverse[i][0];
    y = inverse[i][1];
    z = inverse[i][2];


    f.phi[i] = 0.0; // cstrab*(1.0 + 0.1*drand48());
    
    f.pifull[i] = 0.0;
	
    f.T[i] = 0.0; // p.Tconst;
    
    if( x < p.Lx/2)
      f.E[i] = El;
    else
      f.E[i] = Er;
	
	
    // For debugging purposes
    // fprintf(stderr,"xc = %lf fi = %lf E = %lf\n", xc[x], phi[x],E[x]);

    f.Z[0][i] = 0.0;	
    f.Z[1][i] = 0.0;
    f.Z[2][i] = 0.0;
    
    
    f.W[i] = 1.0;
	
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

