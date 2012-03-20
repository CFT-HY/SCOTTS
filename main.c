/* main.c
 * 
 * Execution entrypoint.
 */
#include "hydro.h"


int main(int argc, char *argv[])
{
  fprintf(stderr,"1D hydro port\n");

  // Parse params from stdin
  hydro_params params = get_parameters();

  // Calculate terms in potential
  params.alpha = 1.0/sqrt(2.0*params.sigma*
			      pow(params.lcorr, 5.0)/3.0);

  params.gamma = (params.Lheat 
		      + 6.0*params.sigma/params.lcorr) 
    / (6.0*params.sigma*params.lcorr);

  params.lambda = 1.0/(3.0*params.sigma*pow(params.lcorr,3.0));

  // What did we find?
  fprintf(stderr,
	  "-- calculated potential terms\n" \
	  "-- alpha %g, gamma %g, lambda %g\n", \
	  params.alpha, params.gamma, params.lambda);

  // Make these user-modifiable eventually
  params.Tconst = 0.86;
  params.xxwall = -.5;

  // gdeg = 34.25
  params.a = 34.25*M_PI*M_PI/90.0;

  // Ugly but it's a straight conversion of the fortran
  params.T0 = sqrt(1.0-2.0*params.alpha*params.alpha
		       /9.0/params.gamma/params.lambda);


  // time iterate
  int step;
  double t = 0.0;

  // space iterate
  int x;

  // location labels set up in initial condition functions
  double *xe = (double *) malloc(params.N*sizeof(double));
  double *xc = (double *) malloc(params.N*sizeof(double));

  // fields below also set up in initial condition functions
  double *phi = (double *) malloc(params.N*sizeof(double));
  double *pifull = (double *) malloc(params.N*sizeof(double));
  double *T = (double *) malloc(params.N*sizeof(double));
  double *E = (double *) malloc(params.N*sizeof(double));
  double *Z = (double *) malloc(params.N*sizeof(double));
  double *v = (double *) malloc(params.N*sizeof(double));
  double *gb = (double *) malloc(params.N*sizeof(double));
  
  // pi gets initialised in backstep
  double *pi = (double *) malloc(params.N*sizeof(double));

  // phiold initialised by evolve_field
  double *phiold = (double *) malloc(params.N*sizeof(double));

  // equation of state (obtained by eq_of_state first time
  // and used in hydro)
  double *kappa = (double *) malloc(params.N*sizeof(double));

  // pressure (obtained by eq_of_state first time
  // and used in hydro)
  double *p = (double *) malloc(params.N*sizeof(double));

  // (calloc considered harmful)

  // We don't need to do this because create_gaussian_bubble
  // initialises everything, but it keeps valgrind quiet.
  for(x=0;x<params.N;x++) {
    xe[x] = 0.0000;
    xc[x] = 0.0000;

    phi[x] = 0.0000;
    pifull[x] = 0.0000;
    T[x] = 0.0000;
    E[x] = 0.0000;
    Z[x] = 0.0000;
    v[x] = 0.0000;
    gb[x] = 0.0000;

    pi[x] = 0.0000;
    phiold[x] = 0.0000;
    kappa[x] = 0.0000;
    p[x] = 0.0000;
  }



  // initial conditions
  create_gaussian_bubble(xe, xc, phi, pifull, T, E, Z, v, gb, params);


  // set up nb lookup table
  int **nb = (int **) malloc(params.N*sizeof(int *));

  for(x=0; x<params.N; x++) {

    nb[x] = (int *)malloc(2*sizeof(int));

    nb[x][0] = x+1;
    nb[x][1] = x-1;

    // Reflective
    if(nb[x][1] == -1)
      nb[x][1] = 0;

    if(nb[x][0] == params.N)
      nb[x][0] = params.N-1;

    // These would be the periodic bc's...
    // nb[x][0] = (x+1)%N;
    // nb[x][1] = (x+N-1)%N;
    // ... which make no sense on a sphere
  }


  fprintf(stderr, "initial avg energy per site: %g\n", 
	  total_energy(kappa, E, gb, xe, xc, phi, pifull, 
		       T, nb, params));

  // Output headers
  // printf("Step\ttime\tenergy\twallps\n");

  // Output of some field at every time step
  FILE *phi_fh = fopen("field.dat","w");

  // Step back for leapfrog initial conds
  evolve_backstep(phi, pifull, xe, xc,
		  T, pi, nb, params);


  // Record how many lattice sites in output file, to keep self contained
  fwrite(&params.N, sizeof(int), 1, phi_fh);


  for(step = 0; step < params.steps; step++) {

    // Write time step, x coords and velocity field
    // in principle don't need to do it every time step
    fwrite(&t, sizeof(double), 1, phi_fh);
    fwrite(xc, sizeof(double), params.N, phi_fh);
    fwrite(v, sizeof(double), params.N, phi_fh);
    fwrite(phi, sizeof(double), params.N, phi_fh);

    // Do field step
    evolve_field(gb, v, xe, xc, pi, phi,
		 T, phiold, pifull, nb, params);


    // Advection of state variables
    donor_r(v, xe, xc, E, nb, params);
    donor_z(v, xe, xc, Z, nb, params);


    // Calculate EOS
    eq_of_state(E, gb, phi,
		T, p, kappa, params);


    // Do the hydro bits
    evolve_hydro(kappa,
		 phiold, phi, pi, p,
		 xe, xc,
		 E, Z, v, gb,
		 T, nb, params);


    // Solve for T
    find_Ta(E, gb, phi, T, params);

    t += params.dt;

    // Don't write to stdout too often, and don't calculate too much
    if(step % 100 == 0)
      printf("%04d\t%6lf\t%6lf\t%6lf\n",
	     step, t,
	     total_energy(kappa, E, gb, xe, xc, phi, pifull, T, nb, params), 
	     wallpos(xc, phi, T, params));
    
  } // main loop ends here


  fclose(phi_fh);

  // Clean up memory
  free(xe);
  free(xc);
  free(phi);
  free(pifull);
  free(phiold);
  free(T);
  free(E);
  free(Z);
  free(v);
  free(gb);
  free(kappa);
  free(p);
  
  for(x=0; x<params.N; x++) {
    free(nb[x]);
  }
  
  free(nb);
  free(pi);
 
  
  return 0;
}
