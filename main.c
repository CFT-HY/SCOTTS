/* main.c
 * 
 * Execution entrypoint.
 */
#include "hydro.h"


int main(int argc, char *argv[])
{
  fprintf(stderr,"1D hydro port\n");

  // Parse parameters from stdin
  hydro_params parameters = get_parameters();

  // Calculate terms in potential
  parameters.alpha = 1.0/sqrt(2.0*parameters.sigma*
			      pow(parameters.lcorr, 5.0)/3.0);

  parameters.gamma = (parameters.Lheat 
		      + 6.0*parameters.sigma/parameters.lcorr) 
    / (6.0*parameters.sigma*parameters.lcorr);

  parameters.lambda = 1.0/(3.0*parameters.sigma*pow(parameters.lcorr,3.0));

  // What did we find?
  fprintf(stderr,
	  "-- calculated potential terms\n" \
	  "-- alpha %g, gamma %g, lambda %g\n", \
	  parameters.alpha, parameters.gamma, parameters.lambda);

  // Make these user-modifiable eventually
  parameters.Tconst = 0.86;
  parameters.xxwall = -.5;

  // gdeg = 34.25
  parameters.a = 34.25*M_PI*M_PI/90.0;

  // Ugly but it's a straight conversion of the fortran
  parameters.T0 = sqrt(1.0-2.0*parameters.alpha*parameters.alpha
		       /9.0/parameters.gamma/parameters.lambda);


  // time iterate
  int step;
  double t = 0.0;

  // space iterate
  int x;

  // location labels set up in initial condition functions
  double *xe = (double *) malloc(parameters.N*sizeof(double));
  double *xc = (double *) malloc(parameters.N*sizeof(double));

  // fields below also set up in initial condition functions
  double *phi = (double *) malloc(parameters.N*sizeof(double));
  double *pifull = (double *) malloc(parameters.N*sizeof(double));
  double *T = (double *) malloc(parameters.N*sizeof(double));
  double *E = (double *) malloc(parameters.N*sizeof(double));
  double *Z = (double *) malloc(parameters.N*sizeof(double));
  double *v = (double *) malloc(parameters.N*sizeof(double));
  double *gb = (double *) malloc(parameters.N*sizeof(double));
  
  // pi gets initialised in backstep
  double *pi = (double *) malloc(parameters.N*sizeof(double));

  // Q, don't actually care, only used in evolve_hydro
  // may be able to remove it
  double *Q = (double *) malloc(parameters.N*sizeof(double));

  // phiold initialised by evolve_field
  double *phiold = (double *) malloc(parameters.N*sizeof(double));

  // equation of state (obtained by eq_of_state first time
  // and used in hydro)
  double *kappa = (double *) malloc(parameters.N*sizeof(double));

  // pressure (obtained by eq_of_state first time
  // and used in hydro)
  double *p = (double *) malloc(parameters.N*sizeof(double));

  // (calloc considered harmful)

  // We don't need to do this because create_gaussian_bubble
  // initialises everything, but it keeps valgrind quiet.
  for(x=0;x<parameters.N;x++) {
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
    Q[x] = 0.0000;
    phiold[x] = 0.0000;
    kappa[x] = 0.0000;
    p[x] = 0.0000;
  }



  // initial conditions
  create_gaussian_bubble(parameters.N, parameters.xxwall,
			 xe, xc, phi, pifull, T, E, Z, v, gb, 
			 parameters.dx, parameters.a, parameters.alpha,
			 parameters.gamma, parameters.lambda,
			 parameters.Tconst, parameters.T0);

  // set up nb lookup table
  int **nb = (int **) malloc(parameters.N*sizeof(int *));

  for(x=0; x<parameters.N; x++) {

    nb[x] = (int *)malloc(2*sizeof(int));

    nb[x][0] = x+1;
    nb[x][1] = x-1;

    // Reflective
    if(nb[x][1] == -1)
      nb[x][1] = 0;

    if(nb[x][0] == parameters.N)
      nb[x][0] = parameters.N-1;

    // These would be the periodic bc's...
    // nb[x][0] = (x+1)%N;
    // nb[x][1] = (x+N-1)%N;
    // ... which make no sense on a sphere
  }


  fprintf(stderr, "initial avg energy per site: %g\n", 
	  total_energy(parameters.dx, parameters.N,
		       kappa, E, gb, xe, xc, phi, pifull, Q, 
		       parameters.alpha, parameters.gamma, parameters.lambda,
		       T, nb));

  // Output headers
  // printf("Step\ttime\tenergy\twallps\n");

  // Output of some field at every time step
  FILE *phi_fh = fopen("field.dat","w");

  // Step back for leapfrog initial conds
  evolve_backstep(parameters.dt, parameters.dx, parameters.N, 
		  phi, pifull, xe, xc,
		  parameters.alpha, parameters.gamma, parameters.lambda,
		  T, parameters.T0, pi, nb);


  // Record how many lattice sites in output file, to keep self contained
  fwrite(&parameters.N, sizeof(int), 1, phi_fh);


  for(step = 0; step < parameters.steps; step++) {

    // Write time step, x coords and velocity field
    // in principle don't need to do it every time step
    fwrite(&t, sizeof(double), 1, phi_fh);
    fwrite(xc, sizeof(double), parameters.N, phi_fh);
    fwrite(v, sizeof(double), parameters.N, phi_fh);
    fwrite(phi, sizeof(double), parameters.N, phi_fh);

    // Do field step
    evolve_field(parameters.dt, parameters.dx, parameters.C, parameters.N,
		 gb, v, xe, xc, pi, phi,
		 parameters.alpha, parameters.gamma, parameters.lambda,
		 T, parameters.T0, phiold, pifull, nb);


    // Advection of state variables
    donor_r(parameters.dt, parameters.dx, parameters.N, v, xe, xc, E, nb);
    donor_z(parameters.dt, parameters.dx, parameters.N, v, xe, xc, Z, nb);


    // Calculate EOS
    eq_of_state(parameters.a, parameters.N,
		E, gb, phi,
		 parameters.alpha, parameters.gamma, parameters.lambda,
		parameters.T0, T, p, kappa);


    // Do the hydro bits
    evolve_hydro(parameters.dt, parameters.dx, kappa,
		 parameters.C, parameters.Cav,
		 parameters.N, phiold, phi, pi, p,
		 xe, xc,
		 parameters.alpha, parameters.gamma, parameters.lambda,
		 parameters.T0, parameters.a, E, Z, v, gb,
		 T, Q, nb);


    // Solve for T
    find_Ta(E, gb, phi, parameters.a, 
	    parameters.alpha, parameters.gamma, parameters.lambda,
	    parameters.T0, parameters.N, T);

    t += parameters.dt;

    // Don't write to stdout too often, and don't calculate too much
    if(step % 100 == 0)
      printf("%04d\t%6lf\t%6lf\t%6lf\n",
	     step, t,
	     total_energy(parameters.dx, parameters.N,
			  kappa, E, gb, xe, xc, phi, pifull, Q,
			  parameters.alpha, parameters.gamma, 
			  parameters.lambda,
			  T, nb), 
	     wallpos(xc, parameters.N, phi, Q, parameters.dt,
		     parameters.alpha, parameters.gamma, 
		     parameters.lambda,
		     T, parameters.T0));
    
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
  
  for(x=0; x<parameters.N; x++) {
    free(nb[x]);
  }
  
  free(nb);
  free(pi);
  free(Q);
 
  
  return 0;
}
