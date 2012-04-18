/* main.c
 * 
 * Execution entrypoint.
 */
#include "hydro.h"


void alloc_fields(hydro_fields *f, hydro_params p) {
  // location labels set up in initial condition functions
  f->xe = (double *) malloc(p.N*sizeof(double));
  f->xc = (double *) malloc(p.N*sizeof(double));

  // fields below also set up in initial condition functions
  f->phi = (double *) malloc(p.N*sizeof(double));
  f->pifull = (double *) malloc(p.N*sizeof(double));
  f->T = (double *) malloc(p.N*sizeof(double));
  f->E = (double *) malloc(p.N*sizeof(double));
  f->Z = (double *) malloc(p.N*sizeof(double));
  f->v = (double *) malloc(p.N*sizeof(double));
  f->gb = (double *) malloc(p.N*sizeof(double));
  
  // pi gets initialised in backstep
  f->pi = (double *) malloc(p.N*sizeof(double));

  // phiold initialised by evolve_field
  f->phiold = (double *) malloc(p.N*sizeof(double));

  // equation of state (obtained by eq_of_state first time
  // and used in hydro)
  f->kappa = (double *) malloc(p.N*sizeof(double));

  // pressure (obtained by eq_of_state first time
  // and used in hydro)
  f->p = (double *) malloc(p.N*sizeof(double));

  // (calloc considered harmful)
}


void zero_fields(hydro_fields f, hydro_params p) {

  int x;

  // We don't need to do this because create_gaussian_bubble
  // initialises everything, but it keeps valgrind quiet.
  for(x=0;x<p.N;x++) {
    f.xe[x] = 0.0000;
    f.xc[x] = 0.0000;

    f.phi[x] = 0.0000;
    f.pifull[x] = 0.0000;
    f.T[x] = 0.0000;
    f.E[x] = 0.0000;
    f.Z[x] = 0.0000;
    f.v[x] = 0.0000;
    f.gb[x] = 0.0000;

    f.pi[x] = 0.0000;
    f.phiold[x] = 0.0000;
    f.kappa[x] = 0.0000;
    f.p[x] = 0.0000;
  }

}

void free_fields(hydro_fields *f) {
  free(f->xe);
  free(f->xc);
  free(f->phi);
  free(f->pifull);
  free(f->phiold);
  free(f->T);
  free(f->E);
  free(f->Z);
  free(f->v);
  free(f->gb);
  free(f->kappa);
  free(f->p);
  free(f->pi);
}

int **init_nb(hydro_params p) {
  int x;

  int **nb;

  // set up nb lookup table
  nb = (int **) malloc(p.N*sizeof(int *));
  
  for(x=0; x<p.N; x++) {
    
    nb[x] = (int *)malloc(2*sizeof(int));

    nb[x][0] = x+1;
    nb[x][1] = x-1;
    
    // Reflective
    if(nb[x][1] == -1)
      nb[x][1] = 0;
    
    if(nb[x][0] == p.N)
      nb[x][0] = p.N-1;
    
    // These would be the periodic bc's...
    // nb[x][0] = (x+1)%N;
    // nb[x][1] = (x+N-1)%N;
    // ... which make no sense on a sphere
  }

  return nb;
}

void free_nb(int **nb, hydro_params p) {
  int x;

  for(x=0; x<p.N; x++) {
    free(nb[x]);
  }
  
  free(nb);
}


int main(int argc, char *argv[])
{
  fprintf(stderr,"1D hydro port\n");

  // Parse params from stdin
  hydro_params p = get_parameters();

  // Calculate terms in potential
  p.alpha = 1.0/sqrt(2.0*p.sigma*
			      pow(p.lcorr, 5.0)/3.0);

  p.gamma = (p.Lheat 
		      + 6.0*p.sigma/p.lcorr) 
    / (6.0*p.sigma*p.lcorr);

  p.lambda = 1.0/(3.0*p.sigma*pow(p.lcorr,3.0));

  // What did we find?
  fprintf(stderr,
	  "-- calculated potential terms\n" \
	  "-- alpha %g, gamma %g, lambda %g\n", \
	  p.alpha, p.gamma, p.lambda);

  // Make these user-modifiable eventually
  p.Tconst = 0.86;
  p.xxwall = -.5;

  // gdeg = 34.25
  p.a = 34.25*M_PI*M_PI/90.0;

  // Ugly but it's a straight conversion of the fortran
  p.T0 = sqrt(1.0-2.0*p.alpha*p.alpha
		       /9.0/p.gamma/p.lambda);


  // time iterate
  int step;
  double t = 0.0;

  // space iterate
  int x;


  hydro_fields f;

  // Just runs malloc on all the fields therein
  alloc_fields(&f, p);

  fprintf(stderr, "- Allocated fields.\n");

  // For safety, set everything to zero
  zero_fields(f, p);

  fprintf(stderr, "- Zeroed fields.\n");

  // initial conditions
  create_gaussian_bubble(f, p);

  int **nb;

  nb = init_nb(p);

  fprintf(stderr, "- Initialised neighbour lookup table.\n");

  fprintf(stderr, "Initial avg energy per site: %g\n", 
	  total_energy(f, nb, p));

  // Output headers
  // printf("Step\ttime\tenergy\twallps\n");

  // Output of some field at every time step
  FILE *phi_fh = fopen("field.dat","w");

  // Step back for leapfrog initial conds
  evolve_backstep(f, nb, p);


  // Record how many lattice sites in output file, to keep self contained
  fwrite(&p.N, sizeof(int), 1, phi_fh);


  for(step = 0; step < p.steps; step++) {

    if(step % p.interval == 0) {
      // Write time step, x coords and velocity field
      fwrite(&t, sizeof(double), 1, phi_fh);
      fwrite(f.xc, sizeof(double), p.N, phi_fh);
      fwrite(f.v, sizeof(double), p.N, phi_fh);
      fwrite(f.phi, sizeof(double), p.N, phi_fh);
    }

    // Do field step
    evolve_field(f, nb, p);

    // Advection of state variables
    donor_E(f, nb, p);
    donor_Z(f, nb, p);

    // Calculate EOS
    eq_of_state(f, p);

    // Do the hydro bits
    evolve_hydro(f, nb, p);

    // Solve for T
    find_Ta(f, p);

    t += p.dt;

    // Don't write to stdout too often, and don't calculate too much
    if(step % 100 == 0)
      printf("%04d\t%6lf\t%6lf\t%6lf\n",
	     step, t,
	     total_energy(f, nb, p), 
	     wallpos(f, p));
    
  } // main loop ends here

  fclose(phi_fh);

  // Clean up memory
  free_nb(nb, p);
  free_fields(&f);
  
  return 0;
}
