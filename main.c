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
  f->Zx = (double *) malloc(p.N*sizeof(double));
  f->Zy = (double *) malloc(p.N*sizeof(double));
  f->Zz = (double *) malloc(p.N*sizeof(double));
  f->Ux = (double *) malloc(p.N*sizeof(double));
  f->Uy = (double *) malloc(p.N*sizeof(double));
  f->Uz = (double *) malloc(p.N*sizeof(double));
  f->Vx = (double *) malloc(p.N*sizeof(double));
  f->Vy = (double *) malloc(p.N*sizeof(double));
  f->Vz = (double *) malloc(p.N*sizeof(double));
  f->W = (double *) malloc(p.N*sizeof(double));
  
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

  f->V = (double **) malloc(3*sizeof(double *));
  f->V[0] = f->Vx;
  f->V[1] = f->Vy;
  f->V[2] = f->Vz;

  f->U = (double **) malloc(3*sizeof(double *));

  f->U[0] = f->Ux;
  f->U[1] = f->Uy;
  f->U[2] = f->Uz;

  f->deltaM = (double **) malloc(3*sizeof(double *));

  f->deltaM[0] = (double *) malloc(p.N*sizeof(double));
  f->deltaM[1] = (double *) malloc(p.N*sizeof(double));
  f->deltaM[2] = (double *) malloc(p.N*sizeof(double));

  // (calloc considered harmful)
}


void dump(double *field, hydro_params p) {
  int x;

  fprintf(stderr,"%g", field[0]);
  for(x=1;x<p.N;x++) {
    fprintf(stderr,", %g", field[x]);
  }
  fprintf(stderr,"\n");
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
    f.Zx[x] = 0.0000;
    f.Zy[x] = 0.0000;
    f.Zz[x] = 0.0000;
    f.Ux[x] = 0.0000;
    f.Uy[x] = 0.0000;
    f.Uz[x] = 0.0000;
    f.Vx[x] = 0.0000;
    f.Vy[x] = 0.0000;
    f.Vz[x] = 0.0000;
    f.W[x] = 0.0000;

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
  free(f->Zx);
  free(f->Zy);
  free(f->Zz);
  free(f->Ux);
  free(f->Uy);
  free(f->Uz);
  free(f->Vx);
  free(f->Vy);
  free(f->Vz);
  free(f->W);
  free(f->kappa);
  free(f->p);
  free(f->pi);

  free(f->deltaM[0]);
  free(f->deltaM[1]);
  free(f->deltaM[2]);
  free(f->deltaM);

  free(f->U);
  free(f->V);
 
}

int iix(int x, int y, int z, hydro_params p) {
  return ((x+p.Lx)%p.Lx)*p.Ly*p.Lz + ((y+p.Ly)%p.Ly)*p.Lz + ((z+p.Lz)%p.Lz);
}

int **init_nb(hydro_params p) {
  int x, y, z;

  int **nb;

  // set up nb lookup table
  nb = (int **) malloc(p.N*sizeof(int *));
  
  for(x=0; x<p.Lx; x++) {
    for(y=0; y<p.Ly; y++) {
      for(z=0; z<p.Lz; z++) {
    
	nb[x*p.Ly*p.Lz + y*p.Lz + z] = (int *)malloc(6*sizeof(int));

	nb[x*p.Ly*p.Lz + y*p.Lz + z][0] = iix(x+1, y, z, p);
	nb[x*p.Ly*p.Lz + y*p.Lz + z][1] = iix(x-1, y, z, p);
        nb[x*p.Ly*p.Lz + y*p.Lz + z][2] = iix(x, y+1, z, p);
	nb[x*p.Ly*p.Lz + y*p.Lz + z][3] = iix(x, y-1, z, p);
	nb[x*p.Ly*p.Lz + y*p.Lz + z][4] = iix(x, y, z+1, p);
	nb[x*p.Ly*p.Lz + y*p.Lz + z][5] = iix(x, y, z-1, p);
      }
    }
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
  fprintf(stderr,"3D hydro port\n");

  if(argc != 2) {
    fprintf(stderr,"Usage: hydro <parameter file>\n");
    return 100;
  }

  // Parse params from stdin
  hydro_params p = get_parameters(argv[1]);

  

  p.N = p.Lx*p.Ly*p.Lz;

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

  double wmax;

  hydro_fields f;

  double initial_energy, current_energy;

  // Just runs malloc on all the fields therein
  alloc_fields(&f, p);

  fprintf(stderr, "- Allocated fields.\n");

  // For safety, set everything to zero
  zero_fields(f, p);

  fprintf(stderr, "- Zeroed fields.\n");

  // initial conditions
  /*  if(p.initial == INIT_BUBBLE) {
    create_gaussian_bubble(f, p);
  } else if(p.initial == INIT_SHOCK_TUBE) {
    create_shock_tube(f, p);
  } else {
    fprintf(stderr, "Unknown initial condition parameter!\n");
    exit(100);
  }
  */

  initial_3D(f,p);

  fprintf(stderr, "Initial conditions done\n");

  int **nb;

  nb = init_nb(p);

  fprintf(stderr, "- Initialised neighbour lookup table.\n");

  initial_energy = total_energy(f, nb, p);

  fprintf(stderr, "Initial avg energy per site: %g\n", 
	  initial_energy);

  // Output headers
  // printf("Step\ttime\tenergy\twallps\n");

  // Output of some field at every time step
  FILE *phi_fh = fopen("field.dat","w");
#ifdef SILO
  if(p.interval > 0)
    silo_init(p);
#endif // SILO

  // Step back for leapfrog initial conds
  evolve_backstep(f, nb, p);

  // Record how many lattice sites in output file, to keep self contained
  fwrite(&p.N, sizeof(int), 1, phi_fh);


  for(step = 0; step < p.steps; step++) {

    if((p.interval > 0) && (step % p.interval == 0)) {

#ifdef SILO
      write_silo_step(f, p, step);
#endif // SILO
      // Write time step, x coords and velocity field
      fwrite(&t, sizeof(double), 1, phi_fh);
      fwrite(f.E, sizeof(double), p.N, phi_fh);
      fwrite(f.p, sizeof(double), p.N, phi_fh);

      current_energy = total_energy(f, nb, p);

      
      fprintf(stderr,"%04d\t%6lf\t%6lf\t%6lf\t%6lf\n",
	     step, t,
	     current_energy, 
	     field_energy(f, nb, p), wmax);
      

      fflush(stdout);

      fprintf(stderr, "Energy violation: %lf%%\n",
	      100.0*fabs((current_energy-initial_energy)/initial_energy));
    }

    // Do field step
    //    evolve_field(f, nb, p);

    
    
    // Calculate EOS
    eq_of_state(f, p);

    // Do the hydro bits
    evolve_hydro(f, nb, p);


    // Advection of state variables
    donor_E(f, nb, p);
    //    transport_E_dir(f, nb, p, 0);

    // Advection of momentum
    //    donor_Z(f, nb, p);


    //    artificial_viscosity(f, nb, p);

    // Solve for T
    find_Ta(f, p);
    
    t += p.dt;

    wmax = get_gamma_max(f, p);

    
  } // main loop ends here

  
  for(x=0;x<p.Lx;x++) {
    fprintf(stdout,"%d %lf %lf\n", x, f.Vx[iix(x,0,0,p)], f.E[iix(x,0,0,p)]);
  }
  

  fclose(phi_fh);

  // Clean up memory
  free_nb(nb, p);
  free_fields(&f);
  
  return 0;
}
