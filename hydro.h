/* hydro.h
 *
 * Header file for C port of 1D spherical hydro code.
 */
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <strings.h>

#ifdef SILO
#include <silo.h>
#endif // SILO

#define INIT_SHOCK_TUBE 1
#define INIT_BUBBLE 2

typedef struct {
  double dx;
  double dt;
  int Lx, Ly, Lz;
  int steps;

  double Cav;
  double C;

  double Lheat;
  double sigma;
  double lcorr;

  int interval;
  int initial;

  // The following parameters are calculated
  double xxwall;
  double Tconst;

  double alpha;
  double lambda;
  double gamma;

  double a;
  int N;
  double T0;


} hydro_params;

typedef struct {
  double *phi;
  double *pi;
  double *pifull;
  double *xe;
  double *xc;
  double *T;
  double *kappa;
  double *phiold;
  double *p;
  double *E;
  double *Zx;
  double *Zy;
  double *Zz;
  double *Ux;
  double *Uy;
  double *Uz;
  double *Vx;
  double *Vy;
  double *Vz;
  double *W;

  double **V;
  double **U;
  double **deltaM;
} hydro_fields;

// main.c

int iix(int x, int y, int z, hydro_params p);

// evolve.c

void evolve_backstep(hydro_fields f, int **nb, hydro_params p);
void evolve_field(hydro_fields f, int **nb, hydro_params p);
void evolve_hydro(hydro_fields f, int **nb, hydro_params p);
void artificial_viscosity(hydro_fields f, int **nb, hydro_params p);

// potential.c

double Vf(hydro_params p, double T, double this_phi);
double Vdf(hydro_params p, double T, double this_phi);
double VTf(hydro_params p, double T, double this_phi);
double VTTf(hydro_params p, double T, double this_phi);
void Vpot(hydro_params p, double *T, double *phi, double *Vprecalc);
void Vdpot(hydro_params p, double *T, double *phi, double *Vprecalc);


// energy.c

double field_energy(hydro_fields f, int **nb, hydro_params p);
double total_energy(hydro_fields f, int **nb, hydro_params p);


// eos.c

void find_Ta(hydro_fields f, hydro_params p);
void eq_of_state(hydro_fields f, hydro_params p);

// transport.c

void advect_E(hydro_fields f, int **nb, hydro_params p);
void advect_Z(hydro_fields f, int **nb, hydro_params p);

// initial.c
/*
void create_1D_bubble(hydro_fields f, hydro_params p);
double psibar(double x, double lbar);
void create_gaussian_bubble(hydro_fields f, hydro_params p);
void create_shock_tube(hydro_fields f, hydro_params p);
*/
void initial_3D(hydro_fields f, hydro_params p);

// output.c

double wallpos(hydro_fields f, hydro_params p);
double get_gamma_max(hydro_fields f, hydro_params p);

// parameters.c

hydro_params get_parameters(char *filename);

// silage.c

#ifdef SILO
void silo_init(hydro_params p);
void write_silo_step(hydro_fields f, hydro_params p, int step);
#endif // SILO
