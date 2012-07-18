/* hydro.h
 *
 * Header file for C port of 1D spherical hydro code.
 */
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <strings.h>
#include <malloc.h>

#ifdef SILO
#include <silo.h>
#endif // SILO

#ifdef PAPI
#include <papi.h>
#endif // PAPI

#ifdef MPI
#include <mpi.h>
#endif // MPI

#define INIT_SHOCK_TUBE 1
#define INIT_BUBBLE 2

// Composed directions
#define DIR_02 6
#define DIR_04 7
#define DIR_24 8
#define DIR_13 9
#define DIR_15 10
#define DIR_35 11
#define DIR_024 12
#define DIR_135 13


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
  int silointerval;
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

  int fieldN;

  int size;
  int rank;

#ifdef MPI

  int totalN;
  int slicex;
  int slicey;
  

  int inner;

  // Single haloes
  int offset_xM;
  int offset_xP;
  int offset_yM;
  int offset_yP;
  
  // Double haloes
  int offset_xMyM;
  int offset_xMyP;
  int offset_xPyM;
  int offset_xPyP;
    
  // Single haloes
  int inner_xM;
  int inner_xP;
  int inner_yM;
  int inner_yP;
  
  // Double haloes
  int inner_xMyM;
  int inner_xMyP;
  int inner_xPyM;
  int inner_xPyP;

  // Endposts
  int offset_xMyeM;
  int offset_xMyeP;
  int offset_xPyeM;
  int offset_xPyeP;
  int offset_yMxeM;
  int offset_yMxeP;
  int offset_yPxeM;
  int offset_yPxeP;



  // Ranks of neighbours
  int rank_xM;
  int rank_xP;
  int rank_yM;
  int rank_yP;

  int rank_xMyM;
  int rank_xMyP;
  int rank_xPyM;
  int rank_xPyP;





  int myposx;
  int myposy;

  int **inverse;

#endif // MPI
 



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
  double *W;

  double **V;
  double **U;
  double **F;
  double **Z;

} hydro_fields;



// main.c

int iix(int x, int y, int z, hydro_params p);



// arrangement.c
void layout(hydro_params *p);
int **init_nb(hydro_params *p);
int get_x(int n, hydro_params p);
int get_y(int n, hydro_params p);
int get_z(int n, hydro_params p);
void halo_field(double *field, hydro_params p);
double reduce_sum(double result, hydro_params p);
double reduce_max(double result, hydro_params p);

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

// papi.c

#ifdef PAPI
void papi_init();
void papi_finalise();
#endif // PAPI

// parameters.c

void get_parameters(char *filename, hydro_params *p);

// silage.c

#ifdef SILO
void silo_init(hydro_params p);
void write_silo_step(hydro_fields f, hydro_params p, int step);
#endif // SILO


// util.c

double minof3(double a, double b, double c);
double maxof3(double a, double b, double c);
double minof2(double a, double b);
