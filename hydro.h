/* hydro.h
 *
 * Header file for hydro+scalar code.
 */

#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <strings.h>
#include <time.h>

// Apple doesn't seem to believe in malloc.h
#ifndef __APPLE__
#include <malloc.h>
#define HAVE_MALLOC_H
#endif // ! __APPLE__

// Store time slice data in the Silo format
#ifdef SILO
#include <silo.h>
#endif // SILO

// Instrumentation using PAPI
#ifdef PAPI
#include <papi.h>
#endif // PAPI

// Parallelism with MPI
#ifdef MPI
#include <mpi.h>
#endif // MPI

#ifdef FFT
#include <fftw3-mpi.h>
#endif // FFT


/*
 * Not sure these are useful in 3D, but a way of enumerating choices of the
 * initial conditions and communicating that between the parameter input
 * step and the initial conditions step.
 */
#define INIT_SHOCK_TUBE 1
#define INIT_BUBBLE 2



/*
 * Composed directions
 *
 * In the neighbour lookup tables, we have nb[site][2*dim] for positive
 * directions and nb[site][2*dim+1] for negative directions. To handle
 * haloing gracefully, and avoid having invalid values of nb, we explicitly
 * store the diagonal neighbours.
 *
 * It involves storing about 4*sizeof(double) (=8*sizeof(int)) extra
 * information at each site; that is small fry compared to the
 * ~21*sizeof(double) we're already storing (see struct hydro_fields).
 */

#define DIR_02 6
#define DIR_04 7
#define DIR_24 8
#define DIR_13 9
#define DIR_15 10
#define DIR_35 11
#define DIR_024 12
#define DIR_135 13


/*
 * struct hydro_params
 *
 * Stores the parameters, including:
 * - physics (supplied and derived quantities)
 * - lattice geometry
 * - communication (MPI) details
 */

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

  //  double comms_time;

  int slicex;
  int slicey;

  int shiftx;
  int shifty;

  char silodir[500];

#ifdef MPI

  int totalN;
  int inner;

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

#endif // MPI
 



} hydro_params;

typedef struct {
  double ***phi;
  double ***pi;
  double ***pifull;
  double ***T;
  double ***kappa;
  double ***phiold;
  double ***p;
  double ***E;
  double ***W;

  double ****V;
  double ****U;
  double ****F;
  double ****Z;


  double *phi_root;
  double *pi_root;
  double *pifull_root;
  double *T_root;
  double *kappa_root;
  double *phiold_root;
  double *p_root;
  double *E_root;
  double *W_root;

  double **V_root;
  double **U_root;
  double **F_root;
  double **Z_root;

} hydro_fields;



// main.c
void alloc_fields(hydro_fields *f, hydro_params p);
void zero_fields(hydro_fields f, hydro_params p);
void free_fields(hydro_fields *f, hydro_params p);
void free_field(hydro_params p, double ***field);
double ***make_field(hydro_params p);
double ****make_vector(hydro_params p);
void free_vector(hydro_params p, double ****vector);

// arrangement.c
void layout(hydro_params *p);
int get_x(int n, hydro_params p);
int get_y(int n, hydro_params p);
int get_z(int n, hydro_params p);
void halo_field(double ***field, hydro_params p);
double reduce_sum(double result, hydro_params p);
double reduce_max(double result, hydro_params p);
int reduce_and(int result, hydro_params p);
void init_comms_time(hydro_params *p);
double get_comms_time(hydro_params *p);

int **init_inverse(hydro_params *p);
int **init_nb(hydro_params *p);

// evolve.c

void evolve_backstep(hydro_fields f, hydro_params p);
void evolve_field(hydro_fields f, hydro_params p);
void evolve_hydro(hydro_fields f, hydro_params p);
// void artificial_viscosity(hydro_fields f, int **nb, hydro_params p);

// potential.c

double Vf(hydro_params p, double T, double this_phi);
double Vdf(hydro_params p, double T, double this_phi);
double VTf(hydro_params p, double T, double this_phi);
double VTTf(hydro_params p, double T, double this_phi);
void Vpot(hydro_params p, double *T, double *phi, double *Vprecalc);
void Vdpot(hydro_params p, double *T, double *phi, double *Vprecalc);


// energy.c

double field_energy(hydro_fields f, hydro_params p);
double total_energy(hydro_fields f, hydro_params p);
void energy_density(hydro_fields f, hydro_params p, double ***en);


// eos.c

void find_Ta(hydro_fields f, hydro_params p);
void eq_of_state(hydro_fields f, hydro_params p);

// transport.c

void advect_E(hydro_fields f, hydro_params p);
void advect_Z(hydro_fields f, hydro_params p);

// initial.c
/*
void create_1D_bubble(hydro_fields f, hydro_params p);
double psibar(double x, double lbar);
void create_gaussian_bubble(hydro_fields f, hydro_params p);
void create_shock_tube(hydro_fields f, hydro_params p);
*/
void initial_3D(hydro_fields f, hydro_params p);
void initial_blank(hydro_fields f, hydro_params p);
void nucleate_at(hydro_fields f, hydro_params p, int x0, int y0, int z0);
int can_nucleate(hydro_fields f, hydro_params p, int x0, int y0, int z0);

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

// fft.c
#ifdef FFT
void fft(hydro_fields f, hydro_params p);
#endif // FFT
