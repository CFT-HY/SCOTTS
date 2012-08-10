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

#ifndef MPI
#error Cannot use FFTW3 without MPI - local FFTs not implemented!
#endif // !MPI

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
 * For the uij's, there are six components.
 */
#define TENSOR_CPTS 6

#define CPT_11 0
#define CPT_21 1
#define CPT_31 2
#define CPT_22 3
#define CPT_32 4
#define CPT_33 5



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

  // NB: artificial viscosity not implemented yet!
  double Cav;

  // C aka eta, the field viscosity
  double C;

  // Input parameters for the potential
  double Lheat;
  double sigma;
  double lcorr;

  // Calculated potential parameters
  double alpha;
  double lambda;
  double gamma;

  // How often to deal with output
  int interval;
  int silointerval;

  // Initial conditions type (see #defines above)
  int initial;

  // The following parameters are preset, but changeable in main.c
  double Tconst;
  double T0;

  // Degrees of freedom
  double a;

  // Number of physical sites in volume
  int N;

  // Sites in local area (including halo)
  int fieldN;

  // How many nodes, and which are we
  int size;
  int rank;

  //  double comms_time;

  // How many sites in each direction are unique to us
  int slicex;
  int slicey;
  // -- i.e. we have (slicex+2)*(slicey+2)*Lz sites to ourselves

  // Where do our sites start in the physical volume?
  int shiftx;
  int shifty;
  // NB: the physical position of a site x is e.g. (x+shiftx-1)


  // Where the silo files go
  char silodir[500];


#ifdef MPI

  // Ranks of neighbours
  int rank_xM;
  int rank_xP;
  int rank_yM;
  int rank_yP;

  // Ranks of double neighbours
  int rank_xMyM;
  int rank_xMyP;
  int rank_xPyM;
  int rank_xPyP;

  // Where am I in the domain decomposition?
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

  double ****uij;
  double ****udotij;


} hydro_fields;



// main.c just contains main()


// alloc.c
double ***make_field(hydro_params p);
double ****make_vector(hydro_params p);
double ****make_tensor(hydro_params p);

void free_field(hydro_params p, double ***field);
void free_vector(hydro_params p, double ****vector);
void free_tensor(hydro_params p, double ****tensor);

void alloc_fields(hydro_fields *f, hydro_params p);
void zero_fields(hydro_fields f, hydro_params p);
void free_fields(hydro_fields *f, hydro_params p);


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
void evolve_uij(hydro_fields f, hydro_params p);

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
void fft_field(hydro_fields f, hydro_params p);
void fft_tensor(hydro_fields f, hydro_params p);
#endif // FFT



void stress_energy(hydro_fields f, hydro_params p, double ****Tij);
double tzerozero(hydro_fields f, hydro_params p);


