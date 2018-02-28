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
#include <unistd.h>
#include <stdarg.h>

// Apple doesn't seem to believe in malloc.h
#ifndef __APPLE__
#include <malloc.h>
#define HAVE_MALLOC_H
#endif // ! __APPLE__

// Store time slice data in the Silo format
#ifdef SILO
#include <silo.h>

#define CPMODE DB_HDF5
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

// #define FFT_DEBUG

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


#define NUC_OFF 0
#define NUC_EXP 1
#define NUC_LIST 2
#define NUC_FILE 3


#define GW_BOTH 1
#define GW_FIELD 2
#define GW_FLUID 3

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

  float dx;
  float dt;

  int Lx, Ly, Lz;

  int steps;

  // NB: artificial viscosity not implemented yet!
  float Cav;

  // C aka eta, the field viscosity
  float C;

  // Input parameters for the potential
  /*
  float Lheat;
  float sigma;
  float lcorr;
  */

  // Calculated potential parameters
  float alpha;
  float lambda;
  float gamma;
  // Temperature parameters
  float Tconst;
  float T0;

#ifdef INITPS
  // Cutoff time
  float initps;

  float initcutoff;
  float initlength;
#endif

#ifdef CUTOFF
  // Cutoff time
  float tcutoff;
#endif

  // How often to deal with output
  int interval;
  int fftinterval;
  int silointerval;
  int checkpointinterval;

  // x coord to slice through for write_silo_slice_step
  int siloslicecoord;
  
  int uetcstart;

  // Initial conditions type (see #defines above)
  int initial;



  // Scale factor
  float a;

  // Degrees of freedom
  float gstar, gdeg;

  // Number of physical sites in volume
  int N;

  // Sites in local area (including halo)
  int fieldN;

  // How many nodes, and which are we
  int size;
  int rank;

  //  float comms_time;

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

  // Where checkpoint files go
  char checkpointdir[500];

  int bubbles;

  // How we nucleate
  int nucleation;

  int *nucsteps;
  int n_nucsteps;

  float beta;

  // to rescale bubble size
  float scale;

  int Linv;

  int seed;

  int gwsource;

#ifdef MPI

  // Ranks of neighbours
  int rank_xM;
  int rank_xP;
  int rank_yM;
  int rank_yP;

  // Ranks of float neighbours
  int rank_xMyM;
  int rank_xMyP;
  int rank_xPyM;
  int rank_xPyP;

  // Where am I in the domain decomposition?
  int myposx;
  int myposy;

#endif // MPI
 

} hydro_params;



// See alloc_fields in alloc.c for details
typedef struct {

  float *x_inv;
  float *phi_inv;
  float *V_inv;

  // Scalar field
  float ***phi;
  float ***pifull;
  float ***pi;
  float ***phiold;

#ifndef SCALAR
  // Fluid
  float ***T;
  float ***E;
  float ***W;
  float ***kappa;
  float ***p;

  float ****V;
  float ****Z;
  float ****U;
  float ****F;
#endif // SCALAR

  // Gravity
  float ****uij;
  float ****udotij;

  // (Only for UETCs)
  float ****initial_Tij;


} hydro_fields;



// main.c just contains main()


// alloc.c
float ***make_field(hydro_params p);
float ****make_vector(hydro_params p);
float ****make_tensor(hydro_params p);


void free_field(hydro_params p, float ***field);
void free_vector(hydro_params p, float ****vector);
void free_tensor(hydro_params p, float ****tensor);

void alloc_fields(hydro_fields *f, hydro_params p);
void zero_fields(hydro_fields f, hydro_params p);
void free_fields(hydro_fields *f, hydro_params p);


// arrangement.c
void layout(hydro_params *p);
int get_x(int n, hydro_params p);
int get_y(int n, hydro_params p);
int get_z(int n, hydro_params p);
void halo_field(float ***field, hydro_params p);
float reduce_sum(float result, hydro_params p);
int reduce_sum_int(int result, hydro_params p);
float reduce_max(float result, hydro_params p);
int reduce_max_int(int result, hydro_params p);
float reduce_min(float result, hydro_params p);
int reduce_and(int result, hydro_params p);
void init_comms_time(hydro_params *p);
float get_comms_time(hydro_params *p);
void printf0(hydro_params p, char *msg, ...);
void die(int howbad);

// checkpoint.c
void checkpoint(hydro_fields f, hydro_params p, int step);
int usable_checkpoint(hydro_fields f, hydro_params p);
int load_checkpoint(hydro_fields f, hydro_params p);

// evolve.c
void evolve_backstep(hydro_fields f, hydro_params p);
void evolve_field(hydro_fields f, hydro_params p);
void evolve_hydro(hydro_fields f, hydro_params p);
#ifdef CUTOFF
void evolve_uij(hydro_fields f, hydro_params p, float cutoff);
#else
void evolve_uij(hydro_fields f, hydro_params p);
#endif

// Not implemented yet
void artificial_viscosity(hydro_fields f, int **nb, hydro_params p);

// potential.c
float Vf(hydro_params p, float T, float this_phi);
float Vdf(hydro_params p, float T, float this_phi);
float VTf(hydro_params p, float T, float this_phi);
float VTTf(hydro_params p, float T, float this_phi);
void Vpot(hydro_params p, float *T, float *phi, float *Vprecalc);
void Vdpot(hydro_params p, float *T, float *phi, float *Vprecalc);


// energy.c
float field_energy(hydro_fields f, hydro_params p);
float gradient_energy(hydro_fields f, hydro_params p);
float total_energy(hydro_fields f, hydro_params p);
float kinetic_energy(hydro_fields f, hydro_params p);
float rest_energy(hydro_fields f, hydro_params p);
void energy_density(hydro_fields f, hydro_params p, float ***en);
void stress_energy(hydro_fields f, hydro_params p, float ****Tij);
float avg_pressure(hydro_fields f, hydro_params p);
float tzerozero(hydro_fields f, hydro_params p);

// eos.c

void find_Ta(hydro_fields f, hydro_params p);
void eq_of_state(hydro_fields f, hydro_params p);

// transport.c

void advect_E(hydro_fields f, hydro_params p);
void advect_Z(hydro_fields f, hydro_params p);

// initial.c
void initial_scalar_bubble(hydro_fields f, hydro_params *p_ptr);
void initial_blank(hydro_fields f, hydro_params p);
int safe_distance(hydro_fields f, hydro_params p);
int can_nucleate(hydro_fields f, hydro_params p, int x0, int y0, int z0);
void nucleate_at(hydro_fields f, hydro_params p, int x0, int y0, int z0);
void initial_3D(hydro_fields f, hydro_params p);
int do_nucleate(hydro_fields f, hydro_params *p_ptr);
int should_nucleate(hydro_fields f, hydro_params p, float t, int step);
void init_profile(hydro_fields *f, hydro_params *p);

// output.c
float get_gamma_max(hydro_fields f, hydro_params p);
float get_veltot(hydro_fields f, hydro_params p);
void dump(float *field, hydro_params p);
void histo_field(float ***field, hydro_params p, int step);
void didj(float *cpts, hydro_fields f, hydro_params p);

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

void write_silo_slice_step(hydro_fields f, hydro_params *p_ptr, int step);
#endif // SILO


// util.c
float minof3(float a, float b, float c);
int minof3_int(int a, int b, int c);
float maxof3(float a, float b, float c);
float minof2(float a, float b);


#ifdef FFT
// fft.c
void fft_field(hydro_fields f, hydro_params p, float ***field, int step);

// uetc.c
void init_uetc(hydro_fields f, hydro_params p);
void fft_uetc(hydro_fields f, hydro_params p, int step);

// gw.c
float proj(int T, float kx, float ky, float kz);
float fft_tensor(hydro_fields f, hydro_params p, int step,
		  float energydensity);
int indexof(int i, int j);
float lambda(int i, int j, int l, int m, float kx, float ky, float kz);

// velps.c
float vel_proj(int T, float kx, float ky, float kz);
void fft_vel(hydro_fields f, hydro_params p, int step, float ****vectorfield);

#endif // FFT


#ifdef INITPS

// initps.c
void init_ps(hydro_fields f, hydro_params p, float ****field);

#endif // INITPS
