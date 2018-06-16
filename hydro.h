/** @file hydro.h
 *
 * Header file for hydro+scalar code.
 *
 * Contains the structs `hydro_params` and `hydro_fields`.
 * Also contains prototypes for all functions in the code.
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

#ifdef DOUBLE
typedef double Real;
#else
typedef float Real;
#endif

/*
 * Not sure these are useful in 3D, but a way of enumerating choices of the
 * initial conditions and communicating that between the parameter input
 * step and the initial conditions step.
 */
#define INIT_SHOCK_TUBE 1
#define INIT_BUBBLE 2
#define INIT_PS 3

#define NUC_OFF 0
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



/** Struct containing parameters that are not changed during the
 * simulation.
 *
 * Stores the parameters, including:
 * - physics (supplied and derived quantities)
 * - lattice geometry
 * - communication (MPI) details
 */

typedef struct {

  /** Lattice spacing.
   */
  Real dx;

  /** Timestep interval.
   */
  Real dt;

  /** Number of lattice points in `x` direction.
   */
  int Lx;
  /** Number of lattice points in `y` direction.
   */
  int Ly;
  /** Number of lattice points in `z` direction.
   */
  int Lz;

  /* Number of timesteps in the simulation.
   */
  int steps;

  /** Artificial viscosity parameter.
   *
   * NB: artificial viscosity not implemented yet!
   */
  Real Cav;

  /** C aka \f$ \eta \f$, the field viscosity
   */
  Real C;

  /** Cubic potential parameter \f$ A \f$.
   */
  Real alpha;

  /** Quartic potential parameter \f$ \lambda \f$.
   */
  Real lambda;

  /** Quadratic potential parameter \f$ \gamma \f$.
   */
  Real gamma;
  
  /** Nucleation temperature \f$ T_N \f$ used to initialise the
   *  simulation.
   */
  Real Tconst;

  /** Potential parameter \f$ T_0 \f$ which corresponds to the
   *  temperature at which the transition is no longer first order.
   */
  Real T0;

  Real initnorm;
  Real initcutoff;
  Real initlength;

  /** How frequently to write global outputs to stderr.
   */
  int interval;

  /** How often to write gravitational wave and velocity power
   *  spectra.
   */
  int fftinterval;

  /** How frequently to write fields to a silo file.
   */
  int silointerval;

  /** How frequently to write slices through the simulation box to a
      silo file.
   */
  int silosliceinterval;
  int checkpointinterval;
  
  /** `x` coord to slice through for write_silo_slice_step()
   */
  int siloslicecoord;

  /** When to start performing uetc calculations.
   */
  int uetcstart;

  /** Initial conditions type (see #defines above)
   */
  int initial;  

  /** Scale factor.
   */
  Real a;

  /** Effective degrees of freedom.
   */
  Real gstar;

  /** `gdeg` is `gstar` scaled by \f$ (\pi^2/90) \f$.
   */
  Real gdeg;

  /** Number of physical sites in volume.
   */
  int N;

  /** Sites in local area (including halo).
   */
  int fieldN;

  /** How many nodes we are using.
   */
  int size;

  /** Our node identifier.
   */
  int rank;

  //  Real comms_time;

  /** How many sites in the `x` direction are unique to us.
   *
   * -- i.e. we have (slicex+2)*(slicey+2)*Lz sites to ourselves
   */
  int slicex;
  /** How many sites in the `y` direction are unique to us.
   *
   * -- i.e. we have `(slicex+2)*(slicey+2)*Lz` sites to ourselves
   */
  int slicey;
  

  /** `shiftx` corresponds to where the local `x` coordinate starts in
   * the physical volume.
   *
   * NB: the physical position of a site `x` is e.g. `(x+shiftx-1)`
   */
  int shiftx;

  /** `shiftx` corresponds to where the local `y` coordinate starts in
   * the physical volume.
   *
   * NB: the physical position of a site `y` is e.g. `(y+shifty-1)`
   */
  int shifty;

  /** Global output filename.
   *
   * Either a filename or "stdout".
   */
  char output_fname[500];

  /** Global output destination.
   *
   * Points to a file specified at runtime, either a filename or stdout.
   */
  FILE *outputdest;
  
  /** Where the silo files go.
   */
  char silodir[500];

  /** Where checkpoint files go
   */
  char checkpointdir[500];

  /** Number of bubbles to spawn on the first timestep.
   */
  int bubbles;

  /** Nucleation type.
   */
  int nucleation;

  /** Steps on which to perform nucleation.
   *
   * NB: The same step can appear more than once to indicate to try to
   * nucleate multiple bubbles on that timestep.
   */
  int *nucsteps;
  /** Total number of steps on which we perform nucleation. 
   *
   * NB: If a step is repeated it is counted multiple times.
   */
  int n_nucsteps;

  /** Used to rescale bubble size
   */
  Real scale;

  /** Random seed so it is kept the same across nodes.
   */
  int seed;

  /** What is used to source gravitational waves? (fluid, scalar or
   *  both)
   */
  int gwsource;

#ifdef MPI

  /** Rank of neighbour in negative `x` direction.
   */
  int rank_xM;
  /** Rank of neighbour in positive `x` direction.
   */
  int rank_xP;
  /** Rank of neighbour in negative `y` direction.
   */
  int rank_yM;
  /** Rank of neighbour in positive `x` direction.
   */
  int rank_yP;

  /** Rank of neighbour in negative `x` and `y` direction.
   */
  int rank_xMyM;
  /** Rank of neighbour in negative `x` and positive `y` direction.
   */
  int rank_xMyP;
  /** Rank of neighbour in positive `x` and negative `y` direction.
   */
  int rank_xPyM;
  /** Rank of neighbour in positive `x` and `y` direction.
   */
  int rank_xPyP;

  /** Where am I in the domain decomposition along `x` direction?
   */
  int myposx;
  /** Where am I in the domain decomposition along `y` direction?
   */
  int myposy;

#endif // MPI

  /** Surface tension \f$ \sigma \f$ for the bubble.
   *
   *  \f[
   *  \sigma=\frac{2\sqrt{2}}{81}\frac{A^3}{\lambda^{5/2}}T_c^3
   *  \f]
   */
  Real surface_tension;

  /** Value of `phi` at the center of the nucleated bubble.
   *
   * Read from the parameter file, if given as <=0 default to broken 
   * phase value.
   */
  Real phimin;

  /** Critical radius of the bubble \f$ R_c \f$.
   *
   * Read from the parameter file, if given as <=0 default to:
   * \f[
   * R_c= -\dfrac{2\sigma}{-V(\phi_b,T_N)}
   * \f]
   */
  Real R_critical;

  /** Scaled version of the critical radius.
   *
   * Used to stop bubbles from collapsing when necessary.
   */
  Real R_scaled;

  /** Constant in potential term.
   *
   * Calculated s.t either $V(\phi_b,T=0)=0$ or if `TINDEP` flag
   * declared $V(\phi_b)=0$.
   */
  Real V0;

  
} hydro_params;



/** Struct containing the fields defined on the lattice that we track
 * at every timestep.
 * 
 * Scalar field fields:
 *
 * `phi`, `phi_old`, `pi`, `pi_future`.
 *
 * Fluid fields:
 *
 * `T`, `E`, `W`, `kappa`, `p`, `V`, `Z`, `U`, and `F`.
 *
 *
 * Gravity fields:
 *
 * `uij`, `udotij`, and `initial_Tij`.
 * 
 */
typedef struct {

  Real *x_inv;
  Real *phi_inv;
  Real *V_inv;

  // Scalar field

  /** `phi` is the scalar field on the current timestep.
   */
  Real ***phi;

  /** `pi_future` is the value of `pi` on the same timestep as `phi`.
   */
  Real ***pi_future;

  /** `pi` is the conjugate momenta of the scalar field.
   */
  Real ***pi;

  /** `phi_old` is the scalar field on the last timestep.
   */
  Real ***phi_old;

#ifndef SCALAR
  // Fluid
  
  /** `T` is the temperature of the fluid.
   */
  Real ***T;

  /** `E` is the internal energy density of the fluid.
   */
  Real ***E;

  /** `W` is the zone-centred boost.
   */
  Real ***W;

  /** `kappa` is the adiabatic index for the system.
   */
  Real ***kappa;

  /** `p` is the pressure of the fluid. Obtained by eq_of_state() and
   * used in hydro.
   */
  Real ***p;

  /** `V` is the fluid 3-velocity.
   */
  Real ****V;

  /** `Z` is the fluid momentum density.
   */
  Real ****Z;

  /** `U` is the fluid 4-velocity.
   */
  Real ****U;

  /** `F` is a temporary variable used in advection.
   */
  Real ****F;
#endif // SCALAR

  // Gravity

  /** `uij` are the unprojected metric peturbations.
   */
  Real ****uij;

  /** `udotij` are the time derivative of the `uij` and are used for
   *    GW power spectrum.
   */
  Real ****udotij;

  // (Only for UETCs)

  /** `initial_Tij` is only used for calculating unequal time
   *    correlators.
   */
  Real ****initial_Tij;


} hydro_fields;

/** Helper struct for performing MPI_MAXLOC and MPI_MINLOC operations.
 *
 * Useful for debugging, e.g finding max/min value of something on the
 * lattice and the processor containing the lattice site.
 */
typedef struct{
  Real value;
  int rank;
} value_rank;

// main.c just contains main()


// alloc.c
Real ***make_field(hydro_params p);
Real ****make_vector(hydro_params p);
Real ****make_tensor(hydro_params p);


void free_field(hydro_params p, Real ***field);
void free_vector(hydro_params p, Real ****vector);
void free_tensor(hydro_params p, Real ****tensor);

void alloc_fields(hydro_fields *f, hydro_params p);
void zero_fields(hydro_fields f, hydro_params p);
void free_fields(hydro_fields *f, hydro_params p);


// comms.c
void layout(hydro_params *p);
int get_x(int n, hydro_params p);
int get_y(int n, hydro_params p);
int get_z(int n, hydro_params p);
void halo_field(Real ***field, hydro_params p);
Real reduce_sum(Real result, hydro_params p);
int reduce_sum_int(int result, hydro_params p);
Real reduce_max(Real result, hydro_params p);
int reduce_max_int(int result, hydro_params p);
Real reduce_min(Real result, hydro_params p);
int reduce_and(int result, hydro_params p);
value_rank reduce_maxloc(Real result, hydro_params p);
value_rank reduce_minloc(Real result, hydro_params p);
void init_comms_time(hydro_params *p);
Real get_comms_time(hydro_params *p);
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
void evolve_uij(hydro_fields f, hydro_params p);

// Not implemented yet
void artificial_viscosity(hydro_fields f, int **nb, hydro_params p);

// potential.c
Real Vf(hydro_params p, Real T, Real this_phi);
Real Vdf(hydro_params p, Real T, Real this_phi);
Real VTf(hydro_params p, Real T, Real this_phi);
Real VTTf(hydro_params p, Real T, Real this_phi);
void Vpot(hydro_params p, Real ***T, Real ***phi, Real ***Vprecalc);
void Vdpot(hydro_params p, Real ***T, Real ***phi, Real ***Vprecalc);


// energy.c
Real field_energy(hydro_fields f, hydro_params p);
Real gradient_energy(hydro_fields f, hydro_params p);
Real total_energy(hydro_fields f, hydro_params p);
Real kinetic_energy(hydro_fields f, hydro_params p);
Real rest_energy(hydro_fields f, hydro_params p);
void energy_density(hydro_fields f, hydro_params p, Real ***en);
void stress_energy(hydro_fields f, hydro_params p, Real ****Tij);
Real avg_pressure(hydro_fields f, hydro_params p);
Real tzerozero(hydro_fields f, hydro_params p);

// eos.c

void find_Ta(hydro_fields f, hydro_params p);
void eq_of_state(hydro_fields f, hydro_params p);

// transport.c

void advect_E(hydro_fields f, hydro_params p);
void advect_Z(hydro_fields f, hydro_params p);

// initial.c
void initial_blank(hydro_fields f, hydro_params p);
int safe_distance(hydro_fields f, hydro_params p);
int can_nucleate(hydro_fields f, hydro_params p, int x0, int y0, int z0);
void nucleate_at(hydro_fields f, hydro_params p, int x0, int y0, int z0);
void initial_3D(hydro_fields f, hydro_params p);
int try_nucleate(hydro_fields f, hydro_params p);
int bubbles_at_step(hydro_fields f, hydro_params p, Real t, int step);
void init_profile(hydro_fields *f, hydro_params *p);

// output.c
Real get_gamma_max(hydro_fields f, hydro_params p);
Real get_s_max(hydro_fields f, hydro_params p);
Real get_veltot(hydro_fields f, hydro_params p);
void dump(Real *field, hydro_params p);
void histo_field(Real ***field, hydro_params p, int step);
void didj(Real *cpts, hydro_fields f, hydro_params p);

// papi.c
#ifdef PAPI
void papi_init();
void papi_finalise();
#endif // PAPI

// parameters.c
void get_parameters(char *filename, hydro_params *p);
void set_bubble_parameters(hydro_params *p);

// silage.c
#ifdef SILO
void silo_init(hydro_params p);
void write_silo_step(hydro_fields f, hydro_params p, int step);

void write_silo_slice_step(hydro_fields f, hydro_params p, int step);
#endif // SILO


// util.c
Real minof3(Real a, Real b, Real c);
int minof3_int(int a, int b, int c);
Real maxof3(Real a, Real b, Real c);
Real minof2(Real a, Real b);


#ifdef FFT
// fft.c
void fft_field(hydro_fields f, hydro_params p, Real ***field, int step);

// uetc.c
void init_uetc(hydro_fields f, hydro_params p);
void fft_uetc(hydro_fields f, hydro_params p, int step);

// gw.c
Real proj(int T, Real kx, Real ky, Real kz);
Real fft_tensor(hydro_fields f, hydro_params p, int step,
		  Real energydensity);
int indexof(int i, int j);
Real lambda(int i, int j, int l, int m, Real kx, Real ky, Real kz);

// velps.c
Real vel_proj(int T, Real kx, Real ky, Real kz);
void fft_vel(hydro_fields f, hydro_params p, int step, Real ****vectorfield);

#endif // FFT


#ifndef SCALAR
// initps.c
void init_ps(hydro_fields f, hydro_params p, Real ****field);
void norm_power(hydro_fields f, hydro_params p, Real ****field);
Real get_momtot(hydro_fields f, hydro_params p);

#endif // !SCALAR
