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
#ifdef USE_MPI
#include <mpi.h>
#endif // USE_MPI

#ifdef FFT

// #define FFT_DEBUG

#ifndef USE_MPI
#error Cannot use FFTW3 without MPI - local FFTs not implemented!
#endif // !USE_MPI

#include <fftw3-mpi.h>
#endif // FFT

/*
 * Not sure these are useful in 3D, but a way of enumerating choices of the
 * initial conditions and communicating that between the parameter input
 * step and the initial conditions step.
 */
#define INIT_SHOCK_TUBE 1
#define INIT_BUBBLE 2
#define INIT_PS 3
#define INIT_FLUID_SPHERE 4

#define NUC_OFF 0
#define NUC_LIST 1
#define NUC_FILE 2
#define NUC_FILE_LOC 3

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


#define INITPSFILE_DIV 0
#define INITPSFILE_ROT 1
#define INITPSFILE_ALL 2

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
  double dx;

  /** Timestep interval.
   */
  double dt;

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
  double Cav;

  /** C aka \f$ \eta \f$, the field viscosity
   */
  double C;

  /** Cubic potential parameter \f$ A \f$.
   */
  double alpha;

  /** Quartic potential parameter \f$ \lambda \f$.
   */
  double lambda;

  /** Quadratic potential parameter \f$ \gamma \f$.
   */
  double gamma;

  /** Nucleation temperature \f$ T_N \f$ used to initialise the
   *  simulation.
   */
  double Tconst;

  /** Potential parameter \f$ T_0 \f$ which corresponds to the
   *  temperature at which the transition is no longer first order.
   *  (Not used in BAG model).
   */
  double T0;

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
  double a;

  /** Effective degrees of freedom.
   */
  double gstar;

  /** `gdeg` is `gstar` scaled by \f$ (\pi^2/90) \f$.
   */
  double gdeg;

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

  //  double comms_time;

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

  /** Where the silo files go.
   */
  char silodir[500];

  /** Where checkpoint files go
   */
  char checkpointdir[500];

  /** Where to find the initial power spectrum
   */
  char initpsfile[500];
  int initpsfile_type;

  /** Number of bins.
   */
  int initpsbins;

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

  /** Bubble nucleation center locations on the lattice.
   *
   * NB: Only used if nucleation type is `NUC_FILE_LOC`.
   */
  int **nuclocs;

  /** Total number of steps on which we perform nucleation.
   *
   * NB: If a step is repeated it is counted multiple times.
   */
  int n_nucsteps;

  /** Used to rescale bubble size
   */
  double scale;

  /** Random seed so it is kept the same across nodes.
   */
  int seed;

  /** What is used to source gravitational waves? (fluid, scalar or
   *  both)
   */
  int gwsource;

#ifdef USE_MPI

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

#endif // USE_MPI

  /** Surface tension \f$ \sigma \f$ for the bubble.
   *
   *  \f[
   *  \sigma=\frac{2\sqrt{2}}{81}\frac{A^3}{\lambda^{5/2}}T_c^3
   *  \f]
   *
   *  For bag model
   *  \f[
   *  \sigma=\frac{\left(\alpha \left(\alpha + \sqrt{\alpha^2 - 4 \gamma \lambda}
   *                                  \right) -2\gamma \lambda \right)^{3/2}}
   *              {24 \lambda^{5/2}}
   *  \f]
   */
  double surface_tension;

  /** Value of `phi` at the center of the nucleated bubble.
   *
   * Read from the parameter file, if given as <=0 default to broken
   * phase value at \f$T=T_N\f$.
   */
  double phimin;

  /** Critical radius of the bubble \f$ R_c \f$.
   *
   * Read from the parameter file, if given as <=0 default to:
   * \f[
   * R_c= -\dfrac{2\sigma}{-V(\phi_b,T_N)}
   * \f]
   */
  double R_critical;

  /** Scaled version of the critical radius.
   *
   * Used to stop bubbles from collapsing when necessary.
   */
  double R_scaled;

  /** Value of temperature in the centre of the fluid sphere,
   * used in INIT_FLUID_SPHERE ("initfs") only.
   */
  double T_central;

  /** Radius of gaussian sphere of fluid,
   * used in INIT_FLUID_SPHERE ("initfs") only.
   */
  double sphere_radius;

  /** Constant in potential term.
   *
   * Calculated s.t \f$V(\phi_b,T=0)=0\f$
   * or s.t \f$ V_B(\phi_b)=0\f$ for BAG
   */
  double V0;
#ifdef BAG
  /** Broken phase value of phi in the bag model.
   */
  double phi_0;
#endif //BAG

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

  double *x_inv;
  double *phi_inv;
  double *V_inv;

  // Scalar field

  /** `phi` is the scalar field on the current timestep.
   */
  double ***phi;

  /** `pi_future` is the value of `pi` on the same timestep as `phi`.
   */
  double ***pi_future;

  /** `pi` is the conjugate momenta of the scalar field.
   */
  double ***pi;

  /** `phi_old` is the scalar field on the last timestep.
   */
  double ***phi_old;

#ifndef SCALAR
  // Fluid

  /** `T` is the temperature of the fluid.
   */
  double ***T;

  /** `E` is the internal energy density of the fluid.
   */
  double ***E;

  /** `W` is the zone-centred boost.
   */
  double ***W;

  /** `kappa` is the adiabatic index for the system.
   */
  double ***kappa;

  /** `p` is the pressure of the fluid. Obtained by eq_of_state() and
   * used in hydro.
   */
  double ***p;

  /** `V` is the fluid 3-velocity.
   */
  double ****V;

  /** `Z` is the fluid momentum density.
   */
  double ****Z;

  /** `U` is the fluid 4-velocity.
   */
  double ****U;

  /** `F` is a temporary variable used in advection.
   */
  double ****F;
#endif // SCALAR

  // Gravity

  /** `uij` are the unprojected metric peturbations.
   */
  double ****uij;

  /** `udotij` are the time derivative of the `uij` and are used for
   *    GW power spectrum.
   */
  double ****udotij;

} hydro_fields;

#ifdef FFT

/** Helper struct for performing fast fourier transforms.
 * Also contains fields at reference times for UETC's
 *
 * This documentation should be extended.
 */
typedef struct{

  /** Input field for fast fourier transforms.
   *
   */
  fftw_complex *in;

  /** Output field for fast fourier transforms.
   *
   */
  fftw_complex *out;

  /** Reference time field for Tij in momentum space.
   */
  fftw_complex **initial_Tij;

  /** Reference time field for velocity in momentum space.
   */
  fftw_complex **initial_V;

  /** Plan for performing the fast fourier transforms.
   *
   */
  fftw_plan plan;

  /** Map for which node is responsible for each `x` coordinate.
   */
  int *map;

} fft_fields;

#endif //FFT

/** Helper struct for performing MPI_MAXLOC and MPI_MINLOC operations.
 *
 * Useful for debugging, e.g finding max/min value of something on the
 * lattice and the processor containing the lattice site.
 */
typedef struct{
  double value;
  int rank;
} value_rank;

/** Helper struct for passing a value and array specifying its location.
 *
 * Useful for debugging, e.g finding max/min value of something on the
 * lattice and the lattice site where it occurs.
 */
typedef struct{
  double value;
  int loc[3];
} value_loc;

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


// comms.c
void layout(hydro_params *p);
int get_x(int n, hydro_params p);
int get_y(int n, hydro_params p);
int get_z(int n, hydro_params p);
void halo_field(double ***field, hydro_params p);
double reduce_sum(double result, hydro_params p);
int reduce_sum_int(int result, hydro_params p);
double reduce_max(double result, hydro_params p);
int reduce_max_int(int result, hydro_params p);
double reduce_min(double result, hydro_params p);
int reduce_and(int result, hydro_params p);
value_rank reduce_maxloc(double result, hydro_params p);
value_rank reduce_minloc(double result, hydro_params p);
void init_comms_time(hydro_params *p);
double get_comms_time(hydro_params *p);
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
double Vf(hydro_params p, double T, double this_phi);
double Vdf(hydro_params p, double T, double this_phi);
double VTf(hydro_params p, double T, double this_phi);
double VTTf(hydro_params p, double T, double this_phi);
void Vpot(hydro_params p, double ***T, double ***phi, double ***Vprecalc);
void Vdpot(hydro_params p, double ***T, double ***phi, double ***Vprecalc);


// energy.c
double field_energy(hydro_fields f, hydro_params p);
double gradient_energy_field(hydro_fields f, hydro_params p);
double kinetic_energy_field(hydro_fields f, hydro_params p);
double total_energy(hydro_fields f, hydro_params p);
double kinetic_energy_fluid(hydro_fields f, hydro_params p);
double rest_energy(hydro_fields f, hydro_params p);
void energy_density(hydro_fields f, hydro_params p, double ***en);
void stress_energy(hydro_fields f, hydro_params p, double ****Tij);
double avg_pressure(hydro_fields f, hydro_params p);
double tzerozero(hydro_fields f, hydro_params p);
double get_curlJ_tot(hydro_fields f, hydro_params p);
double get_divJ_tot(hydro_fields f, hydro_params p);

// eos.c

void find_Ta(hydro_fields f, hydro_params p);
void eq_of_state(hydro_fields f, hydro_params p);

// transport.c

void advect_E(hydro_fields f, hydro_params p, int adv_order);
void advect_Z(hydro_fields f, hydro_params p, int adv_order);
void donor_E_dir(hydro_fields f, hydro_params p, int dir);
void donor_Z_dir(hydro_fields f, hydro_params p, int dir);
void van_leer_E(hydro_fields f, hydro_params p, int dir);
void van_leer_Z(hydro_fields f, hydro_params p, int dir);



// initial.c
void initial_blank(hydro_fields f, hydro_params p);
int safe_distance(hydro_fields f, hydro_params p);
int can_nucleate(hydro_fields f, hydro_params p, int x0, int y0, int z0);
void nucleate_at(hydro_fields f, hydro_params p, int x0, int y0, int z0);
void initial_3D(hydro_fields f, hydro_params p);
int try_nucleate(hydro_fields f, hydro_params p);
int bubbles_at_step(hydro_fields f, hydro_params p, double t, int step);
void init_profile(hydro_fields *f, hydro_params *p);
void fluid_sphere(hydro_fields f, hydro_params p);


// output.c
void write_global_headers(hydro_fields f, hydro_params p);
void write_globals(hydro_fields f, hydro_params p, double gwen,
		    int bcount, double t_sim, int step);
double get_gamma_max(hydro_fields f, hydro_params p);
double get_s_max(hydro_fields f, hydro_params p);
double get_veltot(hydro_fields f, hydro_params p);
long long get_N_broken(hydro_fields f, hydro_params p);
long long get_broken_links(hydro_fields f, hydro_params p);
void dump(double *field, hydro_params p);
void histo_field(double ***field, hydro_params p, int step);
void didj(double *cpts, hydro_fields f, hydro_params p);
value_loc find_max_loc(double ***f, hydro_params p, int abs_max);
value_loc find_min_loc(double ***f, hydro_params p, int abs_min);
void dump_max_min(hydro_fields f, hydro_params p);

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

// silage_slice.c
void make_kinetic(hydro_fields f, hydro_params p, double ***temp);
void make_slice(hydro_fields f, hydro_params p, int xcoord, double *slice,
		double ***temp);
void make_vort(hydro_fields f, hydro_params p, double ****temp);
void make_divV(hydro_fields f, hydro_params p, double ***temp);
void make_curlJ(hydro_fields f, hydro_params p, double ****temp);
void make_divJ(hydro_fields f, hydro_params p, double ***temp);
void make_vel(hydro_fields f, hydro_params p, double ***temp);
void make_Z(hydro_fields f, hydro_params p, double ***temp);
void write_silo_slice_step(hydro_fields f, hydro_params p, int step,
			   int xcoord);
#endif // SILO


// util.c
double minof3(double a, double b, double c);
int minof3_int(int a, int b, int c);
double maxof3(double a, double b, double c);
double minof2(double a, double b);


#ifdef FFT
// fft_func.c
void fft_init(hydro_params p, fft_fields *fft_f);
void fft_finalise(hydro_params p, fft_fields *fft_f);
void fft_scalar(hydro_params p, fft_fields fft_f, double ***real_field);
void fft_vector(hydro_params p, fft_fields fft_f, double ****vector_field,
		fftw_complex **outcpts);
void fft_tensor(hydro_params p, fft_fields fft_f, double ****tensor_field,
		fftw_complex **outcpts, double norm);
void fft_J(hydro_fields f, hydro_params p, fft_fields fft_f, fftw_complex **outcpts);
void fft_X(hydro_fields f, hydro_params p, fft_fields fft_f, fftw_complex **outcpts);
void fft_e(hydro_fields f, hydro_params p, fft_fields fft_f);
double output_ps_uetcs(hydro_fields f, hydro_params p, fft_fields fft_f, int step);

// projectors.c
double proj(int T, double kx, double ky, double kz);
double lambda(int i, int j, int l, int m, double kx, double ky, double kz);
int indexof(int i, int j);
void tens_proj(hydro_params p, int x_start, int slab,
	       double *product, fftw_complex **tensor);
void uetc_tens_proj(hydro_params p, int x_start, int slab,
		    double *product_re, double *product_im,
		    fftw_complex **tensora, fftw_complex **tensorb);
void split_vector(hydro_params p, int x_start, int slab,
		  double *product_rot, double *product_div, double *product_tot,
		  fftw_complex **vec);
void uetc_split_vector(hydro_params p, int x_start, int slab,
		       double *product_rot_re, double *product_rot_im,
		       double *product_div_re, double *product_div_im,
		       fftw_complex **veca, fftw_complex **vecb);

// uetc.c
void init_uetc(hydro_fields f, hydro_params p, fft_fields fft_f);
void uetc_tensor(hydro_params p, fftw_complex **tensor_then,
		 fftw_complex **tensor_now, int step_then, int step_now,
		 char *label);
void uetc_vector(hydro_params p, fftw_complex **vector_then,
		 fftw_complex **vector_now, int step_then, int step_now,
		 char *label);

// tensorps.c
double tensorps(hydro_params p, fftw_complex **outcpts, int step, char *label);

// vectorps.c
void histogram(hydro_params p, double *slice, char *filename,
	       int slab, int x_start); // Why is this here?
void vectorps(hydro_params p, fftw_complex **outcpts, int step, char *label);

// scalarps.c
void scalarps(hydro_params p, fftw_complex *field, int step, char *label);

#endif // FFT


#if defined(FFT) && defined(BAG) && ! defined(SCALAR)
// initps.c
void init_ps(hydro_fields f, hydro_params p, double ****field);
void spectrum_interp(double ksq, hydro_params p, fftw_complex *res, double *k_bins, double *pow_bins, int n_bins);
void UtoZ(hydro_fields f, hydro_params p);
double get_normal(double mean, double dev);
void init_energy(hydro_params p, hydro_fields f, ptrdiff_t x_start, ptrdiff_t x_thickness, int* map, ptrdiff_t alloc_local,double *k_bins, double *pow_bins);
#endif // FFT && !SCALAR
