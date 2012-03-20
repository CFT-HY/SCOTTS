/* hydro.h
 *
 * Header file for C port of 1D spherical hydro code.
 */
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <strings.h>

typedef struct {
  double dx;
  double dt;
  int N;
  int steps;

  double Cav;
  double C;

  double Lheat;
  double sigma;
  double lcorr;

  // The following parameters are calculated
  double xxwall;
  double Tconst;

  double alpha;
  double lambda;
  double gamma;

  double a;

  double T0;
} hydro_params;

// evolve.c

void evolve_backstep(double *phi,
		     double *pifull, double *xe, double *xc, 
		      double *T, double *pi,
		      int **nb, hydro_params params);

void evolve_field(double *gb, double *v, double *xe, double *xc,
		  double *pi, double *phi,
		  double *T, double *phiold, double *pifull,
		  int **nb, hydro_params params);

void evolve_hydro(double *kappa,
		  double *phiold, double *phi, double *pi, double *p,
		  double *xe, double *xc,
		  double *E, double *Z, double *v, double *gb,
		  double *T,
                  int **nb, hydro_params params);


// potential.c

double Vf(hydro_params params, double T, double this_phi);

double Vdf(hydro_params params,
	   double T, double this_phi);

double VTf(hydro_params params,
	   double T, double this_phi);

double VTTf(hydro_params params,
	    double T, double this_phi);

void Vpot(hydro_params params,
          double *T,
          double *phi, double *Vprecalc);

void Vdpot(hydro_params params,
	   double *T,
	   double *phi, double *Vprecalc);


// energy,c

double total_energy(double *kappa, double *E,
		    double *gb, double *xe, double *xc,
                    double *phi, double *pifull,
		    double *T, int **nb, hydro_params params);


// eos.c

void find_Ta(double *E, double *gb, double *phi,
	     double *T, hydro_params params);

void eq_of_state(double *E, double *gb, double *phi,
		 double *T, double *p, double *kappa, hydro_params params);

// transport.c

void donor_r(double *v,
	     double *xe, double *xc, double *field, int **nb,
	     hydro_params params);

void donor_z(double *v,
	     double *xe, double *xc, double *field, int **nb,
	     hydro_params params);

void transport_r(double *v,
		 double *xe, double *xc, double *field, int **nb,
		 hydro_params params);

void transport_z(double *v,
                 double *xe, double *xc, double *field, int **nb,
		 hydro_params params);

// initial.c

void create_1D_bubble(double *xe, double *xc,
                      double *phi, double *pifull,
                      double *T, double *E, double *Z,
                      double *v, double *gb, 
		      hydro_params params);


double psibar(double x, double lbar);

void create_gaussian_bubble(double *xe, double *xc,
                            double *phi, double *pifull,
                            double *T, double *E, double *Z,
                            double *v, double *gb,
			    hydro_params params);

// output.c

double wallpos(double *xc, double *phi, double *T,
	       hydro_params params);


// parameters.c

hydro_params get_parameters();
