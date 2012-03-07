/* hydro.h
 *
 * Header file for C port of 1D spherical hydro code.
 */
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <strings.h>


// evolve.c

void evolve_backstep(double dt, double dx, int N, double *phi,
		     double *pifull, double *xe, double *xc, 
		     double alpha, double gamma, double lambda,
		     double *T, double T0, double *pi,
		     int **nb);

void evolve_field(double dt, double dx, double C, int N,
		  double *gb, double *v, double *xe, double *xc,
		  double *pi, double *phi,
		  double alpha, double gamma, double lambda,
		  double *T, double T0, double *phiold, double *pifull,
		  int **nb);

void evolve_hydro(double dt, double dx, double *kappa, double C, double Cav,
		  int N, double *phiold, double *phi, double *pi, double *p,
		  double *xe, double *xc,
		  double alpha, double gamma, double lambda,
		  double T0, double a, 
		  double *E, double *Z, double *v, double *gb,
		  double *T, double *Q,
                  int **nb);


// potential.c

double Vf(double alpha, double gamma, double lambda,
	  double T, double T0, double this_phi);

double Vdf(double alpha, double gamma, double lambda,
	   double T, double T0, double this_phi);

double VTf(double alpha, double gamma, double lambda,
	   double T, double T0, double this_phi);

double VTTf(double alpha, double gamma, double lambda,
	    double T, double T0, double this_phi);

void Vpot(int N, double alpha, double gamma,
	  double lambda, double *T, double T0,
	  double *phi, double *Vprecalc);

void Vdpot(int N, double alpha, double gamma,
	   double lambda, double *T, double T0,
	   double *phi, double *Vprecalc);


// energy,c

double total_energy(double dx, int N, double *kappa, double *E,
		    double *gb, double *xe, double *xc,
                    double *phi, double *pifull, double *Q,
		    double alpha, double gamma, double lambda,
                    double *T, int **nb);


// eos.c

void find_Ta(double *E, double *gb, double *phi, double a,
	     double alpha, double gamma, double lambda,
	     double T0, int N, double *T);

void eq_of_state(double a, int N, double *E, double *gb, double *phi,
		 double alpha, double gamma, double lambda,
		 double T0, double *T, double *p, double *kappa);

// transport.c

void donor_r(double dt, double dx, int N, double *v,
	     double *xe, double *xc, double *field, int **nb);

void donor_z(double dt, double dx, int N, double *v,
	       double *xe, double *xc, double *field, int **nb);

void transport_r(double dt, double dx, int N, double *v,
		 double *xe, double *xc, double *field, int **nb);

void transport_z(double dt, double dx, int N, double *v,
                 double *xe, double *xc, double *field, int **nb);

// initial.c

void create_1D_bubble(int N, double xxwall, double *xe, double *xc,
                      double *phi, double *pifull,
                      double *T, double *E, double *Z,
                      double *v, double *gb, double dx, double a,
                      double alpha, double gamma, double lambda,
		      double Tconst, double T0);

double psibar(double x, double lbar);

void create_gaussian_bubble(int N, double xxwall, double *xe, double *xc,
                            double *phi, double *pifull,
                            double *T, double *E, double *Z,
                            double *v, double *gb, double dx, double a,
                            double alpha, double gamma, double lambda,
			    double Tconst, double T0);


// output.c

double wallpos(double *xc, int N, double *phi, double *Q, double dt,
	       double alpha, double gamma, double lambda,
               double *T, double T0);


// parameters.c

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

hydro_params get_parameters();
