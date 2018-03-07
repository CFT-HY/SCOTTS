/** @file alloc.c
 *
 * Routines for memory allocation.
 *
 * Handles memory allocation, freeing and management, for
 * multi-dimensional arrays that are allocated contiguously in memory
 * for improved performance.
 *
 * 2010-2017 David Weir
 */
#include "hydro.h"


/** Allocate memory for a 3D, 1-component field
 *
 * Performs a sneaky trick to return 3D array with contiguous memory.
 * Accessed as: `field[x][y][z]`, with the x- and y-coordinates
 * having haloes at `x=0` and `x=p.slicex+1`, `y=0`, `y=p.slicey+1`.
 *
 * Use free_field() to free this memory again.
 */
float ***make_field(hydro_params p) {
   
  float *true_field = malloc((p.slicex+2)*(p.slicey+2)
			      *(p.Lz)*sizeof(float));


  float ***field = (float ***)malloc((p.slicex+2)*sizeof(float **));
  int x, y;

  for(x=0;x<(p.slicex+2);x++) {
    field[x] = (float **)malloc((p.slicey+2)*sizeof(float *));
    for(y=0;y<(p.slicey+2);y++) {
      field[x][y] = &true_field[x*(p.slicey+2)*(p.Lz) + y*(p.Lz)];
    }
  }

  field[0][0] = true_field;

  return field;
}



/** Allocate memory for a vector field.
 *
 * As for make_field() but for a field with three components.
 * See make_tensor() for a field with six components.
 *
 * Use free_vector() to free memory allocated with this function.
 */
float ****make_vector(hydro_params p) {
   
  float *true_field = malloc(3*(p.slicex+2)*(p.slicey+2)
			      *(p.Lz)*sizeof(float));

  int x, y, i;

  float ****vector = (float ****) malloc(3*sizeof(float***));

  for(i=0;i<3;i++) {
    vector[i] = (float ***)malloc((p.slicex+2)*sizeof(float **));

    for(x=0;x<(p.slicex+2);x++) {
      vector[i][x] = (float **)malloc((p.slicey+2)*sizeof(float *));
      for(y=0;y<(p.slicey+2);y++) {
	vector[i][x][y]
	  = &true_field[i*(p.slicex+2)*(p.slicey+2)*(p.Lz) 
			+ x*(p.slicey+2)*(p.Lz) + y*(p.Lz)];
      }
    }
  }

  vector[0][0][0] = true_field;
  
  
  return vector;
}



/** Allocate memory for a tensor.
 *
 * As for make_field() but for a field with six components (set by
 * `TENSOR_CPTS`). See make_vector() for a field with three
 * components.
 *
 * Use free_tensor() to free memory allocated with this function.
 */
float ****make_tensor(hydro_params p) {
   
  float *true_field = malloc(TENSOR_CPTS*(p.slicex+2)*(p.slicey+2)
			      *(p.Lz)*sizeof(float));

  int x, y, i;

  float ****tensor = (float ****) malloc(TENSOR_CPTS*sizeof(float***));

  for(i=0;i<TENSOR_CPTS;i++) {
    tensor[i] = (float ***)malloc((p.slicex+2)*sizeof(float **));

    for(x=0;x<(p.slicex+2);x++) {
      tensor[i][x] = (float **)malloc((p.slicey+2)*sizeof(float *));
      for(y=0;y<(p.slicey+2);y++) {
	tensor[i][x][y]
	  = &true_field[i*(p.slicex+2)*(p.slicey+2)*(p.Lz) 
			+ x*(p.slicey+2)*(p.Lz) + y*(p.Lz)];
      }
    }
  }

  tensor[0][0][0] = true_field;
  
  
  return tensor;
}



/** Free memory allocated by make_field().
 *
 * Frees the memory associated with a field allocated by make_field():
 * first frees the contiguous blob, then the 'shortcut' arrays, and
 * finally the outermost layer.
 */
void free_field(hydro_params p, float ***field) {

  int x;
  
  free(field[0][0]);

  for(x=0;x<(p.slicex+2);x++) {
      free(field[x]);
  }
  
  free(field);
}



 /** Free memory allocated by make_vector().
 *
 * Frees the memory associated with a field allocated by
 * make_vector(): first frees the contiguous blob, then the 'shortcut'
 * arrays and finally the outermost layer.
 */
void free_vector(hydro_params p, float ****vector) {

  int x, i;
  
  free(vector[0][0][0]);

  for(i=0;i<3;i++) {
    for(x=0;x<(p.slicex+2);x++) {
      free(vector[i][x]);
    }
    free(vector[i]);
  }

  free(vector);
}



/** Free memory allocated by make_tensor().
 *
 * Frees the memory associated with a field allocated by
 * make_tensor(): first frees the contiguous blob, then the 'shortcut'
 * arrays and finally the outermost layer.
 */
void free_tensor(hydro_params p, float ****tensor) {

  int x, i;
  
  free(tensor[0][0][0]);

  for(i=0;i<TENSOR_CPTS;i++) {
    for(x=0;x<(p.slicex+2);x++) {
      free(tensor[i][x]);
    }
    free(tensor[i]);
  }

  free(tensor);
}




/** Allocate all the fields needed for a simulation.
 */
void alloc_fields(hydro_fields *f, hydro_params p) {

  // (everything below also set up in initial condition functions)


  // SCALAR FIELD

  // scalar field and conjugate momentum
  f->phi = make_field(p);
  f->pifull = make_field(p);

  // pi gets initialised in backstep
  f->pi = make_field(p);

  // phiold initialised by evolve_field
  f->phiold = make_field(p);


#ifndef SCALAR
  // FLUID

  // temperature
  f->T = make_field(p);

  // fluid energy density
  f->E = make_field(p);

  // zone-centred boost
  f->W = make_field(p);
  
  // equation of state (obtained by eq_of_state first time
  // and used in hydro)
  f->kappa = make_field(p);

  // pressure (obtained by eq_of_state first time
  // and used in hydro)
  f->p = make_field(p);

  // 3-velocity
  f->V = make_vector(p); 

  // fluid momentum
  f->Z = make_vector(p);

  // 4-velocity
  f->U = make_vector(p);

  // F is a temporary variable used in advection
  f->F = make_vector(p);
#endif // SCALAR


  // GRAVITY

  // unprojected metric perturbations for GW power spectrum
  f->uij = make_tensor(p);
  f->udotij = make_tensor(p);

  // used for calculating unequal time correlators only
  f->initial_Tij = make_tensor(p);

}




/** Zero all the fields needed for a simulation.
 */
void zero_fields(hydro_fields f, hydro_params p) {

  int x, i;

  // We don't need to do this (because create_gaussian_bubble
  // initialises everything), but it keeps valgrind quiet, and
  // is on the safe side.
  for(x=0;x<p.fieldN;x++) {

    f.phi[0][0][x] = 0.0000;
    f.pifull[0][0][x] = 0.0000;
    f.pi[0][0][x] = 0.0000;
    f.phiold[0][0][x] = 0.0000;

#ifndef SCALAR
    f.T[0][0][x] = 0.0000;
    f.E[0][0][x] = 0.0000;
    f.W[0][0][x] = 0.0000;
    f.kappa[0][0][x] = 0.0000;
    f.p[0][0][x] = 0.0000;

    f.V[0][0][0][x] = 0.0000;
    f.V[1][0][0][x] = 0.0000;
    f.V[2][0][0][x] = 0.0000;

    f.Z[0][0][0][x] = 0.0000;
    f.Z[1][0][0][x] = 0.0000;
    f.Z[2][0][0][x] = 0.0000;

    f.U[0][0][0][x] = 0.0000;
    f.U[1][0][0][x] = 0.0000;
    f.U[2][0][0][x] = 0.0000;

    f.F[0][0][0][x] = 0.0000;
    f.F[1][0][0][x] = 0.0000;
    f.F[2][0][0][x] = 0.0000;
#endif // SCALAR


    for(i=0; i<TENSOR_CPTS; i++) {
      f.uij[i][0][0][x] = 0.0;
      f.udotij[i][0][0][x] = 0.0;

      f.initial_Tij[i][0][0][x] = 0.0;
    }
  }

}


/** Free all the fields needed for a simulation.
 */
void free_fields(hydro_fields *f, hydro_params p) {

  free_field(p, f->phi);
  free_field(p, f->pifull);
  free_field(p, f->pi);
  free_field(p, f->phiold);

#ifndef SCALAR
  free_field(p, f->T);
  free_field(p, f->E);
  free_field(p, f->W);
  free_field(p, f->kappa);
  free_field(p, f->p);
  free_vector(p, f->V);
  free_vector(p, f->Z);
  free_vector(p, f->U);
  free_vector(p, f->F);
#endif // SCALAR

  free_tensor(p, f->uij);
  free_tensor(p, f->udotij);
  free_tensor(p, f->initial_Tij);
 
}

